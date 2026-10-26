"""End-to-end screen harmonisation and validation (DESIGN v3 §7.1, §7.2, §7.6, §9).

    python scripts/harmonise_screens.py [--write-draft-config] [--skip-blind-write]

Blinding (DESIGN §9 and the G5a brief):
  * DepMap targets: DEV lines only, via ``essml.guard.load_targets`` (anchors are DepMap dev lines;
    anchor entries whose ModelID is not 'dev' are refused before any read).
  * Non-transformed (set-B) screens: control genes (CEGv2 / NEGv1 / scaling set) through
    ``guard.load_external_screen(..., controls_only=True)``; scaling parameters and QC use ONLY these.
    Full harmonised profiles are written to cache/harmonised/ (git-ignored) by ``blind_write``, which
    returns nothing and logs a 'harmonise_write_only' row in results/test_access_log.csv. Replicate
    ceilings read those files and return only aggregate pairwise statistics ('replicate_ceiling_
    aggregate' log rows). No model is trained or applied anywhere in this script.
  * Nothing gene-level from a non-transformed screen is printed or written to results/.

Outputs: configs/control_split.json, configs/screen_plan.csv, results/harmonisation/*.csv|json,
cache/harmonised/*.parquet, and (with --write-draft-config, only while the file is untracked)
configs/harmonisation_frozen.json.
"""
from __future__ import annotations

import argparse
import gzip
import io
import json
import re
import subprocess
import sys
import time
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from essml import genesets, guard  # noqa: E402
from essml import harmonise as H  # noqa: E402
from check_external_screens import ORCS_DIRECTION, load_hgnc, read_spec  # noqa: E402

CFG_SCREENS = ROOT / "configs" / "external_screens.json"
CFG_HARM = ROOT / "configs" / "harmonisation_screens.json"
CONTROL_SPLIT = ROOT / "configs" / "control_split.json"
FROZEN = ROOT / "configs" / "harmonisation_frozen.json"
PLAN = ROOT / "configs" / "screen_plan.csv"
OUT = ROOT / "results" / "harmonisation"
CACHE = ROOT / "cache" / "harmonised"
ORCS_DIR = ROOT / "data" / "external_screens" / "biogrid_orcs" / "extracted"
ORCS_META = ROOT / "results" / "external_screens" / "orcs_nontransformed_screens.csv"
SIM_SEED_OFFSET = 7  # simulated-batch RNG = split_seed + offset (documented)
ALT = (-1.0, 0.0)    # ALTERNATIVE affine targets (reported only); primary targets live in H.DEFAULTS (set in main)
ALT_NAME = "design_-1_0"
PRIMARY_NAME = "depmap_matched"
CHRONOS_DECISION = ROOT / "results" / "harmonisation" / "chronos_tkov3_decision.json"
CHRONOS_ADOPTED = bool(CHRONOS_DECISION.exists() and json.loads(CHRONOS_DECISION.read_text()).get("adopt_chronos_for_TKOv3"))


# ================================================================ helpers: IDs and tables
class SymbolMap:
    def __init__(self):
        norm, _, by_entrez = load_hgnc()
        self.norm = norm
        self.ent = {s: e for e, s in by_entrez.items()}

    def __call__(self, sym):
        if not isinstance(sym, str):
            return np.nan
        return self.ent.get(self.norm(sym), np.nan)


@lru_cache(maxsize=None)
def screens_cfg():
    return {s["id"]: s for s in json.load(open(CFG_SCREENS, encoding="utf-8"))["screens"]}


_TABLES = {}


def screen_table(screen_id, pools=False, sheet=None):
    key = (screen_id, pools, sheet)
    if key in _TABLES:
        return _TABLES[key]
    sc = screens_cfg()[screen_id]
    d = ROOT / "data" / "external_screens" / sc["deposit"]
    s = dict(sc["sanity"])
    if sheet is not None:
        s["sheet"] = sheet
        s["format"] = "excel"
    if pools:  # read pools separately so that RPM normalisation is within pool
        parts = []
        for sh in s["sheets"]:
            x = pd.read_excel(d / s["file"], sheet_name=sh, header=s.get("header", 0))
            x["__pool"] = sh
            parts.append(x)
        x = pd.concat(parts, ignore_index=True)
    else:
        x = read_spec(d, s)
    x.columns = [str(c).strip() if not str(c).startswith("Unnamed") else str(c) for c in x.columns]
    _TABLES[key] = (x, s, d / s["file"])
    return _TABLES[key]


def screen_path(entry):
    if entry["screen"].startswith("orcs:"):
        return ORCS_DIR / f"BIOGRID-ORCS-SCREEN_{entry['screen'][5:]}-2.0.18.screen.tab.txt"
    sc = screens_cfg()[entry["screen"]]
    return ROOT / "data" / "external_screens" / sc["deposit"] / sc["sanity"]["file"]


def guide_genes(x, s, symmap):
    if s["format"] == "zip_excel_mapped":
        g = x["__gene"]
    elif "gene_from_guide_regex" in s:
        rx = re.compile(s["gene_from_guide_regex"])
        g = x[s["guide_col"]].astype(str).map(lambda v: rx.search(v).group(1) if rx.search(v) else np.nan)
    else:
        g = x[s["gene_col"]]
    return g.map(symmap)


def raw_arm_scores(entry, symmap):
    """Raw gene score (Entrez-indexed Series) of one arm + an aggregate report (no gene values)."""
    lvl = entry["data_level"]
    rep = {"data_level": lvl}
    if entry["screen"].startswith("orcs:"):
        x = pd.read_csv(screen_path(entry), sep="\t", dtype={"IDENTIFIER_ID": str}, low_memory=False)
        v = pd.to_numeric(x["SCORE.1"], errors="coerce")
        s = v.groupby(x["IDENTIFIER_ID"].astype(str)).mean().dropna()
        rep.update(score="ORCS SCORE.1", n_genes=int(len(s)))
        return s, rep
    if lvl == "counts" and entry.get("chronos_screen") and CHRONOS_ADOPTED:
        f = CACHE / f"chronos_TKOv3__{entry['chronos_screen']}.parquet"
        sc = pd.read_parquet(f)["chronos_raw"].astype(float).dropna()
        rep.update(score="Chronos gene effect (joint TKOv3 fit)", t0_filter="reference reads >= 30 (late timepoints masked)",
                   n_genes=int(len(sc)))
        return sc, rep
    if lvl == "counts":
        x, s, _ = screen_table(entry["screen"], pools=entry.get("pools", False))
        if "arm" in entry:
            t0c, endc = s["arms"][entry["arm"]]
        else:
            t0c, endc = entry["t0"], entry["end"]
        genes = guide_genes(x, s, symmap)
        normalised = bool(entry.get("normalised", s.get("is_normalised", False)))
        lfc, r = H.guide_lfc(x[t0c], x[endc], normalised=normalised,
                             groups=x["__pool"] if entry.get("pools") else None)
        sc = H.gene_scores(lfc, genes)
        rep.update(r, score="median guide LFC", t0=" + ".join(t0c), end=" + ".join(endc),
                   n_guides_unmapped=int(genes.reindex(lfc.index).isna().sum()), n_genes=int(len(sc)))
        return sc, rep
    if lvl == "guide_lfc":
        x, s, _ = screen_table(entry["screen"])
        genes = guide_genes(x, s, symmap)
        v = x[entry["lfc_cols"]].apply(pd.to_numeric, errors="coerce").mean(axis=1)
        sc = H.gene_scores(v, genes)
        rep.update(score="median published guide LFC", n_guides=int(v.notna().sum()),
                   t0_filter="as published (no counts)", n_genes=int(len(sc)))
        return sc, rep
    if lvl == "gene_score":
        if "sheet" in entry:  # one sheet per arm (Tian 2019)
            sc_cfg = screens_cfg()[entry["screen"]]
            sp = sc_cfg["sanity"]
            y = pd.read_excel(ROOT / "data" / "external_screens" / sc_cfg["deposit"] / sp["file"],
                              sheet_name=entry["sheet"], header=sp.get("header", 0))
            g = y[sp["gene_col"]].map(symmap)
            v = pd.to_numeric(y[sp["score_col"]], errors="coerce")
        else:
            x, s, _ = screen_table(entry["screen"])
            g = x[s["gene_col"]].map(symmap)
            cols = entry.get("score_cols") or [s.get("score_col")]
            v = x[cols].apply(pd.to_numeric, errors="coerce").mean(axis=1)
        sc = v.groupby(g).mean().dropna()
        sc.index = sc.index.astype(str)
        rep.update(score="published gene score", n_genes=int(len(sc)))
        return sc, rep
    raise ValueError(lvl)


# ================================================================ controls and DepMap
def make_control_split(seed):
    ceg, _ = genesets.load_hart_ceg()
    neg, _ = genesets.load_hart_neg()
    gc = genesets.load_gene_classes()
    ce = set(gc.index[gc["common_essential"]])
    cs = H.control_split(ceg, neg, ce, seed)
    meta = {"_comment": "DESIGN v3 §7.2 control split (scripts/harmonise_screens.py). Fixed-seed random halves "
                        "(numpy default_rng(split_seed).permutation of sorted Entrez IDs) of Hart CEGv2 and NEGv1 "
                        "(Entrez via HGNC, essml.genesets). scaling_set = DepMap dev common essentials "
                        "(configs/gene_classes.csv: fraction of the 1022 dev lines < -0.5 above 0.95) minus the "
                        "CEGv2 evaluation half and minus NEGv1. The NEGv1 scaling half is mapped to 0, the scaling "
                        "set to -1; evaluation halves are held out for QC.",
            "seed": seed, "n": {k: len(v) for k, v in cs.items()}}
    out = {**meta, **cs}
    if CONTROL_SPLIT.exists():
        old = json.loads(CONTROL_SPLIT.read_text())
        if any(old.get(k) != v for k, v in cs.items()):
            raise RuntimeError("configs/control_split.json differs from the recomputed split")
    else:
        CONTROL_SPLIT.write_text(json.dumps(out, indent=1))
    return cs


def load_dev():
    sp = pd.read_csv(ROOT / "configs" / "splits.csv", index_col=0)
    dev = sp.index[sp["role"] == "dev"].tolist()
    Y = guard.load_targets(dev).astype(np.float64)
    return Y, Y.mean(axis=0, skipna=True)


# ================================================================ guarded set-B access
def guarded_controls(entry, symmap):
    """Control-gene raw scores of a non-transformed arm (guard controls_only) + aggregate report."""
    holder = {}

    def reader(_path):
        s, rep = raw_arm_scores(entry, symmap)
        if entry.get("rank_only"):
            s = H.rank_scores(s, entry["essential_is"])
        holder["rep"] = rep
        return s.to_frame(entry["id"]).T

    df = guard.load_external_screen(screen_path(entry), entry["id"], reader, controls_only=True)
    return df.iloc[0].dropna(), holder["rep"]


def blind_write(entry, symmap, params, corr):
    """Full harmonised profile -> cache/harmonised/<id>.parquet. Returns None (nothing is shown)."""
    s, _ = raw_arm_scores(entry, symmap)
    if entry.get("rank_only"):
        s = H.rank_scores(s, entry["essential_is"])
    guard._log("harmonise_write_only", pd.Series(["external"], index=[entry["id"]]), True,
               caller="harmonise_screens.py:blind_write",
               reason="full harmonised profile written to cache/harmonised (not returned/printed)",
               scope=f"setB:{entry['id']}", n_genes=int(len(s)))
    df = pd.DataFrame({"raw": s})
    if params is not None:
        df["affine"] = H.affine_apply(s, *params)
        df["final"] = H.apply_genewise(df["affine"], corr) if corr is not None else df["affine"]
        df["affine_alt"] = H.affine_apply(s, *params, *ALT)
    else:
        df["affine"] = df["final"] = df["affine_alt"] = s
    CACHE.mkdir(parents=True, exist_ok=True)
    df.index.name = "entrez"
    df.astype("float32").to_parquet(CACHE / f"{entry['id']}.parquet")
    return None


# ================================================================ anchor metrics
def agree(x, y, mdev, Ydev, own, ctrl_free):
    out = {"dev_corr": H.deviation_corr(x, y, mdev),
           "dev_corr_ctrlfree": H.deviation_corr(x, y, mdev, genes=ctrl_free),
           "raw_corr": float(pd.concat([x, y], axis=1).dropna().corr().iloc[0, 1])}
    out.update(H.nn_rank(x, Ydev, mdev, own))
    return out


def run_anchors(entries, symmap, cs, Ydev, mdev, ctrl_free, p53_ent, P):
    roles = guard.roles_of([e["depmap_id"] for e in entries])
    bad = roles[roles != "dev"]
    if len(bad):
        raise guard.BlindingError(f"anchor entries with non-dev ModelIDs: {bad.to_dict()}")
    prof, rows, reps = {}, [], []
    for e in entries:
        raw, rep = raw_arm_scores(e, symmap)
        y = Ydev.loc[e["depmap_id"]]
        ms, mn = H.affine_params(raw, cs["scaling_set"], cs["neg_scaling"])
        aff = H.affine_apply(raw, ms, mn)
        aff_dm = H.affine_apply(raw, ms, mn, *ALT)
        prof[e["id"]] = {"raw": raw, "affine": aff, "affine_dm": aff_dm}
        reps.append({"id": e["id"], "role": e["role"], "library": e["library"], "group": e["group"],
                     "depmap_id": e["depmap_id"], "m_scaling": ms, "m_neg": mn, **rep})
        for step, x in (("0_raw", raw), ("1_affine", aff), (f"1b_affine_alt_{ALT_NAME}", aff_dm)):
            q = H.qc_metrics(x, cs["ceg_evaluation"], cs["neg_evaluation"])
            rows.append({"id": e["id"], "role": e["role"], "group": e["group"], "library": e["library"],
                         "depmap_id": e["depmap_id"], "step": step, **q,
                         **agree(x, y, mdev, Ydev, e["depmap_id"], ctrl_free)})
    return prof, pd.DataFrame(rows), pd.DataFrame(reps)


def group_matrices(ids, prof, entries_by_id, Ydev, key="affine"):
    X = pd.DataFrame({i: prof[i][key] for i in ids}).T
    Y = Ydev.loc[[entries_by_id[i]["depmap_id"] for i in ids]].copy()
    Y.index = ids
    genes = X.columns.intersection(Y.columns)
    return X[genes], Y[genes]


def run_genewise(anchor_entries, prof, Ydev, mdev, ctrl_free, p53_ent, P, groups_cfg, seed, key="affine",
                 controls=True):
    by_id = {e["id"]: e for e in anchor_entries}
    loao_rows, kappa_rows, fits, loao_prof = [], [], {}, {}
    for grp in sorted({e["group"] for e in anchor_entries}):
        ids = [e["id"] for e in anchor_entries if e["group"] == grp]
        X, Y = group_matrices(ids, prof, by_id, Ydev, key)

        def metric(xc, a):
            return agree(xc, Y.loc[a], mdev, Ydev, by_id[a]["depmap_id"], ctrl_free)

        variants = [("location", {"scale": False})]
        if len(ids) >= P["min_anchors_scale"] + 1:  # LOAO fits then still have >= 8 anchors
            variants.append(("location_scale", {"scale": True}))
        for vname, kw in variants:
            res, pr, ft = H.loao(X, Y, metric, exclude=p53_ent, fit_kwargs={"kappa": P["kappa"], **kw})
            for a, c in ft.items():
                assert a not in c.anchors
            res["group"], res["variant"] = grp, vname
            res["fit_anchors"] = [";".join(ft[a].anchors) for a in res["anchor"]]
            loao_rows.append(res)
            loao_prof[(grp, vname)] = pr
        # controls on the step-(iii) fit (location variant):
        #  * within-set permutation of the DepMap references: a location offset depends only on the
        #    fit-set mean of Y, so this should reproduce the real gain (it is a gene-level library
        #    offset, not line matching) - reported as a mechanism check, NOT a negative control;
        #  * wrong-line references (negative control): each fit anchor is paired with a random
        #    non-anchor DepMap dev line; the offset then absorbs lineage biology and must not help.
        if controls and len(ids) >= 3:
            rng = np.random.default_rng(seed + sorted(groups_cfg).index(grp))
            pool = Ydev.index.difference([by_id[i]["depmap_id"] for i in ids])
            for vname in ("location_within_set_permutation", "location_wrong_line_reference"):
                rows = []
                for a in ids:
                    keep = [i for i in ids if i != a]
                    Yp = Y.loc[keep].copy()
                    if vname == "location_within_set_permutation":
                        Yp.iloc[:] = np.roll(Yp.to_numpy(), 1, axis=0)
                    else:
                        Yp.iloc[:] = Ydev.loc[rng.choice(pool, len(keep), replace=False), Yp.columns].to_numpy()
                    c = H.fit_genewise(X.loc[keep], Yp, exclude=p53_ent, kappa=P["kappa"])
                    rows.append({"anchor": a, "n_fit_anchors": len(keep), "estimable": c.estimable,
                                 "tau2": c.tau2, **metric(H.apply_genewise(X.loc[a], c), a)})
                r = pd.DataFrame(rows)
                r["group"], r["variant"], r["fit_anchors"] = grp, vname, ""
                loao_rows.append(r)
        # shrinkage sensitivity (informational; kappa is fixed a priori)
        def light(xc, a):
            return {"dev_corr": H.deviation_corr(xc, Y.loc[a], mdev)}

        for k in ((1.0, 2.0, 4.0) if controls else ()):
            res, _, _ = H.loao(X, Y, light, exclude=p53_ent, fit_kwargs={"kappa": k})
            kappa_rows.append({"group": grp, "kappa": k, "mean_dev_corr": res["dev_corr"].mean(),
                               "n": len(res)})
        for vname, kw in variants:
            fits[(grp, vname)] = H.fit_genewise(X, Y, exclude=p53_ent, kappa=P["kappa"], **kw)
    return pd.concat(loao_rows, ignore_index=True), pd.DataFrame(kappa_rows), fits, loao_prof


def decide(loao_df, step1, groups_cfg, anchor_entries, P):
    """Per-library decision: apply step (iii) only if LOAO improves deviation correlation."""
    margin = P["improvement_margin"]
    lib_of = {e["id"]: e["library"] for e in anchor_entries}
    base = step1.set_index("id")["dev_corr"]
    d = loao_df.copy()
    d["delta_dev_corr"] = d["dev_corr"] - d["anchor"].map(base)
    d["library"] = d["anchor"].map(lib_of)
    rows = []
    for grp, gcfg in groups_cfg.items():
        gd = d[d["group"] == grp]
        for lib in gcfg["libraries"]:
            row = {"group": grp, "library": lib, "n_group_anchors": int(gd["anchor"].nunique()),
                   "n_library_anchors": int(gd.loc[gd["library"] == lib, "anchor"].nunique())}
            if gd.empty or not gd["estimable"].any():
                rows.append({**row, "variant": "none", "applied": False,
                             "reason": "fewer than 2 anchors in the lab/library group: tau2 not identifiable, LOAO impossible"})
                continue
            best = None
            for v in ("location", "location_scale"):
                vd = gd[gd["variant"] == v]
                if vd.empty:
                    continue
                own = vd[vd["library"] == lib]
                if len(own):
                    gain = own["delta_dev_corr"].mean()
                    ok = gain > margin and (own["delta_dev_corr"] > 0).mean() >= 0.5
                    basis = f"own-library anchors ({len(own)})"
                else:
                    gain = vd["delta_dev_corr"].mean()
                    ok = gain > margin and (vd["delta_dev_corr"] > margin).all()
                    basis = f"no own-library anchor: all {len(vd)} group anchors must gain > margin"
                row[f"loao_gain_{v}"] = gain
                row[f"loao_ok_{v}"] = bool(ok)
                row[f"basis_{v}"] = basis
                if ok and (best is None or (v == "location_scale" and gain > row["loao_gain_location"] + margin)):
                    best = v
            for v, col in (("location_within_set_permutation", "within_set_perm_gain"),
                           ("location_wrong_line_reference", "wrong_line_ref_gain")):
                vv = gd[gd["variant"] == v]
                row[col] = vv["delta_dev_corr"].mean() if len(vv) else np.nan
            row["variant"] = best or "none"
            row["applied"] = best is not None
            row["reason"] = (f"LOAO gain > {margin}" if best else f"LOAO gain <= {margin} (or mixed): not applied")
            rows.append(row)
    return pd.DataFrame(rows), d


# ================================================================ simulated batch
def simulated_batch(Ydev, mdev, cs, p53_ent, P, seed):
    out = []
    for conv, tg in (("primary_" + PRIMARY_NAME, (H.DEFAULTS["affine_target_scaling"], H.DEFAULTS["affine_target_neg"])),
                     ("alt_" + ALT_NAME, ALT)):
        r = _simulated_batch(Ydev, mdev, cs, p53_ent, P, seed, tg)
        r.insert(0, "convention", conv)
        out.append(r)
    return pd.concat(out, ignore_index=True)


def _simulated_batch(Ydev, mdev, cs, p53_ent, P, seed, tg):
    rng = np.random.default_rng(seed)
    ids = Ydev.index
    genes = Ydev.columns[Ydev.notna().all(axis=0)]
    rows = []
    for scen, n_anchor, gw, scale in (("A3_location", 3, 0.25, False), ("A11_location_scale", 11, 0.25, True),
                                      ("A3_no_genewise_distortion", 3, 0.0, False),
                                      ("A11_no_genewise_distortion", 11, 0.0, True)):
        pick = rng.choice(len(ids), n_anchor + 20, replace=False)
        anc, tgt = ids[pick[:n_anchor]], ids[pick[n_anchor:]]
        Ya, Yt = Ydev.loc[anc, genes], Ydev.loc[tgt, genes]
        delta = pd.Series(rng.normal(0, gw, len(genes)), index=genes)
        lam = pd.Series(np.exp(rng.normal(0, 0.25, len(genes))) if (scale and gw > 0) else 1.0, index=genes)
        cen = Ya.mean(axis=0)
        a_true, b_true, noise = 3.0, 0.4, 0.15

        def distort(Yb):
            Z = cen + (Yb - cen) * lam + delta
            return a_true * Z + b_true + rng.normal(0, noise, Z.shape)

        Xa, Xt = distort(Ya), distort(Yt)

        def scale_rows(Xb):
            out, par = [], []
            for i in Xb.index:
                ms, mn = H.affine_params(Xb.loc[i], cs["scaling_set"], cs["neg_scaling"])
                out.append(H.affine_apply(Xb.loc[i], ms, mn, *tg))
                par.append((ms, mn))
            return pd.DataFrame(out, index=Xb.index), par

        Sa, _ = scale_rows(Xa)
        St, par_t = scale_rows(Xt)
        corr = H.fit_genewise(Sa, Ya, exclude=p53_ent, kappa=P["kappa"], scale=scale)
        Ct = H.apply_genewise(St, corr)
        nonp = genes.difference(p53_ent)
        for step, M in (("distorted_raw", Xt), ("affine", St), ("affine+genewise", Ct)):
            err = (M[nonp] - Yt[nonp]).to_numpy()
            dc = [H.deviation_corr(M.loc[i], Yt.loc[i], mdev) for i in tgt]
            q = [H.qc_metrics(M.loc[i], cs["ceg_evaluation"], cs["neg_evaluation"]) for i in tgt]
            rows.append({"scenario": scen, "n_anchors": n_anchor, "genewise_sd": gw, "scale_distortion": scale,
                         "step": step, "rmse_vs_truth": float(np.sqrt(np.nanmean(err ** 2))),
                         "mean_dev_corr_vs_truth": float(np.mean(dc)),
                         "median_ceg_eval": float(np.median([v["median_ceg_eval"] for v in q])),
                         "median_neg_eval": float(np.median([v["median_neg_eval"] for v in q])),
                         "corr_delta_hat_vs_planted": (float(np.corrcoef(corr.delta[nonp], delta[nonp])[0, 1])
                                                       if (step == "affine+genewise" and gw > 0) else np.nan),
                         "tau2_hat": corr.tau2 if step == "affine+genewise" else np.nan,
                         "p53_genes_untouched": bool(np.allclose(Ct[list(set(p53_ent) & set(genes))].to_numpy(),
                                                                St[list(set(p53_ent) & set(genes))].to_numpy()))
                         if step == "affine+genewise" else np.nan,
                         "affine_scale_recovered_median": float(np.median([mn - ms for ms, mn in par_t])) if step == "affine" else np.nan})
    return pd.DataFrame(rows)


# ================================================================ replicate ceilings
def replicate_ceilings(nt_entries, qc_pass, cs, mdev, gc, ctrl_free):
    """Pairwise agreement among QC-passing KO effect-scale screens of the same background."""
    frac = gc["frac_below"] if "frac_below" in gc else None
    rows = []
    for bg in ("RPE1", "hPSC"):
        es = [e for e in nt_entries if e.get("background") == bg and qc_pass.get(e["id"])
              and not e.get("rank_only") and e.get("modality", "KO") == "KO" and not e.get("focused")]
        prof = {}
        for e in es:
            p = CACHE / f"{e['id']}.parquet"
            prof[e["id"]] = pd.read_parquet(p)
        guard._log("replicate_ceiling_aggregate", pd.Series(["external"] * len(es), index=[e["id"] for e in es]),
                   True, caller="harmonise_screens.py:replicate_ceilings",
                   reason="pairwise aggregate agreement between observed screens (no model)", scope=f"setB:ceiling_{bg}")
        for i, ei in enumerate(es):
            for ej in es:
                if ei["id"] == ej["id"]:
                    continue
                for step in ("affine", "final", "affine_alt"):
                    xi, xj = prof[ei["id"]][step].astype(float), prof[ej["id"]][step].astype(float)
                    common = xi.dropna().index.intersection(xj.dropna().index).intersection(mdev.dropna().index)
                    cf = common.intersection(pd.Index(sorted(ctrl_free)))
                    lab = (xj[cf] < -0.5).to_numpy()
                    # Q3b-style: screen-specific = obs_j < -0.5 and dev dependency fraction <= 0.5
                    spec = cf[(xj[cf] < -0.5).to_numpy() & (frac.reindex(cf).fillna(1.0) <= 0.5).to_numpy()] if frac is not None else cf[:0]
                    negs = common.intersection(pd.Index(cs["neg_evaluation"]))
                    rows.append({
                        "background": bg, "predictor": ei["id"], "target": ej["id"], "step": step,
                        "same_study": ei["study"] == ej["study"],
                        "tp53_match": ei.get("tp53", "")[:2] == ej.get("tp53", "")[:2],
                        "n_genes": len(common), "n_ctrlfree": len(cf),
                        "dev_corr": H.deviation_corr(xi, xj, mdev, genes=common),
                        "dev_corr_ctrlfree": H.deviation_corr(xi, xj, mdev, genes=cf),
                        "dev_spearman_ctrlfree": float(pd.Series(xi[cf] - mdev[cf]).corr(pd.Series(xj[cf] - mdev[cf]), method="spearman")),
                        "auprc_ctrlfree_lt-0.5": H.auprc(-xi[cf].to_numpy(), lab),
                        "prevalence_ctrlfree_lt-0.5": float(lab.mean()),
                        "n_screen_specific": int(len(spec)),
                        "recall_at_5pct_fpr_screen_specific": H.recall_at_fpr(-xi, spec, negs) if len(spec) else np.nan})
    return pd.DataFrame(rows)


# ================================================================ main
def git_tracked(rel):
    return subprocess.run(["git", "ls-files", "--error-unmatch", rel], cwd=ROOT,
                          capture_output=True).returncode == 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write-draft-config", action="store_true")
    ap.add_argument("--skip-blind-write", action="store_true")
    a = ap.parse_args()
    t_start = time.time()

    def stage(msg):
        print(f"[{time.time() - t_start:6.0f}s] {msg}", flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    P = dict(H.DEFAULTS)
    frozen_targets = None
    if FROZEN.exists():
        fr = json.loads(FROZEN.read_text())
        frozen_targets = (fr.get("parameters", {}).get("affine_target_scaling"),
                          fr.get("parameters", {}).get("affine_target_neg"))
        frozen_targets = None if None in frozen_targets or fr.get("parameters", {}).get(
            "affine_target_mode") != P["affine_target_mode"] else frozen_targets
        for k, v in fr.get("parameters", {}).items():
            if k.startswith("affine_target_"):
                continue  # targets are recomputed from dev lines and verified below
            if k in P and P[k] != v:
                print(f"[config] {k}: default {P[k]} -> frozen {v}")
                P[k] = v
    H.DEFAULTS.update(P)
    seed = json.loads((ROOT / "configs" / "seeds.json").read_text())["split_seed"]
    hcfg = json.loads(CFG_HARM.read_text(encoding="utf-8"))
    groups_cfg = hcfg["groups"]
    symmap = SymbolMap()

    cs = make_control_split(seed)
    gc_full = pd.read_csv(ROOT / "configs" / "gene_classes.csv", dtype={"entrez": str}).set_index("entrez")
    ceg_all, _ = genesets.load_hart_ceg()
    neg_all, _ = genesets.load_hart_neg()
    p53_ent = sorted({symmap(s) for s in H.p53_ddr_symbols()} - {np.nan})
    print(f"control split: " + ", ".join(f"{k}={len(v)}" for k, v in cs.items()) + f"; p53/DDR genes={len(p53_ent)}")

    Ydev, mdev = load_dev()
    global ALT, ALT_NAME, PRIMARY_NAME
    DM = (float(Ydev[Ydev.columns.intersection(cs["scaling_set"])].median(axis=1).median()),
          float(Ydev[Ydev.columns.intersection(cs["neg_scaling"])].median(axis=1).median()))
    print(f"DepMap-matched affine targets (median over dev lines of per-line medians): scaling {DM[0]:.3f}, NEG half {DM[1]:.3f}")
    if P["affine_target_mode"] == "depmap_matched":   # PI decision D10
        primary, ALT, ALT_NAME, PRIMARY_NAME = DM, (-1.0, 0.0), "design_-1_0", "depmap_matched"
    else:
        primary, ALT, ALT_NAME, PRIMARY_NAME = (-1.0, 0.0), DM, "depmap_matched", "design_-1_0"
    if frozen_targets is not None and not np.allclose(frozen_targets, primary, atol=1e-6):
        raise RuntimeError(f"recomputed affine targets {primary} differ from frozen {frozen_targets}")
    P["affine_target_scaling"], P["affine_target_neg"] = primary
    H.DEFAULTS.update(P)
    print(f"primary affine targets ({PRIMARY_NAME}): {primary}; alternative reported ({ALT_NAME}): {ALT}")
    ctrl_free = set(Ydev.columns) - set(ceg_all) - set(neg_all)
    print(f"DepMap dev targets: {Ydev.shape}")

    # ---- DepMap dev QC distribution (Chronos) on the evaluation halves: threshold context
    dq = pd.DataFrame([H.qc_metrics(Ydev.loc[i], cs["ceg_evaluation"], cs["neg_evaluation"]) for i in Ydev.index],
                      index=Ydev.index)
    dq.describe(percentiles=[0.01, 0.05, 0.5]).T.to_csv(OUT / "depmap_dev_qc_summary.csv")

    entries = hcfg["entries"]
    anchors = [e for e in entries if e["role"] == "anchor"]
    durs = [e for e in entries if e["role"] == "duration_check"]

    stage("anchors steps 0/1")
    # ---- anchors: steps 0/1
    prof, anchor_rows, anchor_reps = run_anchors(anchors, symmap, cs, Ydev, mdev, ctrl_free, p53_ent, P)
    step1 = anchor_rows[anchor_rows["step"] == "1_affine"]

    stage("gene-wise LOAO")
    # ---- step (iii): LOAO per group, decisions
    loao_df, kappa_df, fits, loao_prof = run_genewise(anchors, prof, Ydev, mdev, ctrl_free, p53_ent, P, groups_cfg, seed)
    decisions, loao_detail = decide(loao_df, step1, groups_cfg, anchors, P)
    # alternative convention (DepMap-matched targets): reported, not used for the DESIGN decisions
    loao_dm, _, fits_dm, loao_prof_dm = run_genewise(anchors, prof, Ydev, mdev, ctrl_free, p53_ent, P, groups_cfg,
                                                     seed, key="affine_dm", controls=False)
    step1b = anchor_rows[anchor_rows["step"] == f"1b_affine_alt_{ALT_NAME}"]
    decisions_dm, loao_detail_dm = decide(loao_dm, step1b, groups_cfg, anchors, P)
    loao_detail_dm.to_csv(OUT / f"genewise_loao_alt_{ALT_NAME}.csv", index=False)
    decisions_dm.to_csv(OUT / f"genewise_decisions_alt_{ALT_NAME}.csv", index=False)
    loao_detail.to_csv(OUT / "genewise_loao.csv", index=False)
    kappa_df.to_csv(OUT / "genewise_kappa_sensitivity.csv", index=False)
    decisions.to_csv(OUT / "genewise_decisions.csv", index=False)
    # anchor agreement table incl. LOAO step(s)
    add = []
    qc_of = {}
    for (grp, v), pr in loao_prof.items():
        for aid, xc in pr.items():
            q = H.qc_metrics(xc, cs["ceg_evaluation"], cs["neg_evaluation"])
            r = loao_detail[(loao_detail["anchor"] == aid) & (loao_detail["variant"] == v)].iloc[0]
            add.append({"id": aid, "role": "anchor", "group": grp, "library": r["library"],
                        "depmap_id": next(e["depmap_id"] for e in anchors if e["id"] == aid),
                        "step": f"2_genewise_LOAO_{v}", **q,
                        **{k: r[k] for k in ("dev_corr", "dev_corr_ctrlfree", "raw_corr", "nn_rank", "nn_recall1",
                                             "nn_recall5", "nn_nauc", "n_ref_lines")}})
            CACHE.mkdir(parents=True, exist_ok=True)
            xc.to_frame("loao_corrected").astype("float32").to_parquet(CACHE / f"anchor_LOAO_{v}__{aid}.parquet")
    anchor_tab = pd.concat([anchor_rows, pd.DataFrame(add)], ignore_index=True)
    for g in groups_cfg:
        av = groups_cfg[g]["avana_like"]
        m = anchor_tab["group"] == g
        anchor_tab.loc[m, "qc_pass"] = anchor_tab.loc[m].apply(lambda r: H.qc_pass(r, avana_like=av), axis=1)
    anchor_tab.sort_values(["group", "id", "step"]).to_csv(OUT / "anchor_agreement.csv", index=False)
    anchor_reps.to_csv(OUT / "anchor_scoring_report.csv", index=False)

    # ---- correction to apply per library (final fit on all group anchors)
    corr_for_lib = {}
    for r in decisions.itertuples():
        if r.applied:
            corr_for_lib[r.library] = fits[(r.group, r.variant)]
    fit_summary = [{"group": g, "variant": v, "anchors": ";".join(c.anchors), "n_anchors": len(c.anchors),
                    "estimable": c.estimable, "tau2": c.tau2, "kappa": c.kappa, "tau2_scale": c.tau2_scale,
                    "sigma2": json.dumps({k: round(v_, 5) for k, v_ in c.sigma2.items()}),
                    "n_genes_corrected": int((c.delta != 0).sum()), "sd_delta_hat": float(c.delta.std()),
                    "n_p53ddr_excluded": len(c.excluded)} for (g, v), c in fits.items()]
    pd.DataFrame(fit_summary).to_csv(OUT / "genewise_fits.csv", index=False)

    stage("duration check")
    # ---- duration check on anchors (HAP1 T6..T18, HT-29 D7..D25 vs plasmid)
    drows = []
    for e in durs:
        raw, rep = raw_arm_scores(e, symmap)
        ms, mn = H.affine_params(raw, cs["scaling_set"], cs["neg_scaling"])
        aff = H.affine_apply(raw, ms, mn)
        y = Ydev.loc[e["depmap_id"]]
        q = H.qc_metrics(aff, cs["ceg_evaluation"], cs["neg_evaluation"])
        d = (aff - y).dropna()
        ce_dev = d.index.intersection(pd.Index(cs["scaling_set"]))
        drows.append({"id": e["id"], "library": e["library"], "days": e["days"], "raw_gap_neg_minus_scaling": mn - ms,
                      **{k: q[k] for k in ("median_ceg_eval", "median_neg_eval", "nnmd", "ssmd")},
                      **agree(aff, y, mdev, Ydev, e["depmap_id"], ctrl_free),
                      "median_affine_minus_depmap_all": float(d.median()),
                      "sd_affine_minus_depmap_all": float(d.std())})
    pd.DataFrame(drows).to_csv(OUT / "duration_check.csv", index=False)

    stage("set-B controls")
    # ---- non-transformed screens (set B): controls only through the guard
    nt = [dict(e) for e in entries if e["role"] == "nontransformed"]
    orcs_plan = orcs_entries(hcfg["orcs_rules"])
    nt += [e for e in orcs_plan if e["role"] == "nontransformed"]
    srows, qc_pass_map, params_of, per_screen_gw = [], {}, {}, []
    for e in nt:
        av = groups_cfg.get(e.get("group"), {}).get("avana_like", False)
        try:
            ctrl, rep = guarded_controls(e, symmap)
        except Exception as ex:  # report and continue
            srows.append({"id": e["id"], "step": "error", "note": repr(ex)[:200]})
            print("ERROR", e["id"], repr(ex)[:200])
            continue
        base = {"id": e["id"], "background": e.get("background"), "library": e.get("library"),
                "group": e.get("group"), "modality": e.get("modality", "KO"),
                "treatment": "rank-only" if e.get("rank_only") else "effect-scale",
                **{k: rep.get(k) for k in ("score", "n_genes", "n_guides", "n_guides_kept", "t0_filter",
                                           "pseudocount", "n_end_reps", "replicate_guide_lfc_r")}}
        if e.get("rank_only"):
            q = H.qc_metrics(ctrl, cs["ceg_evaluation"], cs["neg_evaluation"])
            ok = H.qc_pass(q, avana_like=False)
            srows.append({**base, "step": "rank", **q, "qc_pass": ok})
            qc_pass_map[e["id"]] = ok
            params_of[e["id"]] = None
            continue
        ms, mn = H.affine_params(ctrl, cs["scaling_set"], cs["neg_scaling"])
        params_of[e["id"]] = (ms, mn)
        aff = H.affine_apply(ctrl, ms, mn)
        steps = [("0_raw", ctrl), ("1_affine", aff), (f"1b_affine_alt_{ALT_NAME}", H.affine_apply(ctrl, ms, mn, *ALT))]
        corr = corr_for_lib.get(e.get("library")) if e.get("modality", "KO") == "KO" else None
        if e.get("genewise_override") is False:  # per-screen PI override (D12)
            corr = None
        if not e.get("rank_only"):
            per_screen_gw.append({"screen": e["id"], "library": e.get("library"), "modality": e.get("modality", "KO"),
                                  "gene_wise_applied": corr is not None,
                                  "reason": (e.get("genewise_override_reason") if e.get("genewise_override") is False
                                             else ("library decision: LOAO gain > margin" if corr is not None
                                                   else "library decision: not applied / no anchors / not KO")),
                                  "decision": "D12" if e.get("genewise_override") is False else "LOAO rule"})
        if corr is not None:
            steps.append(("2_genewise", H.apply_genewise(aff, corr)))
        for step, x in steps:
            q = H.qc_metrics(x, cs["ceg_evaluation"], cs["neg_evaluation"])
            ok = H.qc_pass(q, avana_like=av)
            srows.append({**base, "step": step, "m_scaling": ms, "m_neg": mn, **q, "qc_pass": ok})
        qc_pass_map[e["id"]] = H.qc_pass(H.qc_metrics(aff, cs["ceg_evaluation"], cs["neg_evaluation"]), avana_like=av)
        e["_corr"] = corr
    screen_qc = pd.DataFrame(srows)
    screen_qc.to_csv(OUT / "screen_qc.csv", index=False)

    stage("blind write")
    # ---- blind write of full harmonised profiles (cache only)
    if not a.skip_blind_write:
        for e in nt:
            if e["id"] in params_of:
                blind_write(e, symmap, params_of[e["id"]], e.get("_corr"))

    stage("replicate ceilings")
    # ---- replicate ceilings (aggregate only)
    ceil = replicate_ceilings(nt, qc_pass_map, cs, mdev, gc_full, ctrl_free)
    ceil.to_csv(OUT / "replicate_ceilings.csv", index=False)
    if len(ceil):
        summ = (ceil.groupby(["background", "step", "same_study", "tp53_match"])
                [["dev_corr", "dev_corr_ctrlfree", "dev_spearman_ctrlfree", "auprc_ctrlfree_lt-0.5",
                  "prevalence_ctrlfree_lt-0.5", "recall_at_5pct_fpr_screen_specific"]]
                .agg(["median", "min", "max", "count"]))
        summ.columns = ["__".join(c) for c in summ.columns]
        summ.reset_index().to_csv(OUT / "replicate_ceilings_summary.csv", index=False)
        # study-level: average both directions and all arm pairs within each unordered study pair first
        study_of = {e["id"]: e.get("study") for e in nt}
        cc = ceil[~ceil["same_study"]].copy()
        cc["pair"] = [" | ".join(sorted([study_of[i], study_of[j]])) for i, j in zip(cc["predictor"], cc["target"])]
        metr = ["dev_corr", "dev_corr_ctrlfree", "dev_spearman_ctrlfree", "auprc_ctrlfree_lt-0.5",
                "recall_at_5pct_fpr_screen_specific"]
        sp_ = cc.groupby(["background", "step", "pair", "tp53_match"])[metr].median().reset_index()
        sp_.to_csv(OUT / "replicate_ceilings_study_pairs.csv", index=False)
        sl = sp_.groupby(["background", "step"])[metr].agg(["median", "min", "max", "count"])
        sl.columns = ["__".join(c) for c in sl.columns]
        sl.reset_index().to_csv(OUT / "replicate_ceilings_study_level.csv", index=False)

    stage("simulated batch")
    # ---- simulated batch control
    sim = simulated_batch(Ydev, mdev, cs, p53_ent, P, seed + SIM_SEED_OFFSET)
    sim.to_csv(OUT / "simulated_batch.csv", index=False)

    # ---- screen plan + exclusions + draft config
    plan, excl = screen_plan(entries, orcs_plan, screen_qc, qc_pass_map, decisions, anchors, groups_cfg, ceil)
    plan.to_csv(PLAN, index=False)
    excl.to_csv(OUT / "exclusions.csv", index=False)
    draft = draft_config(P, decisions, cs, fits, seed)
    ov = json.loads(CHRONOS_DECISION.read_text()) if CHRONOS_DECISION.exists() else {"status": "Chronos not run"}
    ov["applied_in_this_run"] = CHRONOS_ADOPTED
    draft["scoring"]["TKOv3_override"] = ov
    draft["affine_targets"] = {
        "mode": PRIMARY_NAME, "scaling": P["affine_target_scaling"], "neg": P["affine_target_neg"],
        "provenance": "PI decision D10 (DESIGN '-1/0' read as 'map to DepMap scale'). Values = median over the 1022 "
                      "DepMap dev lines (configs/splits.csv role=dev, guard.load_targets) of each line's median "
                      "gene effect over configs/control_split.json scaling_set (target for the scaling set) and "
                      "neg_scaling (target for the NEGv1 scaling half). Recomputed and checked on every run.",
        "alternative_reported": {"mode": ALT_NAME, "targets": list(ALT),
                                 "genewise_decisions": decisions_dm.replace({np.nan: None}).to_dict("records")}}
    draft["genewise"]["per_screen"] = per_screen_gw
    (OUT / "harmonisation_choices.json").write_text(json.dumps(draft, indent=1))
    if a.write_draft_config:
        if git_tracked("configs/harmonisation_frozen.json"):
            print("configs/harmonisation_frozen.json is committed: not overwritten")
        else:
            FROZEN.write_text(json.dumps(draft, indent=1))
    runtime = time.time() - t_start
    (OUT / "run_info.json").write_text(json.dumps({"runtime_s": round(runtime, 1), "params": P,
                                                   "n_anchor_arms": len(anchors), "n_setB_arms": len(nt)}, indent=1))
    print_summary(anchor_tab, decisions, screen_qc, ceil, sim)
    print(f"runtime {runtime:.0f} s")


def orcs_entries(rules):
    meta = pd.read_csv(ORCS_META)
    meta = meta[meta["fitness_like"] == True]  # noqa: E712
    out = []
    for r in meta.to_dict("records"):
        sid = str(r["#SCREEN_ID"])
        cell = r["CELL_LINE"]
        meth = r["LIBRARY_METHODOLOGY"]
        e = {"id": f"ORCS_{sid}", "screen": f"orcs:{sid}", "data_level": "gene_score", "rank_only": True,
             "library": r["LIBRARY"], "group": "orcs_gene_level", "background": rules["background"].get(cell, cell),
             "study": f"orcs_{r['AUTHOR']}", "tp53": "KO" if "TP53" in str(r["CONDITION_NAME"]) else "WT/unknown",
             "duration": r["DURATION"], "readout": f"{r['PHENOTYPE']}; {r['EXPERIMENTAL_SETUP']}; condition: {r['CONDITION_NAME']}",
             "expression": rules["expression"].get(cell, "see inventory"), "score_type": r["SCORE.1_TYPE"],
             "modality": {"Knockout": "KO", "Inhibition": "CRISPRi", "Activation": "CRISPRa"}.get(meth, meth),
             "citation": r["AUTHOR"], "cell": cell}
        direction = ORCS_DIRECTION.get(r["SCORE.1_TYPE"], "unknown")
        e["essential_is"] = direction
        if sid in rules["duplicates_of_count_level"]:
            e.update(role="excluded", reason=f"duplicate of count-level screen {rules['duplicates_of_count_level'][sid]}")
        elif sid in rules["exclude"]:
            e.update(role="excluded", reason=rules["exclude"][sid])
        elif e["modality"] == "CRISPRa":
            e.update(role="excluded", reason=rules["activation_excluded"])
        elif direction not in ("low", "high"):
            e.update(role="excluded", reason=f"score direction unknown ({r['SCORE.1_TYPE']}): rank treatment impossible")
        elif not isinstance(r.get("extracted_file"), str):
            e.update(role="excluded", reason="ORCS file not extracted")
        else:
            e["role"] = "nontransformed"
        out.append(e)
    return out


def screen_plan(entries, orcs_plan, screen_qc, qc_pass_map, decisions, anchors, groups_cfg, ceil):
    lib_anchor = {}
    for e in anchors:
        lib_anchor.setdefault(e["group"], []).append(f"{e['id'][2:]}({e['depmap_id']})")
    applied = {r.library: r.variant for r in decisions.itertuples() if r.applied}
    sq = screen_qc[screen_qc["step"].isin(["1_affine", "rank"])].set_index("id")
    nmin = sq[["n_ceg_eval", "n_neg_eval"]].min(axis=1).to_dict() if len(sq) else {}
    # backgrounds with a cross-study replicate ceiling
    ceiling_bgs = set()
    if len(ceil):
        cc = ceil[(ceil["same_study"] == False) & (ceil["step"] == "final")]  # noqa: E712
        ceiling_bgs = set(cc["background"])
    rows, excl = [], []
    allE = [e for e in entries if e["role"] == "nontransformed"
            or (e["role"] == "excluded" and e.get("background"))] + orcs_plan
    seen = set()
    for e in allE:
        if e["id"] in seen:
            continue
        seen.add(e["id"])
        if e["role"] == "excluded":
            excl.append({"id": e["id"], "kind": "screen/arm", "reason": e.get("reason", "")})
        mod = e.get("modality", "KO")
        grp = e.get("group", "")
        anc = lib_anchor.get(grp, [])
        if grp == "orcs_gene_level" or e.get("same_pipeline_anchor") is False:
            anc = []  # rank-only ORCS, or anchors scored by a different pipeline (D12: Yilmaz vs Wang CS anchors)
        qc = qc_pass_map.get(e["id"])
        treat = ("excluded" if e["role"] == "excluded" else "descriptive (CRISPRi; QC only)" if mod == "CRISPRi"
                 else "rank-only (gene-level published score)" if e.get("rank_only")
                 else "effect-scale (median guide LFC + affine)" if e.get("data_level") in ("counts", "guide_lfc")
                 else "effect-scale (published gene score + affine)")
        expr = e.get("expression", "")
        no_expr = expr.startswith("none")
        has_ceiling = e.get("background") in ceiling_bgs
        if e["role"] == "excluded":
            q3a = q3b = f"none (excluded: {e.get('reason', '')})"
        elif mod == "CRISPRi":
            q3a = q3b = "descriptive only (CRISPRi stratum; no K562 CRISPRi anchor obtained)"
        elif qc is False and nmin.get(e["id"], 0) < H.DEFAULTS["qc_min_controls"]:
            q3a = q3b = "not evaluable (QC not assessable: < 10 evaluation-half CEG/NEG genes in the library)"
        elif qc is False:
            q3a = q3b = "not evaluable (QC fail)"
        elif no_expr:
            q3a = q3b = "not evaluable (no matched expression profile)"
        else:
            if not anc:
                q3a = "inconclusive only (no same-lab/library anchor: anchor-skill conditions of §7.6 not evaluable)"
            elif not has_ceiling:
                q3a = "transfers | inconclusive (no same-background replicate ceiling: MDE undefined, 'fails' unreachable)"
            else:
                q3a = "transfers | fails | inconclusive"
            if e.get("rank_only"):
                q3a = "rank label only; " + ("inconclusive only (no same-pipeline anchor)" if not anc else q3a)
                q3b = "not evaluable (rank-only: no effect-scale label < -0.5)"
            else:
                q3b = q3a
            if e.get("focused"):
                q3b = q3b + "; focused library (<= 2.5k genes): few screen-specific genes, wide CIs"
        rows.append({
            "screen": e["id"], "citation_or_source": e.get("citation", e.get("screen", "")),
            "background": e.get("background", ""), "category": e.get("category", e.get("background", "")),
            "modality_stratum": mod, "library": e.get("library", ""), "correction_group": grp,
            "data_level": e.get("data_level", ""), "treatment": treat, "readout": e.get("readout", ""),
            "tp53": e.get("tp53", ""), "duration": e.get("duration", ""),
            "qc_pass": "" if qc is None else bool(qc),
            "anchors_for_library_lab": ";".join(anc) if anc else ("none same-pipeline (D12: Wang anchors are published CS)" if e.get("same_pipeline_anchor") is False else "none"),
            "genewise_correction": ("none (PI D12 override)" if e.get("genewise_override") is False else
                                    applied.get(e.get("library"), "none")) if e["role"] != "excluded" and not e.get("rank_only") and mod == "KO" else "n/a",
            "matched_expression_source": expr, "replicate_ceiling_available": bool(has_ceiling) if e["role"] != "excluded" else False,
            "reachable_Q3a": q3a, "reachable_Q3b": q3b})
    for e in entries:
        if e["role"] == "excluded" and e["id"].startswith("X_") and not e.get("background"):
            excl.append({"id": e["id"], "kind": "anchor", "reason": e.get("reason", "")})
    return pd.DataFrame(rows), pd.DataFrame(excl)


def draft_config(P, decisions, cs, fits, seed):
    return {
        "_status": "DRAFT written by scripts/harmonise_screens.py --write-draft-config. Becomes binding only when "
                   "the PI commits it (DESIGN v3 §7.2, §9: tag freeze-transfer). Choices fixed on controls and "
                   "transformed anchors only.",
        "parameters": P,
        "scoring": {"counts": "RPM per sample (within pool for multi-pool libraries) + 0.5-read pseudocount; guides "
                              "with mean T0/plasmid reads < min_t0_reads dropped (RPM tables: < min_t0_rpm); LFC per end "
                              "replicate vs mean T0; mean over replicates; gene = median guide LFC over >= min_guides "
                              "guides; Entrez via HGNC (approved/previous/alias symbol).",
                    "guide_lfc_tables": "gene = median of the replicate-mean published guide LFC",
                    "gene_level_published": "affine-scaled published score (CS, CRISPRi phenotype)",
                    "orcs_bf_mageck": "rank-only (percentile, low = essential); rank metrics and AUPRC only",
                    "chronos_refit": "DESIGN §7.2.1 trigger met for TKOv3 (single anchor HAP1 fails NN recall under median LFC). "
                                     "Joint per-library Chronos fit by scripts/chronos_tkov3.py; adoption rule fixed in "
                                     "advance by the PI: HAP1 deviation-correlation gain >= +0.02 (see TKOv3_override)."},
        "control_split": {"file": "configs/control_split.json", "seed": seed,
                          **{f"n_{k}": len(v) for k, v in cs.items()}},
        "affine": "per screen arm: scaling-set median -> affine_targets.scaling, NEGv1 scaling-half median -> "
                  "affine_targets.neg (DepMap-matched, PI decision D10)",
        "genewise": {"model": "EB location (shrink to 0; tau2 from cross-anchor covariance; noise x kappa); scale "
                              "term only with >= 8 anchors; p53/DDR genes excluded from fit and application",
                     "decisions": decisions.replace({np.nan: None}).to_dict("records"),
                     "fits": {f"{g}|{v}": {"anchors": c.anchors, "tau2": c.tau2, "estimable": c.estimable}
                              for (g, v), c in fits.items()}},
        "qc": {"rule": "pass if SSMD(eval halves) < qc_ssmd_max for non-Avana libraries; NNMD < qc_nnmd_max for "
                       "Avana-type; rank-only screens: SSMD on percentile ranks",
               "nnmd_max": P["qc_nnmd_max"], "ssmd_max": P["qc_ssmd_max"]},
        "p53_ddr_genes": H.p53_ddr_symbols(),
    }


def print_summary(anchor_tab, decisions, screen_qc, ceil, sim):
    pd.set_option("display.width", 220)
    print("\n== anchors: deviation correlation with DepMap (same line) by step ==")
    t = anchor_tab.pivot_table(index=["group", "id"], columns="step", values="dev_corr")
    print(t.round(3).to_string())
    print("\n== anchors: NN rank among 1022 dev lines ==")
    print(anchor_tab.pivot_table(index=["group", "id"], columns="step", values="nn_rank").to_string())
    print("\n== gene-wise decisions ==")
    cols = [c for c in decisions.columns if c.startswith(("group", "library", "n_", "loao_gain", "within", "wrong", "variant", "applied"))]
    print(decisions[cols].round(4).to_string())
    print("\n== set-B screen QC (affine step; controls only) ==")
    s = screen_qc[screen_qc["step"].isin(["1_affine", "rank"])]
    print(s[["id", "treatment", "n_ceg_eval", "n_neg_eval", "median_ceg_eval", "median_neg_eval", "nnmd", "ssmd",
             "auroc", "qc_pass"]].round(3).to_string())
    if len(ceil):
        print("\n== replicate ceilings (cross-study, final step) ==")
        c = ceil[(~ceil["same_study"]) & (ceil["step"] == "final")]
        print(c.groupby(["background", "tp53_match"])[["dev_corr", "dev_corr_ctrlfree", "auprc_ctrlfree_lt-0.5",
                                                       "recall_at_5pct_fpr_screen_specific"]].median().round(3).to_string())
    print("\n== simulated batch ==")
    print(sim[["convention", "scenario", "step", "rmse_vs_truth", "mean_dev_corr_vs_truth", "median_ceg_eval", "median_neg_eval",
               "corr_delta_hat_vs_planted"]].round(3).to_string())


if __name__ == "__main__":
    main()
