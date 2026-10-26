"""G3 kernel-track report (D16, D17): results table with within-library / within-lineage P1 and
Nadeau-Bengio SE, controls vs expected outcomes, M1b tuned weights, M2, secondary splits.

    python scripts/g3_report.py        # -> docs/g3_kernel_results.md, results/g3_kernel_summary.json

Uses the latest registry row with status ok for each run name. Headline metrics (P1, P2b, runtime)
come from the frozen evaluation outputs (results/runs/<id>/per_fold.csv, summary.json). The added
columns are computed here from the cached out-of-fold predictions (cache/runs/<id>/pred.parquet)
with the frozen ``essml.metrics.per_gene_pearson``, dev targets read through ``essml.guard``:

* **P1 within library group** -- reviews/g2_review/recompute_p1.py (``group_centre``, lines
  68-81; groups lines 94-98): within each outer fold, subtract the per-gene mean over the observed
  held-out lines of each library group (configs/library_groups.csv: the dev row NaN-pattern groups,
  4 largest + other) from BOTH truth and prediction, then per-gene Pearson r. Same genes as P1
  (eval genes, dev-wide variable class, >= 40 observed held-out lines), mean over genes, then folds.
* **P1 within lineage** -- the lineage-partialled P1 of reviews/g2_review/neg_genes.py
  (``partial_out`` with lineage dummies, lines 31-39 and 84-89) applied to variable genes. With
  lineage dummies only, the OLS residual on [1, dummies] over a gene's observed rows equals
  subtracting the per-lineage observed means, i.e. ``group_centre`` by lineage
  (recompute_p1.py "within_lineage"); that is what is computed.
* **Nadeau-Bengio corrected SE** (review finding 12) for a paired fold-wise difference d_k over K
  folds: SE_NB = sqrt((1/K + n_test/n_train) * s_d^2), n_test/n_train averaged over folds
  (Nadeau & Bengio 2003, Mach. Learn. 52:239). Disclosure only: the selection rule (DESIGN §4:
  dP1 > 0.005 and > 1 fold SE) is unchanged.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

for _v in ["OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"]:
    os.environ.setdefault(_v, "4")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
from g2_report import fold_vals, holm, latest_ok, load_run, paired, subset_eval  # noqa: E402

DELTA = 0.005
HEADLINE = "m1_E_resB0"
MIN_OBS = 40
REF = ["b1", "b2", "r1_E", "m1_E", "t_E_status", HEADLINE]
CONTROLS = ["ctl_resB0_permuted", "ctl_resB0_perm_lineage", "ctl_resB0_perm_linlib", "m1_E_permuted"]
MD = ROOT / "docs" / "g3_kernel_results.md"
OUT_JSON = ROOT / "results" / "g3_kernel_summary.json"
EC_DEV = ["ACH-001578", "ACH-002291", "ACH-001374", "ACH-001668"]


# ============================================================================ statistics
def nb_se(d, ratio):
    """Nadeau-Bengio corrected SE of the mean of K fold-wise differences."""
    d = np.asarray(d, np.float64)
    return float(np.sqrt((1.0 / len(d) + ratio) * d.var(ddof=1)))


def paired_nb(a, b, ratio):
    """g2_report.paired plus the Nadeau-Bengio SE and its t-test (df K-1; disclosure only)."""
    out = paired(a, b)
    se = nb_se(np.asarray(a) - np.asarray(b), ratio)
    out["se_nb"] = se
    k = len(a)
    out["p_nb"] = float(2 * stats.t.sf(abs(out["mean"] / se), k - 1)) if se > 0 else (
        0.0 if out["mean"] != 0 else 1.0)
    return out


def test_train_ratio(folds: pd.Series, excluded=()):
    r = []
    for f in sorted(folds.unique()):
        n_te = int((folds == f).sum())
        n_tr = int(((folds != f) & ~folds.index.isin(list(excluded))).sum())
        r.append(n_te / n_tr)
    return float(np.mean(r))


# ============================================================================ within-group P1
def group_centre(P, T, labels):
    """reviews/g2_review/recompute_p1.py:68-81 (verbatim logic): within each label group subtract
    the per-gene mean over observed (truth non-NaN) rows from truth and prediction."""
    P, T = P.copy(), T.copy()
    for g in np.unique(labels):
        m = labels == g
        obs = ~np.isnan(T[m])
        n = obs.sum(0)
        with np.errstate(invalid="ignore", divide="ignore"):
            mt = np.where(obs, T[m], 0).sum(0) / n
            mp = np.where(obs, P[m], 0).sum(0) / n
        T[m] -= np.nan_to_num(mt)
        P[m] -= np.nan_to_num(mp)
    return P, T


def within_p1(pred: pd.DataFrame, Y: pd.DataFrame, folds: pd.Series, groupings: dict,
              cls_mask: np.ndarray, min_obs=MIN_OBS) -> dict:
    """Per-fold P1 on the class genes: 'raw' and centred within each grouping (name -> Series)."""
    from essml import metrics
    out = {k: [] for k in ["raw", *groupings]}
    genes = Y.columns
    for f in sorted(folds.unique()):
        ids = folds.index[folds == f]
        T = Y.loc[ids].to_numpy(np.float64)
        P = pred.loc[ids, genes].to_numpy(np.float64)
        ev = cls_mask & ((~np.isnan(T)).sum(0) >= min_obs)
        out["raw"].append(float(np.nanmean(metrics.per_gene_pearson(P, T, min_obs)[ev])))
        for name, lab in groupings.items():
            Pc, Tc = group_centre(P, T, lab.loc[ids].to_numpy())
            out[name].append(float(np.nanmean(metrics.per_gene_pearson(Pc, Tc, min_obs)[ev])))
    return out


# ============================================================================ data
def dev_frame():
    sp = pd.read_csv(ROOT / "configs" / "splits.csv", index_col=0)
    dev = sp[sp["role"] == "dev"]
    lib = pd.read_csv(ROOT / "configs" / "library_groups.csv", index_col=0)["lib_group"]
    return dev, lib.reindex(dev.index)


def load_targets_eval(dev):
    from essml import guard
    genes = pd.Index(pd.read_csv(ROOT / "configs" / "eval_genes.csv", dtype={"entrez": str})["entrez"])
    Y = guard.load_targets(dev.index).astype(np.float64)[genes]
    return Y, genes


def h3_genes():
    d = json.loads((ROOT / "configs" / "transfer_criteria.json").read_text())
    return {k: v for k, v in d["hypotheses"]["H3"]["genes"].items() if v}


# ============================================================================ main
def main(argv=None):
    reg = pd.read_csv(ROOT / "results" / "registry.csv", dtype=str)
    names = list(dict.fromkeys(REF + CONTROLS + sorted(
        n for n in reg["name"].dropna().unique()
        if str(n).startswith(("m1bt_", "m2n_")) or str(n).endswith(("_noEC", "_lolo")))))
    ids = {n: latest_ok(reg, n) for n in names}
    missing = [n for n, r in ids.items() if r is None]
    ids = {n: r for n, r in ids.items() if r is not None}
    runs = {n: load_run(r) for n, r in ids.items()}
    dev, lib = dev_frame()
    folds_cv = dev["cv_fold"].astype(int)
    lineage = dev["OncotreeLineage"].astype(str)
    Y, genes = load_targets_eval(dev)
    from essml import genesets
    gc = genesets.load_gene_classes().reindex(genes).fillna(False).astype(bool)
    var_mask = gc["variable"].to_numpy()
    groupings = {"within_library": lib, "within_lineage": lineage}

    cv_runs = [n for n in runs if runs[n][0].get("fold_scheme", "cv") == "cv"
               and not runs[n][0]["config"].get("line_subset")]
    head_pf = runs[HEADLINE][1] if HEADLINE in runs else None
    b1_pf = runs["b1"][1] if "b1" in runs else None
    rows, within = [], {}
    for n in cv_runs:
        s, pf = runs[n]
        sm = s["summary"]
        excl = s["config"].get("exclude_from_training", [])
        ratio = test_train_ratio(folds_cv, excl)
        p1v = fold_vals(pf, "variable", "P1")
        pred_path = ROOT / "cache" / "runs" / ids[n] / "pred.parquet"
        w = None
        if pred_path.exists():
            w = within_p1(pd.read_parquet(pred_path), Y, folds_cv, groupings, var_mask)
            if not np.allclose(w["raw"], p1v, atol=1e-6):  # same definition as the frozen P1
                raise AssertionError(f"{n}: recomputed P1 {w['raw']} != per_fold {p1v.tolist()}")
            within[n] = w
        row = {"name": n, "run_id": ids[n], "model": s["config"]["model"],
               "P1_var": sm["variable"]["P1"], "P1_var_sd": sm["variable"]["P1_sd"],
               "P1_all": sm["all"]["P1"], "P2b_nap": sm["selective"]["P2b_nap_prev05"],
               "P2a": sm["all"]["P2a"], "fit_s": s["fit_s"], "predict_s": s["predict_s"],
               "threads": s["threads"], "tree_clean": s.get("tree_clean"), "ratio_nb": ratio}
        if w is not None:
            for k in groupings:
                row[k] = float(np.mean(w[k]))
                row[k + "_sd"] = float(np.std(w[k], ddof=1))
        if head_pf is not None and n != HEADLINE:
            d = paired_nb(p1v, fold_vals(head_pf, "variable", "P1"), ratio)
            row.update({"dP1_head": d["mean"], "dP1_head_se": d["se"], "dP1_head_se_nb": d["se_nb"],
                        "dP1_head_pos": d["n_folds_pos"], "dP1_head_p": d["p"]})
        if b1_pf is not None and n != "b1":
            d = paired_nb(p1v, fold_vals(b1_pf, "variable", "P1"), ratio)
            row.update({"dP1_B1": d["mean"], "dP1_B1_se": d["se"], "dP1_B1_se_nb": d["se_nb"]})
        rep = s.get("model_report") or {}
        if rep:
            is_m2 = s["config"]["model"] == "M2"
            m1rep = {f: (v.get("m1", {}) if is_m2 else v) for f, v in rep.items()}
            r0 = next(iter(m1rep.values()))
            row["lam"] = json.dumps({k: sorted({float(v.get("lam_by_group", {}).get(k, np.nan))
                                                for v in m1rep.values()})
                                     for k in r0.get("lam_by_group", {})})
            if "weights" in r0:
                row["weights"] = {f: v["weights"] for f, v in rep.items()}
                row["n_eigh"] = {f: v.get("tune_n_eval") for f, v in rep.items()}
            if is_m2:
                row["lam_local"] = {f: v.get("lam_by_group") for f, v in rep.items()}
        rows.append(row)
    tab = {r["name"]: r for r in rows}

    # within-group paired deltas vs headline and B1 (disclosure)
    wdel = {}
    for n, w in within.items():
        for ref in [HEADLINE, "b1"]:
            if ref in within and n != ref:
                for k in groupings:
                    d = paired_nb(w[k], within[ref][k], tab[n]["ratio_nb"])
                    wdel.setdefault(n, {})[f"{k}_vs_{ref}"] = d

    # ---- M1b family (declared contrasts vs headline; Holm)
    m1b = sorted(n for n in tab if n.startswith("m1bt_"))
    m1b_con = []
    for n in m1b:
        d = paired_nb(fold_vals(runs[n][1], "variable", "P1"), fold_vals(head_pf, "variable", "P1"),
                      tab[n]["ratio_nb"])
        m1b_con.append({"contrast": f"{n} vs {HEADLINE}", "A": n, "B": HEADLINE, **d})
    if m1b_con:
        for c, pa in zip(m1b_con, holm([c["p"] for c in m1b_con])):
            c["p_holm"] = float(pa)
            c["helps"] = bool(c["mean"] > 0 and pa < 0.05)
    # C common-subset rule (lines with CN data)
    c_runs = [n for n in [HEADLINE, "b1"] + m1b if n in ids]
    csub = subset_eval({n: ids[n] for n in c_runs}, "C") if any(n.startswith("m1bt_") for n in c_runs) else {}
    c_con = []
    for n in m1b:
        if n in csub and HEADLINE in csub:
            d = paired(csub[n]["P1_variable_folds"], csub[HEADLINE]["P1_variable_folds"])
            c_con.append({"contrast": f"{n} vs {HEADLINE} (C-covered lines only)", "A": n, "B": HEADLINE,
                          **d})

    # ---- M2 family
    m2 = sorted(n for n in tab if n.startswith("m2n_"))
    m2_prior = []
    for a, b in [("paralogs", "none"), ("partners", "paralogs")]:
        pa, pb = f"m2n_E_{a}", f"m2n_E_{b}"
        if pa in runs and pb in runs:
            d = paired_nb(fold_vals(runs[pa][1], "variable", "P1"), fold_vals(runs[pb][1], "variable", "P1"),
                          tab[pa]["ratio_nb"])
            m2_prior.append({"contrast": f"M2 {a} vs {b}", "A": pa, "B": pb, **d})
    if m2_prior:
        for c, pa in zip(m2_prior, holm([c["p"] for c in m2_prior])):
            c["p_holm"] = float(pa)

    # ---- selection rule vs the headline (complexity order M1b < M2)
    sel = []
    for n in m1b + m2:
        r = tab[n]
        win = r["dP1_head"] > DELTA and r["dP1_head"] > r["dP1_head_se"]
        sel.append({"candidate": n, "vs": HEADLINE, "dP1": r["dP1_head"], "se": r["dP1_head_se"],
                    "se_nb": r["dP1_head_se_nb"], "replaces": bool(win)})
    winners = [s for s in sel if s["replaces"]]
    selected = HEADLINE
    for s in winners:  # complexity order: first M1b winner, then M2 must beat it by the same rule
        if selected == HEADLINE:
            selected = s["candidate"]
        else:
            d = paired(fold_vals(runs[s["candidate"]][1], "variable", "P1"),
                       fold_vals(runs[selected][1], "variable", "P1"))
            if d["mean"] > DELTA and d["mean"] > d["se"]:
                selected = s["candidate"]

    # ---- controls vs expected outcomes
    ctl = control_checks(tab, within)

    # ---- secondary splits
    lolo = lolo_section(runs, ids)
    ec = ec_section(runs, ids, Y, dev)

    out = {"runs": rows, "missing_runs": missing, "within_deltas": wdel, "m1b_contrasts": m1b_con,
           "m1b_c_subset": {"runs": csub, "contrasts": c_con}, "m2_prior_contrasts": m2_prior,
           "selection_log": sel, "selected_model": selected, "controls": ctl, "lolo": lolo,
           "ec": ec, "delta": DELTA}
    OUT_JSON.write_text(json.dumps(out, indent=1, default=float) + "\n")
    write_md(out, tab)
    print(json.dumps({"selected": selected, "missing": missing,
                      "controls": {c["run"]: c["pass"] for c in ctl}}, indent=1))


def control_checks(tab, within):
    """Numerical checks of the expected outcomes stated in the control configs (make_g3_configs)."""
    out = []
    b1 = tab.get("b1", {}).get("P1_var")
    head = tab.get(HEADLINE, {}).get("P1_var")
    b1_lib = tab.get("b1", {}).get("within_library")
    for n in CONTROLS:
        if n not in tab:
            continue
        r = tab[n]
        p1, wl, wb = r["P1_var"], r.get("within_lineage"), r.get("within_library")
        checks = {}
        if n == "ctl_resB0_permuted":
            checks["abs(P1) < 0.01"] = abs(p1) < 0.01
        elif n == "ctl_resB0_perm_lineage":
            checks["B1-0.04 <= P1 <= B1+0.01"] = b1 - 0.04 <= p1 <= b1 + 0.01
            checks["abs(within-lineage P1) < 0.01"] = wl is not None and abs(wl) < 0.01
            checks["headline - P1 > 0.08"] = head - p1 > 0.08
        elif n == "ctl_resB0_perm_linlib":
            checks["B1-0.04 <= P1 <= B1+0.03"] = b1 - 0.04 <= p1 <= b1 + 0.03
            checks["within-lineage P1 < 0.03"] = wl is not None and wl < 0.03
            checks["within-library P1 <= B1 within-library + 0.01"] = (
                wb is not None and b1_lib is not None and wb <= b1_lib + 0.01)
            checks["headline - P1 > 0.08"] = head - p1 > 0.08
        elif n == "m1_E_permuted":  # G2 control on the B1 path (reference)
            checks["abs(P1 - B1) < 0.01"] = abs(p1 - b1) < 0.01
        out.append({"run": n, "P1_var": p1, "within_lineage": wl, "within_library": wb,
                    "dP1_B1": r.get("dP1_B1"), "checks": checks, "pass": all(checks.values()),
                    "lam": r.get("lam")})
    return out


def lolo_section(runs, ids):
    """Per-lineage P1 (variable) and P2b (selective, prevalence <= 0.5) of leave-one-lineage-out runs."""
    lolo = {n: runs[n] for n in runs if n.endswith("_lolo")}
    if not lolo:
        return {}
    per = {}
    for n, (s, pf) in lolo.items():
        names = {int(k): v for k, v in (s.get("fold_names") or {}).items()}
        v = pf[pf["class"] == "variable"].set_index("fold")
        se = pf[pf["class"] == "selective"].set_index("fold")
        for f, lin in names.items():
            per.setdefault(lin, {"n_lines": int(v.loc[f, "n_lines"]) if f in v.index else None})
            per[lin][n] = {"P1": float(v.loc[f, "P1"]) if f in v.index else np.nan,
                           "n_genes": int(v.loc[f, "n_genes"]) if f in v.index else 0,
                           "P2b": float(se.loc[f, "P2b_nap_prev05"]) if f in se.index else np.nan}
    summ = {}
    for n, (s, pf) in lolo.items():
        p = np.array([per[l][n]["P1"] for l in per if n in per[l]], float)
        q = np.array([per[l][n]["P2b"] for l in per if n in per[l]], float)
        summ[n] = {"run_id": ids[n], "P1_mean_scored": float(np.nanmean(p)),
                   "n_scored": int(np.isfinite(p).sum()), "P2b_mean_scored": float(np.nanmean(q)),
                   "fit_s": s["fit_s"]}
    # paired over scored lineages (descriptive): model vs B2
    pairs = {}
    for a in lolo:
        for b in lolo:
            if a != b and b.startswith("b2"):
                d = [per[l][a]["P1"] - per[l][b]["P1"] for l in per
                     if a in per[l] and b in per[l] and np.isfinite(per[l][a]["P1"])
                     and np.isfinite(per[l][b]["P1"])]
                if d:
                    pairs[f"{a} - {b}"] = {"mean": float(np.mean(d)), "n": len(d),
                                           "n_pos": int(np.sum(np.array(d) > 0))}
    return {"per_lineage": per, "summary": summ, "pairs": pairs}


def ec_section(runs, ids, Y, dev):
    """EC leave-out: per-line skill on the dev EC lines and the H3 pluripotency genes."""
    noec = {n: runs[n] for n in runs if n.endswith("_noEC")}
    if not noec:
        return {}
    name_of = pd.read_csv(ROOT / "data" / "Model.csv", index_col=0, usecols=["ModelID", "CellLineName"]
                          )["CellLineName"].reindex(EC_DEV)
    h3 = h3_genes()
    pairs = {}
    for n in noec:
        b = n[: -len("_noEC")]
        if b in ids:
            pairs[n] = b
    lines = {}
    for n, b in pairs.items():
        for tag, rn in [("with_EC", b), ("no_EC", n)]:
            pl = pd.read_parquet(ROOT / "cache" / "runs" / ids[rn] / "per_line.parquet").set_index("ModelID")
            pl = pl.reindex(EC_DEV)
            for lid in EC_DEV:
                e = lines.setdefault(lid, {"name": name_of.get(lid), "lineage": dev.loc[lid, "OncotreeLineage"],
                                           "fold": int(dev.loc[lid, "cv_fold"])})
                e.setdefault(b, {})[tag] = {"r_line": float(pl.loc[lid, "r_line"]),
                                            "ap": float(pl.loc[lid, "ap"]),
                                            "dAP_B0": float(pl.loc[lid, "ap"] - pl.loc[lid, "ap_B0"])}
    hg = [g for g in h3.values() if g in Y.columns]
    sym = {v: k for k, v in h3.items()}
    obs = Y.loc[EC_DEV, hg]
    genes_tab = {}
    for n, b in pairs.items():
        for tag, rn in [("with_EC", b), ("no_EC", n)]:
            P = pd.read_parquet(ROOT / "cache" / "runs" / ids[rn] / "pred.parquet").loc[EC_DEV, hg]
            for g in hg:
                for lid in EC_DEV:
                    genes_tab.setdefault(sym[g], {}).setdefault(lid, {"obs": None if np.isnan(obs.loc[lid, g])
                                                                      else float(obs.loc[lid, g])})
                    genes_tab[sym[g]][lid][f"{b}:{tag}"] = float(P.loc[lid, g])
    overall = {}
    for n, b in pairs.items():
        overall[n] = {"P1_var": runs[n][0]["summary"]["variable"]["P1"],
                      "P1_var_base": runs[b][0]["summary"]["variable"]["P1"],
                      "d_vs_base": paired(fold_vals(runs[n][1], "variable", "P1"),
                                          fold_vals(runs[b][1], "variable", "P1"))}
    return {"pairs": pairs, "lines": lines, "h3_genes": genes_tab, "overall": overall,
            "h3_genes_missing": [k for k, v in h3.items() if v not in Y.columns]}


# ============================================================================ markdown
def f4(x, sign=False):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "—"
    return f"{x:+.4f}" if sign else f"{x:.4f}"


def write_md(out, tab):
    L = ["# G3 kernel track — report columns, controls, M1b tuning, M2 (D14 base), secondary splits",
         "", "Generated by `scripts/g3_report.py` from `results/registry.csv` (latest ok run per name).",
         "Dev lines only, 5 outer dev folds unless stated. P1 = mean per-gene Pearson r on variable",
         "genes within each fold, then averaged (± fold SD). ΔP1 = paired fold-wise difference; SE =",
         "SD/√5; **NB SE** = Nadeau–Bengio corrected SE (disclosure only; the selection rule uses the",
         "fold SE). Within-library / within-lineage P1: truth and prediction centred per gene within",
         "library group (configs/library_groups.csv) or lineage inside each held-out fold (definitions:",
         "reviews/g2_review/recompute_p1.py `group_centre`; lineage partialling of neg_genes.py).",
         "P2b = normalised AP on selective genes with held-out prevalence ≤ 0.5. Fit s = sum over folds",
         "(incl. inner tuning), measured with 6 threads on a machine shared with the nonlinear track:",
         "indicative only (LESSONS 14).", "", "## Results table (5-fold dev CV)", "",
         "| run | P1 var ± SD | ΔP1 vs m1_E_resB0 (SE / NB SE; folds +) | ΔP1 vs B1 | within-library P1 | "
         "within-lineage P1 | P2b nAP | fit s | λ (other/var) |",
         "|---|---|---|---|---|---|---|---|---|"]
    for r in out["runs"]:
        dh = (f"{r['dP1_head']:+.4f} ({r['dP1_head_se']:.4f} / {r['dP1_head_se_nb']:.4f}; "
              f"{r['dP1_head_pos']}/5)") if "dP1_head" in r else "—"
        L.append(f"| {r['name']} | {r['P1_var']:.4f} ± {r['P1_var_sd']:.4f} | {dh} | "
                 f"{f4(r.get('dP1_B1'), True)} | {f4(r.get('within_library'))} | "
                 f"{f4(r.get('within_lineage'))} | {f4(r['P2b_nap'])} | {r['fit_s']:.0f} | "
                 f"{r.get('lam', '—')} |")
    L += ["", "Run IDs: " + "; ".join(f"{r['name']} `{r['run_id']}`" for r in out["runs"]), ""]
    L += ["## Negative controls vs expected outcomes (stated in the configs before running)", "",
          "| control | P1 var | ΔP1 vs B1 | within-lineage | within-library | λ | checks | pass |",
          "|---|---|---|---|---|---|---|---|"]
    for c in out["controls"]:
        ch = "; ".join(f"{k}: {'ok' if v else 'FAIL'}" for k, v in c["checks"].items())
        L.append(f"| {c['run']} | {c['P1_var']:.4f} | {f4(c['dP1_B1'], True)} | {f4(c['within_lineage'])} | "
                 f"{f4(c['within_library'])} | {c.get('lam') or '—'} | {ch} | {'yes' if c['pass'] else 'NO'} |")
    if out["m1b_contrasts"]:
        L += ["", "## M1b block-weight tuning (B0 residual; training LOO coordinate descent)", "",
              "| contrast | ΔP1 | SE | NB SE | folds + | p | p Holm | helps |", "|---|---|---|---|---|---|---|---|"]
        for c in out["m1b_contrasts"]:
            L.append(f"| {c['contrast']} | {c['mean']:+.4f} | {c['se']:.4f} | {c['se_nb']:.4f} | "
                     f"{c['n_folds_pos']}/5 | {c['p']:.3g} | {c['p_holm']:.3g} | {'yes' if c['helps'] else 'no'} |")
        for c in out["m1b_c_subset"]["contrasts"]:
            L.append(f"- C common-subset rule: {c['contrast']}: ΔP1 {c['mean']:+.4f} ± {c['se']:.4f} "
                     f"({c['n_folds_pos']}/5 folds, p {c['p']:.3g})")
        L += ["", "Tuned block weights per outer fold (multiplier θ_b of the default 1/√p_b kernel; E fixed at 1):", ""]
        for r in out["runs"]:
            if r.get("weights"):
                for f, w in r["weights"].items():
                    L.append(f"- {r['name']} fold {f}: " + ", ".join(f"{k} {v:g}" for k, v in w.items())
                             + f" ({r['n_eigh'][f]} eigh)")
    if out["m2_prior_contrasts"] or any(s["candidate"].startswith("m2n_") for s in out["selection_log"]):
        L += ["", "## M2 (D14 base, pass-through for non-variable genes, nested local-λ CV)", ""]
        for c in out["m2_prior_contrasts"]:
            L.append(f"- {c['contrast']}: ΔP1 {c['mean']:+.4f} (SE {c['se']:.4f}, NB {c['se_nb']:.4f}; "
                     f"Holm p {c['p_holm']:.3g})")
        for r in out["runs"]:
            if r.get("lam_local"):
                L.append(f"- {r['name']} local λ per fold: {json.dumps(r['lam_local'])}")
    L += ["", "## Selection rule vs m1_E_resB0 (DESIGN §4: ΔP1 > 0.005 and > 1 fold SE)", "",
          "| candidate | ΔP1 | SE | NB SE | replaces |", "|---|---|---|---|---|"]
    for s in out["selection_log"]:
        L.append(f"| {s['candidate']} | {s['dP1']:+.4f} | {s['se']:.4f} | {s['se_nb']:.4f} | "
                 f"{'yes' if s['replaces'] else 'no'} |")
    L += ["", f"**Selected model after the G3 kernel track: `{out['selected_model']}`.**", ""]
    lo = out.get("lolo") or {}
    if lo:
        ms = list(lo["summary"])
        L += ["## Secondary split: leave-one-lineage-out (dev lines; lineages with ≥ 20 lines scored)", "",
              "| lineage | n | " + " | ".join(f"{m} P1 / P2b" for m in ms) + " |",
              "|---|---|" + "---|" * len(ms)]
        for lin, d in sorted(lo["per_lineage"].items(), key=lambda kv: -(kv[1]["n_lines"] or 0)):
            if not any(np.isfinite(d[m]["P1"]) for m in ms if m in d):
                continue
            L.append(f"| {lin} | {d['n_lines']} | " + " | ".join(
                f"{f4(d[m]['P1'])} / {f4(d[m]['P2b'])}" if m in d else "—" for m in ms) + " |")
        L.append("| **mean (scored)** | | " + " | ".join(
            f"{lo['summary'][m]['P1_mean_scored']:.4f} / {lo['summary'][m]['P2b_mean_scored']:.4f}" for m in ms) + " |")
        for k, v in lo["pairs"].items():
            L.append(f"- {k}: mean ΔP1 {v['mean']:+.4f} over {v['n']} lineages, positive in {v['n_pos']}")
        L.append("- Run IDs: " + "; ".join(f"{m} `{lo['summary'][m]['run_id']}`" for m in ms))
        L.append("")
    ec = out.get("ec") or {}
    if ec:
        L += ["## Secondary split: EC leave-out (4 dev germ-cell/EC lines removed from training)", "",
              "| no-EC run | P1 var | P1 var with EC | Δ (SE) |", "|---|---|---|---|"]
        for n, o in ec["overall"].items():
            L.append(f"| {n} | {o['P1_var']:.4f} | {o['P1_var_base']:.4f} | "
                     f"{o['d_vs_base']['mean']:+.4f} ({o['d_vs_base']['se']:.4f}) |")
        L += ["", "Per EC line (all eval genes): per-line r and ΔAP vs B0, with EC lines in training → without.", "",
              "| line | lineage | " + " | ".join(ec["pairs"].values()) + " |",
              "|---|---|" + "---|" * len(ec["pairs"])]
        for lid, e in ec["lines"].items():
            cells = []
            for b in ec["pairs"].values():
                w, n_ = e[b]["with_EC"], e[b]["no_EC"]
                cells.append(f"r {w['r_line']:.3f}→{n_['r_line']:.3f}; ΔAP {w['dAP_B0']:+.3f}→{n_['dAP_B0']:+.3f}")
            L.append(f"| {e['name']} ({lid}) | {e['lineage']} | " + " | ".join(cells) + " |")
        L += ["", "H3 pluripotency genes (observed effect; predicted with EC → without EC):", "",
              "| gene | line | observed | " + " | ".join(ec["pairs"].values()) + " |",
              "|---|---|---|" + "---|" * len(ec["pairs"])]
        for g, d in ec["h3_genes"].items():
            for lid, v in d.items():
                cells = [f"{v[f'{b}:with_EC']:+.2f} → {v[f'{b}:no_EC']:+.2f}" for b in ec["pairs"].values()]
                L.append(f"| {g} | {ec['lines'][lid]['name']} | {f4(v['obs'])} | " + " | ".join(cells) + " |")
        if ec["h3_genes_missing"]:
            L.append(f"- H3 genes not among the evaluable genes: {', '.join(ec['h3_genes_missing'])}")
        L.append("")
    if out["missing_runs"]:
        L += [f"Missing (no ok run yet): {', '.join(out['missing_runs'])}", ""]
    notes = ""
    if MD.exists():
        t = MD.read_text(encoding="utf-8")
        a, b = t.find("<!-- notes -->"), t.find("<!-- /notes -->")
        if a >= 0 and b > a:
            notes = t[a:b + len("<!-- /notes -->")]
    L[1:1] = ["", notes or "<!-- notes -->\n<!-- /notes -->"]
    MD.write_text("\n".join(L) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
