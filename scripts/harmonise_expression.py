"""External expression harmonisation, end to end (DESIGN v3 section 7.3).

Reads
  configs/expression_harmonisation.json   sources, conversion parameters, pre-registered selection rule
  configs/external_expr_profiles.csv      which samples are used, their state / screen / DepMap match
  configs/splits.csv                      roles (dev lines = reference, NN pool and fitting population)
  cache/expression.parquet                DepMap 26Q3 log2(TPM+1), protein-coding (via essml.io)
  Ensembl 116 GTF, HGNC                   gene lengths and ID mapping
Never reads screen data and never applies a prediction model.

Writes
  cache/external_expr/log2tpm_unharmonised.parquet     all profiles (samples + state means), DepMap columns
  cache/external_expr/harmonised_<method>.parquet      each candidate applied to the non-transformed profiles
  cache/external_expr/harmonised.parquet               chosen method (rows = profile id, cols = Entrez)
  results/expression_harmonisation/
      source_qc.csv            per sample: genes observed, PC TPM sum, counts with a length, ...
      count_path_check.csv     ENCODE expected_count -> our TPM vs ENCODE RSEM TPM
      method_comparison.csv    per validation state x candidate x path: NN rank/recall, profile agreement
      method_summary.csv       aggregated per candidate x stratum x path
      selection.json           rule outcome
      profiles.csv             all profiles: AD distance, nearest lineage, flags
      ad_reference.json        AD distribution of DepMap dev lines (vs other dev lines)

Usage: python scripts/harmonise_expression.py [--record]   (--record writes chosen_method into the config)
"""
import argparse
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from essml import expr_harmonise as H  # noqa: E402
from essml import priors  # noqa: E402
from essml.io import load_gene_map, load_matrix  # noqa: E402

CFG = ROOT / "configs" / "expression_harmonisation.json"
PROFILES = ROOT / "configs" / "external_expr_profiles.csv"
SPLITS = ROOT / "configs" / "splits.csv"
OUT = ROOT / "results" / "expression_harmonisation"
CACHE = ROOT / "cache" / "external_expr"
METHODS = ["none", "quantile", "quantile_location", "rank"]


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def tpm_mean_log(x: pd.DataFrame) -> pd.Series:
    """State profile: mean TPM across replicates, back to log2(TPM+1) (NaN if any replicate is NaN)."""
    return np.log2((np.power(2.0, x) - 1.0).mean(axis=0, skipna=False) + 1.0)


def build_profiles(cfg, prof, dep, mapper, genes, pc_target, fragment_length, encode_counts_path=False):
    """Convert every catalogued sample -> log2(TPM+1) on DepMap columns. Returns (X, qc rows, encode count tables)."""
    cols = dep.columns
    proj_ens = mapper.ensembl_projection(cols)
    X, qc, extra = [], [], {}
    for src, grp in prof.groupby("source", sort=False):
        if src == "depmap_native":
            x = dep.loc[grp["sample"]].astype(float)
            x.index = grp["profile_id"].to_numpy()
            X.append(x)
            for pid in x.index:
                qc.append(dict(profile_id=pid, source=src, n_observed=int(x.loc[pid].notna().sum()),
                               pc_tpm_sum=float((2 ** x.loc[pid] - 1).sum())))
            continue
        spec = cfg["sources"][src]
        samples = list(grp["sample"])
        log(f"parse {src} ({len(samples)} samples)")
        raw = H.parse_source(spec, samples)
        tpm, info = H.to_tpm(raw, genes, mapper, fragment_length)
        if "expected_count" in raw.extra:
            extra[src] = raw.extra["expected_count"]
        if raw.id_type == "symbol" and raw.unit != "counts":
            proj = mapper.symbol_projection(tpm.index, cols)
        else:
            proj = proj_ens
        restricted = spec["universe"] == "restricted"
        x = H.project_to_depmap(tpm, proj, cols, absent=spec["absent"], pc_target=pc_target if restricted else None)
        x.index = grp["profile_id"].to_numpy()
        X.append(x)
        pcsum = (2 ** x - 1).sum(axis=1)
        for s, pid in zip(samples, x.index):
            row = dict(profile_id=pid, source=src, unit=raw.unit, id_type=raw.id_type,
                       universe=spec["universe"], tpm_compatible=spec["tpm_compatible"],
                       n_rows_raw=int(raw.values.shape[0]), n_observed=int(x.loc[pid].notna().sum()),
                       pc_tpm_sum=float(pcsum[pid]))
            for k, v in info.items():
                row[k] = v.get(s) if isinstance(v, dict) else v
            qc.append(row)
    return pd.concat(X), pd.DataFrame(qc), extra


def add_states(X: pd.DataFrame, prof: pd.DataFrame) -> pd.DataFrame:
    """Append one 'state:<state>' row per state with >= 2 samples (mean TPM)."""
    rows = {}
    for st, g in prof.groupby("state", sort=False):
        if len(g) >= 2:
            rows[f"state:{st}"] = tpm_mean_log(X.loc[g["profile_id"]])
    return pd.concat([X, pd.DataFrame(rows).T]) if rows else X


def validation_states(X, prof):
    """One row per validation state (mean TPM over replicates) + its DepMap model and match class."""
    v = prof[prof["kind"] == "validation"]
    rows, meta = {}, []
    for st, g in v.groupby("state", sort=False):
        rows[st] = tpm_mean_log(X.loc[g["profile_id"]]) if len(g) > 1 else X.loc[g["profile_id"].iloc[0]]
        meta.append(dict(state=st, source=g["source"].iloc[0], depmap_model=g["depmap_model"].iloc[0],
                         match=g["match"].iloc[0], n_rep=len(g)))
    return pd.DataFrame(rows).T, pd.DataFrame(meta).set_index("state")


def evaluate(V, vmeta, dep, dev, hvg, method, path_label, loso=False):
    """Held-out evaluation of one candidate on validation states. Returns per-state rows."""
    exact = vmeta.index[(vmeta["match"] == "exact") & vmeta["depmap_model"].isin(dev)]
    D = dep.loc[vmeta["depmap_model"]].astype(float)
    D.index = vmeta.index
    pool = dep.loc[dev].astype(float)
    sorted_pool = np.sort(pool.to_numpy(), axis=1)
    ref_all = sorted_pool.mean(axis=0)
    pos = {m: i for i, m in enumerate(dev)}
    # ---------- transform external states (+ pool and DepMap counterpart in the same space)
    if method == "none":
        E = V.copy()
    elif method in ("quantile", "quantile_location"):
        rows = []
        for st in V.index:
            m = vmeta.at[st, "depmap_model"]
            ref = ref_all
            if m in pos:   # quantile reference without the matched line's own DepMap profile
                n = sorted_pool.shape[0]
                ref = (ref_all * n - sorted_pool[pos[m]]) / (n - 1)
            rows.append(H.quantile_normalise(V.loc[[st]], ref))
        E = pd.concat(rows)
        if method == "quantile_location":
            if loso:
                groups = vmeta["source"]
            else:
                groups = vmeta["depmap_model"]
            # offsets are fitted on exact-match dev states only, always excluding the held-out group
            audit = []
            E = H.heldout_location(E, D, groups.loc[E.index], fit_rows=E.index.isin(exact), audit=audit)
            for g, tr in audit:     # the held-out line's own DepMap profile never enters its fit
                held = set(vmeta.loc[groups == g, "depmap_model"])
                assert not held & set(vmeta.loc[tr, "depmap_model"]), g
    elif method == "rank":
        E = H.rank_transform(V)
        D = H.rank_transform(D)
        pool = H.rank_transform(pool)
    else:
        raise ValueError(method)
    mu = pool.mean(axis=0)
    rows = []
    for st in E.index:
        m = vmeta.at[st, "depmap_model"]
        P = pool if m in pos else pd.concat([pool, (H.rank_transform(dep.loc[[m]].astype(float))
                                                   if method == "rank" else dep.loc[[m]].astype(float))])
        r = H.corr_to_pool(E.loc[[st], hvg], P[hvg])
        nn = H.nn_recall(r, pd.Series({st: m})).iloc[0].to_dict()
        rc = H.corr_to_pool((E.loc[[st], hvg] - mu[hvg]), (P[hvg] - mu[hvg]))
        nnc = H.nn_recall(rc, pd.Series({st: m})).iloc[0]
        e, d = E.loc[st], D.loc[st]
        ok = e.notna() & d.notna()
        rows.append(dict(state=st, method=method + ("_loso" if loso else ""), path=path_label,
                         source=vmeta.at[st, "source"], match=vmeta.at[st, "match"], depmap_model=m,
                         pool_size=P.shape[0], true_rank=nn["true_rank"], recall1=nn["recall1"],
                         recall5=nn["recall5"], recall10=nn["recall10"], r_true_hvg=nn["r_true"],
                         nearest=nn["nearest"], true_rank_centred=int(nnc["true_rank"]),
                         recall1_centred=bool(nnc["recall1"]),
                         r_profile_all=float(np.corrcoef(e[ok], d[ok])[0, 1]),
                         r_profile_centred=float(np.corrcoef((e - mu)[ok], (d - mu)[ok])[0, 1]),
                         mae=float((e[ok] - d[ok]).abs().mean()), n_genes=int(ok.sum())))
    return rows


def summarise(cmp):
    g = cmp.groupby(["path", "match", "method"])
    s = g.agg(n=("state", "size"), recall1=("recall1", "mean"), recall5=("recall5", "mean"),
              recall10=("recall10", "mean"), median_true_rank=("true_rank", "median"),
              recall1_centred=("recall1_centred", "mean"), r_profile_all=("r_profile_all", "mean"),
              r_profile_centred=("r_profile_centred", "mean"), mae=("mae", "mean"))
    return s.reset_index()


def select(summary, rule):
    s = summary[(summary["path"] == "reported") & (summary["match"] == "exact")].set_index("method")
    cand = [m for m in METHODS if m in s.index]
    excluded = {}
    if "quantile_location_loso" in s.index and "quantile_location" in s.index:
        drop = s.at["quantile_location", "recall1"] - s.at["quantile_location_loso", "recall1"]
        if drop > rule["tie_margin"]:
            cand.remove("quantile_location")
            excluded["quantile_location"] = f"LOSO recall@1 drop {drop:.3f} > tie_margin"
    best = max(s.loc[cand, "recall1"])
    ties = [m for m in cand if s.at[m, "recall1"] >= best - rule["tie_margin"]]
    order = rule["simplicity_order"]
    ties.sort(key=order.index)
    choice = ties[0]
    for m in ties[1:]:
        if s.at[m, "r_profile_all"] > s.at[choice, "r_profile_all"] + 0.01:
            choice = m
    return dict(best_recall1=float(best), ties=ties, excluded=excluded, chosen=choice,
                table=s[["n", "recall1", "recall1_centred", "median_true_rank", "r_profile_all",
                         "r_profile_centred", "mae"]].round(4).reset_index().to_dict(orient="records"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--record", action="store_true", help="write the rule's choice into the config")
    args = ap.parse_args()
    t0 = time.time()
    cfg = json.load(open(CFG, encoding="utf-8"))
    ev = cfg["evaluation"]
    prof = pd.read_csv(PROFILES, dtype=str).fillna("")
    splits = pd.read_csv(SPLITS)
    dev = list(splits.loc[splits["role"] == "dev", "ModelID"])
    lineage = pd.read_csv(ROOT / "data" / "Model.csv", usecols=["ModelID", "OncotreeLineage"]).set_index("ModelID")["OncotreeLineage"]
    dep = load_matrix("expression")
    dev = [m for m in dev if m in dep.index]
    assert not dep.loc[dev].isna().any().any(), "DepMap dev expression has NaN"
    pc_target = float((2.0 ** dep.loc[dev].astype(float) - 1).sum(axis=1).median())
    log(f"{len(dev)} dev lines; DepMap dev median PC TPM sum {pc_target:.0f}")
    gm = load_gene_map()
    genes = H.load_gene_table()
    mapper = H.GeneMapper(priors.load_hgnc(), genes, dict(zip(gm["entrez"], gm["symbol"])))
    hvg = H.top_variable_genes(dep.loc[dev].astype(float), ev["n_hvg"])
    OUT.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)

    # ---------------------------------------------------------------- conversion
    frag = float(cfg["conversion"]["counts_fragment_length"])
    X, qc, extra = build_profiles(cfg, prof, dep, mapper, genes, pc_target, frag)
    proj = mapper.ensembl_projection(dep.columns)
    qc["n_depmap_cols"] = dep.shape[1]
    qc.to_csv(OUT / "source_qc.csv", index=False)
    Xs = add_states(X, prof)
    Xs.astype("float32").to_parquet(CACHE / "log2tpm_unharmonised.parquet")
    log(f"converted {X.shape[0]} samples; {Xs.shape[0] - X.shape[0]} state means")

    # ---------------------------------------------------------------- count-path check (ENCODE)
    enc = prof[prof["source"] == "encode_geo"]
    cnt = extra["encode_geo"]
    rows = []
    paths = {}
    for fl in (0.0, 200.0):
        tpm_c, used = H.counts_to_tpm(cnt.loc[cnt.index.str.match(r"^ENSG\d{11}$")], genes["length"], fl)
        xc = H.project_to_depmap(tpm_c, proj, dep.columns)
        xc.index = enc["profile_id"].to_numpy()
        paths[fl] = xc
        for pid in xc.index:
            xa, xb = xc.loc[pid], X.loc[pid]
            ok = xa.notna() & xb.notna()
            rows.append(dict(profile_id=pid, fragment_length=fl, r_vs_rsem_tpm=float(np.corrcoef(xa[ok], xb[ok])[0, 1]),
                             mean_diff_log2=float((xa[ok] - xb[ok]).mean()), mae_log2=float((xa[ok] - xb[ok]).abs().mean()),
                             frac_counts_with_length=float(used[enc.set_index('profile_id').at[pid, 'sample']])))
    cpc = pd.DataFrame(rows)
    cpc.to_csv(OUT / "count_path_check.csv", index=False)
    log("count path check:\n" + cpc.groupby("fragment_length")[["r_vs_rsem_tpm", "mae_log2", "mean_diff_log2"]].mean().round(4).to_string())

    # ---------------------------------------------------------------- unit check by gene-length dependence
    # slope of log2(profile / ENCODE H1 RSEM TPM) on log2(gene length), genes > 1 log2 unit in both:
    # ~0 for TPM-like units, ~+1 for CPM-like (length-proportional) units.
    ens2len = genes["length"]
    ent_len = pd.Series({ent: ens2len.get(ens) for ens, ent in proj.items()}).reindex(dep.columns)
    ref_h1 = Xs.loc["state:ENCODE_H1"]
    lrows = []
    for pid in [i for i in Xs.index if i.startswith("state:") and i != "state:ENCODE_H1"]:
        x_ = Xs.loc[pid]
        ok = (x_ > 2) & (ref_h1 > 2) & ent_len.notna()
        sl = np.polyfit(np.log2(ent_len[ok].astype(float)), (x_ - ref_h1)[ok], 1)[0]
        lrows.append(dict(profile_id=pid, slope_log2ratio_vs_log2length=round(float(sl), 3), n_genes=int(ok.sum())))
    pd.DataFrame(lrows).to_csv(OUT / "length_bias_check.csv", index=False)

    # ---------------------------------------------------------------- method comparison
    V, vmeta = validation_states(X, prof)
    Xc = X.copy()
    Xc.loc[paths[frag].index] = paths[frag]
    Vc, _ = validation_states(Xc, prof)
    cmp_rows = []
    for label, VV in (("reported", V), ("encode_from_counts", Vc)):
        for m in METHODS:
            cmp_rows += evaluate(VV, vmeta, dep, dev, hvg, m, label)
        cmp_rows += evaluate(VV, vmeta, dep, dev, hvg, "quantile_location", label, loso=True)
    cmp = pd.DataFrame(cmp_rows)
    cmp.to_csv(OUT / "method_comparison.csv", index=False)
    summary = summarise(cmp)
    summary.to_csv(OUT / "method_summary.csv", index=False)
    log("summary (exact dev matches):\n" + summary[summary["match"] == "exact"].round(3).to_string())
    sel = select(summary, cfg["selection_rule"])
    sel["note"] = "rule from configs/expression_harmonisation.json selection_rule; population = exact dev matches, 'reported' path"
    rule_outcome = sel["chosen"]
    sel["rule_outcome"] = rule_outcome
    log(f"selection rule -> {rule_outcome} (ties {sel['ties']}, excluded {sel['excluded']})")
    ov = cfg.get("pi_override")
    if ov:
        # explicit, recorded PI override: the rule and its outcome stay untouched and are re-checked here
        if ov["rule_outcome"] != rule_outcome:
            raise SystemExit(f"pi_override was made against rule outcome {ov['rule_outcome']!r}, "
                             f"but the rule now selects {rule_outcome!r}: the PI must revisit")
        if cfg.get("chosen_method") != ov["chosen_method"]:
            raise SystemExit("config chosen_method must equal pi_override.chosen_method")
        if args.record:
            raise SystemExit("--record is not allowed while a pi_override is present")
        chosen = ov["chosen_method"]
        sel["pi_override"] = ov
        log(f"PI override ({ov['decision']}): using {chosen} instead of rule outcome {rule_outcome}")
    else:
        chosen = rule_outcome
    sel["chosen"] = chosen
    json.dump(sel, open(OUT / "selection.json", "w"), indent=1)
    if ov:
        pass
    elif args.record:
        txt = open(CFG, encoding="utf-8", newline="").read()   # textual edit keeps formatting/line endings
        txt, n = re.subn(r'"chosen_method": [^\n]*', f'"chosen_method": "{chosen}"', txt)
        assert n == 1
        json.loads(txt)
        open(CFG, "w", encoding="utf-8", newline="").write(txt)
        cfg["chosen_method"] = chosen
    elif cfg.get("chosen_method") not in (None, chosen):
        raise SystemExit(f"config chosen_method={cfg['chosen_method']} but the rule selects {chosen}")

    # ---------------------------------------------------------------- apply to all profiles
    exact = vmeta.index[(vmeta["match"] == "exact") & vmeta["depmap_model"].isin(dev)]
    pool = dep.loc[dev].astype(float)
    ref_q = H.reference_quantiles(pool.to_numpy())
    native = set(prof.loc[prof["source"] == "depmap_native", "profile_id"]) | {
        f"state:{s}" for s in prof.loc[prof["source"] == "depmap_native", "state"]}
    ext_ids = [i for i in Xs.index if i not in native]
    Q = H.quantile_normalise(Xs.loc[ext_ids], ref_q)
    Vq = H.quantile_normalise(V.loc[exact], ref_q)
    Dx = dep.loc[vmeta.loc[exact, "depmap_model"]].astype(float)
    Dx.index = exact
    offsets = H.fit_location(Vq, Dx)
    offsets.rename("offset").to_frame().to_parquet(CACHE / "location_offsets.parquet")
    spaces = {"none": Xs.loc[ext_ids], "quantile": Q, "quantile_location": H.apply_location(Q, offsets),
              "rank": H.rank_transform(Xs.loc[ext_ids])}
    for m, M in spaces.items():
        nat = Xs.loc[[i for i in Xs.index if i in native]]
        full = pd.concat([M, H.rank_transform(nat) if m == "rank" else nat])
        full.astype("float32").to_parquet(CACHE / f"harmonised_{m}.parquet")
    final = pd.concat([spaces[chosen], Xs.loc[[i for i in Xs.index if i in native]]])
    final.astype("float32").to_parquet(CACHE / "harmonised.parquet")

    # ---------------------------------------------------------------- applicability domain
    poolm = H.rank_transform(pool) if chosen == "rank" else pool
    pr = np.corrcoef(poolm[hvg].to_numpy())
    pr = pd.DataFrame(pr, index=dev, columns=dev)
    self_ad = H.pool_self_ad(pr, ev["ad_k"])
    r = H.corr_to_pool(final[hvg], poolm[hvg])
    ad = H.applicability(r, lineage, ev["ad_k"])
    ad["ad_percentile_vs_dev"] = [float((self_ad < d).mean() * 100) for d in ad["ad_distance"]]
    # sensitivity: the same statistic after quantile normalisation (pool = DepMap dev, unchanged)
    qspace = pd.read_parquet(CACHE / "harmonised_quantile.parquet").loc[final.index]
    rq = H.corr_to_pool(qspace[hvg], pool[hvg])
    ad["ad_distance_quantile"] = H.applicability(rq, lineage, ev["ad_k"])["ad_distance"]
    json.dump({"method": chosen, "k": ev["ad_k"], "n_dev": len(dev),
               "dev_self_ad_quantiles": {str(q): float(self_ad.quantile(q)) for q in (0.05, 0.25, 0.5, 0.75, 0.95, 0.99)},
               "dev_self_ad_max": float(self_ad.max())}, open(OUT / "ad_reference.json", "w"), indent=1)

    # ---------------------------------------------------------------- profiles table
    meta = prof.set_index("profile_id")
    qci = qc.set_index("profile_id")
    srcl = {k: v["label"] for k, v in cfg["sources"].items()}
    srcl["depmap_native"] = "DepMap 26Q3 (native pipeline)"
    out = []
    for pid in final.index:
        if pid.startswith("state:"):
            st = pid[len("state:"):]
            g = meta[meta["state"] == st]
            base = g.iloc[0]
            n = len(g)
            obs = int(final.loc[pid].notna().sum())
            pcs = float((2 ** Xs.loc[pid] - 1).sum())
        else:
            base = meta.loc[pid]
            n = 1
            obs = int(qci.at[pid, "n_observed"])
            pcs = float(qci.at[pid, "pc_tpm_sum"])
        spec = cfg["sources"].get(base["source"], {"tpm_compatible": True, "universe": "full"})
        out.append(dict(profile_id=pid, level="state" if pid.startswith("state:") else "sample", n_samples=n,
                        kind=base["kind"], state=base["state"], screens=base["screens"], cell_type=base["cell_type"],
                        source=srcl.get(base["source"], base["source"]), source_id=base["source"],
                        tpm_compatible=bool(spec["tpm_compatible"]), universe=spec["universe"],
                        depmap_match=base["depmap_model"] if base["kind"] == "validation" else "",
                        match=base["match"], n_genes_observed=obs, pc_tpm_sum_before_harmonisation=round(pcs),
                        method=("native" if pid in native else chosen),
                        ad_distance=round(ad.at[pid, "ad_distance"], 4),
                        ad_percentile_vs_dev=round(ad.at[pid, "ad_percentile_vs_dev"], 1),
                        ad_distance_quantile=round(ad.at[pid, "ad_distance_quantile"], 4),
                        nearest_lineage=ad.at[pid, "nearest_lineage"],
                        nearest_lineage_n_of_10=int(ad.at[pid, "nearest_lineage_n"]),
                        nearest_line=ad.at[pid, "nearest_line"],
                        nearest_line_lineage=ad.at[pid, "nearest_line_lineage"],
                        note=base["note"]))
    pt = pd.DataFrame(out)
    pt.to_csv(OUT / "profiles.csv", index=False)
    log(f"profiles.csv: {len(pt)} rows; runtime {time.time() - t0:.0f} s")
    print(pt[pt["level"] == "state"][["profile_id", "tpm_compatible", "ad_distance", "ad_percentile_vs_dev",
                                       "nearest_lineage", "nearest_lineage_n_of_10", "nearest_line"]].to_string())


if __name__ == "__main__":
    main()
