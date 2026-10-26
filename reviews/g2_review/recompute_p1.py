"""Independent G2 metric recomputation and line-offset diagnostic (reviewer code; dev lines only).

Does NOT import essml.evaluate / essml.metrics (LESSONS 10). CRISPR is read only through
essml.guard.load_targets with dev IDs.

    python reviews/g2_review/recompute_p1.py <dir holding cache/runs/<run_id>/pred.parquet copies>

Writes reviews/g2_review/recompute_p1.json (numbers only; no data).
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

for v in ["OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"]:
    os.environ[v] = "6"

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from essml import guard  # noqa: E402

RUNS = {
    "b0": "20261026-150901_b0_831b9415",
    "b1": "20261026-163644_b1_0e2a7a27",
    "b2": "20261026-164759_b2_289f5170",
    "m1_E": "20261026-170015_m1_E_e66eff5f",
    "r1_E": "20261026-171800_r1_E_1df71977",
    "t_E_status": "20261026-193048_t_E_status_24b46e2a",
    "m1_E_resB0": "20261026-202213_m1_E_resB0_91e52f9d",
    "m1b_E_C": "20261026-173206_m1b_E_C_09602fea",
    "m1b_full": "20261026-191141_m1b_full_7a4e925d",
}
MIN_OBS = 40


def masked_r(P, T, min_obs=MIN_OBS):
    """Per-column Pearson r over rows where T is observed; constant -> 0; n < min_obs -> NaN."""
    obs = ~np.isnan(T)
    n = obs.sum(0)
    P0 = np.where(obs, P, 0.0)
    T0 = np.where(obs, T, 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        mp = P0.sum(0) / n
        mt = T0.sum(0) / n
        dp = np.where(obs, P - mp, 0.0)
        dt = np.where(obs, T - mt, 0.0)
        den = np.sqrt((dp ** 2).sum(0) * (dt ** 2).sum(0))
        r = (dp * dt).sum(0) / den
    r = np.where(den > 1e-12, r, 0.0)
    r[n < min_obs] = np.nan
    return r


def line_centre(P, T):
    """Subtract each line's mean over its observed genes from both T and P (same genes)."""
    obs = ~np.isnan(T)
    n = obs.sum(1, keepdims=True)
    mt = np.where(obs, T, 0).sum(1, keepdims=True) / n
    mp = np.where(obs, P, 0).sum(1, keepdims=True) / n
    return P - mp, T - mt


def group_centre(P, T, labels):
    """Within each label group, subtract the per-gene group mean (observed rows) from T and P."""
    P, T = P.copy(), T.copy()
    for g in np.unique(labels):
        m = labels == g
        obs = ~np.isnan(T[m])
        n = obs.sum(0)
        with np.errstate(invalid="ignore", divide="ignore"):
            mt = np.where(obs, T[m], 0).sum(0) / n
            mp = np.where(obs, P[m], 0).sum(0) / n
        mt, mp = np.nan_to_num(mt), np.nan_to_num(mp)
        T[m] -= mt
        P[m] -= mp
    return P, T


def main(run_dir):
    run_dir = Path(run_dir)
    sp = pd.read_csv(ROOT / "configs" / "splits.csv", index_col=0)
    dev = sp[sp["role"] == "dev"]
    genes = pd.Index(pd.read_csv(ROOT / "configs" / "eval_genes.csv", dtype={"entrez": str})["entrez"])
    gc = pd.read_csv(ROOT / "configs" / "gene_classes.csv", dtype={"entrez": str}).set_index("entrez")
    gc = gc.reindex(genes)
    cls = {c: gc[c].fillna(False).astype(bool).to_numpy() for c in
           ["all", "variable", "selective", "nonessential", "common_essential"]}
    Yall = guard.load_targets(dev.index).astype(np.float64)  # dev only
    # library proxy: row NaN pattern over ALL genes; 4 largest patterns + "other"
    key = Yall.isna().apply(lambda r: hash(r.to_numpy().tobytes()), axis=1)
    top = key.value_counts().index[:4]
    lib = np.where(key.isin(top), key.map({k: i for i, k in enumerate(top)}), 4).astype(int)
    lib = pd.Series(lib, index=Yall.index)
    Y = Yall[genes]
    del Yall
    folds = dev["cv_fold"].astype(int)
    lineage = dev["OncotreeLineage"].astype(str)

    # training-fold variable sets (sensitivity: dev-wide classes use held-out lines)
    var_tr = {}
    for f in range(5):
        Ytr = Y.loc[folds.index[folds != f]].to_numpy()
        sd = np.nanstd(Ytr, axis=0, ddof=1)
        var_tr[f] = sd > 0.25

    B0 = pd.read_parquet(run_dir / RUNS["b0"] / "pred.parquet")
    out = {"lib_group_sizes": lib.value_counts().sort_index().to_dict(), "runs": {}}
    per_fold = {}
    for name, rid in RUNS.items():
        Pdf = pd.read_parquet(run_dir / rid / "pred.parquet")
        res = {}
        for f in range(5):
            ids = folds.index[folds == f]
            T = Y.loc[ids].to_numpy()
            P = Pdf.loc[ids, genes].to_numpy(np.float64)
            b0 = B0.loc[ids, genes].to_numpy(np.float64)
            obs = ~np.isnan(T)
            ev = obs.sum(0) >= MIN_OBS
            variants = {"raw": masked_r(P, T)}
            Pc, Tc = line_centre(P, T)
            variants["line_centred"] = masked_r(Pc, Tc)
            # pure line-offset model: B0 + the run's per-line mean shift over observed genes
            s = (np.where(obs, P - b0, 0).sum(1) / obs.sum(1))
            variants["offset_only"] = masked_r(b0 + s[:, None], T)
            # oracle offset: observed per-line mean (upper bound of offset-only skill)
            so = np.where(obs, T, 0).sum(1) / obs.sum(1)
            variants["oracle_offset"] = masked_r(np.broadcast_to(so[:, None], T.shape).copy(), T)
            Pl, Tl = group_centre(P, T, lib.loc[ids].to_numpy())
            variants["within_library"] = masked_r(Pl, Tl)
            Pn, Tn = group_centre(P, T, lineage.loc[ids].to_numpy())
            variants["within_lineage"] = masked_r(Pn, Tn)
            Pcl, Tcl = group_centre(Pc, Tc, lib.loc[ids].to_numpy())
            variants["line_centred_within_library"] = masked_r(Pcl, Tcl)
            for vname, r in variants.items():
                for c, m in cls.items():
                    res.setdefault(vname, {}).setdefault(c, []).append(float(np.nanmean(r[m & ev])))
            # class-wise per-line centring: line mean over the evaluated class's own genes
            for c, m in cls.items():
                mm = m & ev
                Pk, Tk = line_centre(P[:, mm], T[:, mm])
                res.setdefault("line_centred_classwise", {}).setdefault(c, []).append(
                    float(np.nanmean(masked_r(Pk, Tk))))
                sk = (np.where(obs[:, mm], P[:, mm] - b0[:, mm], 0).sum(1) / obs[:, mm].sum(1))
                res.setdefault("offset_only_classwise", {}).setdefault(c, []).append(
                    float(np.nanmean(masked_r(b0[:, mm] + sk[:, None], T[:, mm]))))
            res.setdefault("raw", {}).setdefault("variable_trainfold", []).append(
                float(np.nanmean(variants["raw"][var_tr[f] & ev])))
            per_fold.setdefault(name, {})[f] = variants["raw"]
            # line-mean diagnostics
            pm = np.where(obs, P, 0).sum(1) / obs.sum(1)
            neg = cls["nonessential"] & ev
            nm_obs = np.nanmean(np.where(obs[:, neg], T[:, neg], np.nan), axis=1)
            nm_pred = np.nanmean(np.where(obs[:, neg], P[:, neg], np.nan), axis=1)
            res.setdefault("diag", {}).setdefault("r_linemean_pred_obs_all", []).append(
                float(np.corrcoef(pm, so)[0, 1]))
            res["diag"].setdefault("r_linemean_pred_obs_neg", []).append(
                float(np.corrcoef(nm_pred, nm_obs)[0, 1]))
            res["diag"].setdefault("sd_linemean_pred_all", []).append(float(np.std(pm)))
            res["diag"].setdefault("sd_linemean_obs_all", []).append(float(np.std(so)))
        summ = {v: {c: {"mean": float(np.mean(x)), "sd": float(np.std(x, ddof=1)), "folds": x}
                    for c, x in d.items()} for v, d in res.items()}
        out["runs"][name] = summ
        print(name, {v: round(summ[v]["variable"]["mean"], 4) for v in summ if "variable" in summ[v]},
              "| NEG raw", round(summ["raw"]["nonessential"]["mean"], 4),
              "NEG line-centred", round(summ["line_centred"]["nonessential"]["mean"], 4), flush=True)

    # paired fold differences on raw variable P1
    def delta(a, b, cl="variable", v="raw"):
        x = np.array(out["runs"][a][v][cl]["folds"]) - np.array(out["runs"][b][v][cl]["folds"])
        return {"d": float(x.mean()), "se": float(x.std(ddof=1) / np.sqrt(len(x))),
                "folds_pos": int((x > 0).sum()), "folds": x.round(5).tolist()}
    pairs = [("m1_E", "b1"), ("m1_E", "b2"), ("b2", "b1"), ("m1_E_resB0", "m1_E"), ("r1_E", "m1_E"),
             ("t_E_status", "m1_E"), ("t_E_status", "m1_E_resB0"), ("m1_E_resB0", "r1_E"),
             ("m1_E_resB0", "b1"), ("r1_E", "b1"), ("t_E_status", "b1"), ("m1b_E_C", "m1_E")]
    out["deltas_raw_variable"] = {f"{a}-{b}": delta(a, b) for a, b in pairs}
    out["deltas_line_centred_variable"] = {f"{a}-{b}": delta(a, b, v="line_centred") for a, b in pairs}
    out["deltas_within_library_variable"] = {f"{a}-{b}": delta(a, b, v="within_library") for a, b in pairs}
    out["deltas_within_lineage_variable"] = {f"{a}-{b}": delta(a, b, v="within_lineage") for a, b in pairs}
    # cross-check against the committed aggregated per-gene table (mean r over folds per gene)
    agg = {}
    for name, rid in RUNS.items():
        pg = pd.read_csv(ROOT / "results" / "runs" / rid / "per_gene.csv.gz", dtype={"gene": str})
        pg = pg.set_index("gene").reindex(genes)
        mine = np.nanmean(np.vstack([per_fold[name][f] for f in range(5)]), axis=0)
        ok = np.isfinite(mine) & np.isfinite(pg["r_mean"].to_numpy())
        agg[name] = {"max_abs_diff_gene_mean_r": float(np.max(np.abs(mine[ok] - pg["r_mean"].to_numpy()[ok]))),
                     "P1var_from_per_gene_csv": float(pg.loc[cls["variable"], "r_mean"].mean())}
    out["per_gene_csv_check"] = agg
    (Path(__file__).parent / "recompute_p1.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out["deltas_raw_variable"], indent=0))
    print(json.dumps(agg, indent=0))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ROOT / "cache" / "runs")
