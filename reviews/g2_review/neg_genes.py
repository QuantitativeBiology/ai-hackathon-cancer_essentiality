"""Hart NEGv1 puzzle (D16), dev lines only: where does M1's per-gene r on non-essential genes come from?

    python reviews/g2_review/neg_genes.py <dir with run pred.parquet copies>

Checks (within each outer fold, then averaged):
  - distribution of per-gene r (broad shift vs a few genes);
  - r(own-gene expression, gene effect) and r(own relative CN, gene effect) on observed lines;
  - M1 skill left after partialling own expression + own CN + lineage means out of BOTH the
    observed and predicted values (per gene, within fold).
Writes reviews/g2_review/neg_genes.json.
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
sys.path.insert(0, str(Path(__file__).parent))
from essml import guard  # noqa: E402
from recompute_p1 import RUNS, masked_r  # noqa: E402


def partial_out(A, Xs, obs):
    """Per column g: residual of A[:, g] on [1, Xs[k][:, g] ...] over observed rows (OLS)."""
    out = np.full_like(A, np.nan)
    for g in range(A.shape[1]):
        o = obs[:, g]
        D = np.column_stack([np.ones(o.sum())] + [x[o, g] for x in Xs])
        beta, *_ = np.linalg.lstsq(D, A[o, g], rcond=None)
        out[o, g] = A[o, g] - D @ beta
    return out


def main(run_dir):
    run_dir = Path(run_dir)
    sp = pd.read_csv(ROOT / "configs" / "splits.csv", index_col=0)
    dev = sp[sp["role"] == "dev"]
    gc = pd.read_csv(ROOT / "configs" / "gene_classes.csv", dtype={"entrez": str}).set_index("entrez")
    eg = pd.Index(pd.read_csv(ROOT / "configs" / "eval_genes.csv", dtype={"entrez": str})["entrez"])
    neg = gc.index[gc["nonessential"].astype(bool)].intersection(eg)
    var = gc.index[gc["variable"].astype(bool)].intersection(eg)
    Y = guard.load_targets(dev.index)[neg].astype(np.float64)
    E = pd.read_parquet(ROOT / "cache" / "features" / "E.parquet").reindex(index=dev.index)
    CN = pd.read_parquet(ROOT / "cache" / "features" / "C.parquet").reindex(index=dev.index)
    have_e = neg.intersection(E.columns)
    out = {"n_neg_eval": len(neg), "n_neg_with_expr": len(have_e),
           "n_neg_variable": int(len(neg.intersection(var))),
           "neg_sd_median": float(gc.loc[neg, "sd"].median()),
           "var_sd_median": float(gc.loc[var, "sd"].median())}
    Ef = E.reindex(columns=neg).to_numpy(np.float64)
    Cf = CN.reindex(columns=neg).fillna(0).to_numpy(np.float64)
    Ef = np.where(np.isnan(Ef), np.nanmean(Ef, 0), Ef)
    folds = dev["cv_fold"].astype(int)
    lin = dev["OncotreeLineage"].astype(str)
    res = {}
    for name in ["b1", "m1_E", "m1_E_resB0"]:
        P_all = pd.read_parquet(run_dir / RUNS[name] / "pred.parquet")[neg]
        for f in range(5):
            ids = folds.index[folds == f]
            ix = dev.index.get_indexer(ids)
            T = Y.loc[ids].to_numpy()
            P = P_all.loc[ids].to_numpy(np.float64)
            obs = ~np.isnan(T)
            e, c = Ef[ix], Cf[ix]
            r = masked_r(P, T)
            d = res.setdefault(name, {})
            d.setdefault("P1_neg", []).append(float(np.nanmean(r)))
            d.setdefault("median_r", []).append(float(np.nanmedian(r)))
            d.setdefault("frac_r_gt_0", []).append(float(np.nanmean(r > 0)))
            d.setdefault("frac_r_gt_0.2", []).append(float(np.nanmean(r > 0.2)))
            if name == "m1_E":
                d.setdefault("r_obs_ownexpr", []).append(float(np.nanmean(masked_r(e, T))))
                d.setdefault("r_obs_owncn", []).append(float(np.nanmean(masked_r(c, T))))
                d.setdefault("r_pred_ownexpr", []).append(float(np.nanmean(masked_r(e, np.where(obs, P, np.nan)))))
                d.setdefault("r_pred_owncn", []).append(float(np.nanmean(masked_r(c, np.where(obs, P, np.nan)))))
            L = pd.get_dummies(lin.loc[ids]).to_numpy(np.float64)
            Lcols = [np.repeat(L[:, [k]], T.shape[1], axis=1) for k in range(L.shape[1] - 1)]
            for tag, Xs in [("own_expr_cn", [e, c]), ("own_expr_cn_lineage", [e, c] + Lcols)]:
                Tr = partial_out(np.where(obs, T, 0), Xs, obs)
                Pr = partial_out(P, Xs, obs)
                d.setdefault(f"P1_neg_after_partial_{tag}", []).append(float(np.nanmean(masked_r(Pr, Tr))))
    out["runs"] = {n: {k: {"mean": float(np.mean(v)), "folds": np.round(v, 4).tolist()} for k, v in d.items()}
                   for n, d in res.items()}
    (Path(__file__).parent / "neg_genes.json").write_text(json.dumps(out, indent=1))
    print(json.dumps({k: v for k, v in out.items() if k != "runs"}))
    for n, d in out["runs"].items():
        print(n, {k: round(v["mean"], 4) for k, v in d.items()})


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ROOT / "cache" / "runs")
