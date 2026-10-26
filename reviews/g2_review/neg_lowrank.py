"""Is M1's skill on Hart NEGv1 genes explained by a few shared line-level factors? (dev only)

Per outer fold: PCA loadings of the NEG-gene block are computed on the TRAINING lines' observed
gene effects (NaN -> training gene mean). Held-out observed values and predictions are centred by
the training gene means and the top-k loadings are projected out of both; P1 on NEG genes is
recomputed. Also reports variance explained by the top PCs and how much of the held-out PC1-3
scores is explained by library group (row NaN pattern) and lineage, and how well M1 predicts them.

    python reviews/g2_review/neg_lowrank.py <dir with run pred.parquet copies>
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

KS = [1, 2, 3, 5, 10, 20]


def r2(y, labels):
    """R^2 of y on group means (in-sample)."""
    s = pd.Series(y)
    fit = s.groupby(np.asarray(labels)).transform("mean").to_numpy()
    return float(1 - ((y - fit) ** 2).sum() / ((y - y.mean()) ** 2).sum())


def main(run_dir, cls_name="nonessential"):
    run_dir = Path(run_dir)
    sp = pd.read_csv(ROOT / "configs" / "splits.csv", index_col=0)
    dev = sp[sp["role"] == "dev"]
    gc = pd.read_csv(ROOT / "configs" / "gene_classes.csv", dtype={"entrez": str}).set_index("entrez")
    eg = pd.Index(pd.read_csv(ROOT / "configs" / "eval_genes.csv", dtype={"entrez": str})["entrez"])
    genes = gc.index[gc[cls_name].astype(bool)].intersection(eg)
    Yall = guard.load_targets(dev.index)
    key = Yall.isna().apply(lambda r: hash(r.to_numpy().tobytes()), axis=1)
    top = key.value_counts().index[:4]
    lib = pd.Series(np.where(key.isin(top), key.map({k: i for i, k in enumerate(top)}), 4), index=dev.index)
    Y = Yall[genes].astype(np.float64)
    del Yall
    folds = dev["cv_fold"].astype(int)
    lin = dev["OncotreeLineage"].astype(str)
    preds = {n: pd.read_parquet(run_dir / RUNS[n] / "pred.parquet")[genes] for n in ["b1", "m1_E", "m1_E_resB0"]}
    res = {n: {k: [] for k in [0] + KS} for n in preds}
    diag = {"var_expl_train": [], "lib_r2_pc": [], "lineage_r2_pc": [], "m1_r_pc_scores": []}
    for f in range(5):
        tr, ho = folds.index[folds != f], folds.index[folds == f]
        Ytr = Y.loc[tr].to_numpy()
        mu = np.nanmean(Ytr, 0)
        A = np.where(np.isnan(Ytr), mu, Ytr) - mu
        U, s, Vt = np.linalg.svd(A, full_matrices=False)
        diag["var_expl_train"].append((s[:20] ** 2 / (s ** 2).sum()).round(4).tolist())
        T = Y.loc[ho].to_numpy()
        obs = ~np.isnan(T)
        Tc = np.where(obs, T, mu) - mu
        sc = Tc @ Vt[:3].T
        diag["lib_r2_pc"].append([r2(sc[:, j], lib.loc[ho]) for j in range(3)])
        diag["lineage_r2_pc"].append([r2(sc[:, j], lin.loc[ho]) for j in range(3)])
        for n, Pdf in preds.items():
            P = Pdf.loc[ho].to_numpy(np.float64)
            Pc = P - mu
            if n == "m1_E":
                psc = np.where(obs, Pc, 0) @ Vt[:3].T
                diag["m1_r_pc_scores"].append([float(np.corrcoef(psc[:, j], sc[:, j])[0, 1]) for j in range(3)])
            for k in [0] + KS:
                Vk = Vt[:k]
                Tk = Tc - (Tc @ Vk.T) @ Vk if k else Tc
                Pk = Pc - (Pc @ Vk.T) @ Vk if k else Pc
                res[n][k].append(float(np.nanmean(masked_r(Pk, np.where(obs, Tk, np.nan)))))
    out = {"class": cls_name, "n_genes": len(genes),
           "P1_after_removing_top_k_train_PCs": {n: {str(k): round(float(np.mean(v)), 4) for k, v in d.items()}
                                                 for n, d in res.items()},
           "var_expl_train_top5_mean": np.mean([v[:5] for v in diag["var_expl_train"]], 0).round(4).tolist(),
           "heldout_pc_score_r2_library": np.mean(diag["lib_r2_pc"], 0).round(3).tolist(),
           "heldout_pc_score_r2_lineage": np.mean(diag["lineage_r2_pc"], 0).round(3).tolist(),
           "m1_r_predicted_vs_observed_pc_scores": np.mean(diag["m1_r_pc_scores"], 0).round(3).tolist()}
    (Path(__file__).parent / f"neg_lowrank_{cls_name}.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ROOT / "cache" / "runs",
         sys.argv[2] if len(sys.argv) > 2 else "nonessential")
