"""Feasibility probe for DESIGN.md: dev CV fold 0 only (test set untouched).

Measures per-gene Pearson (held-out lines) for gene mean, lineage mean and ridge on expression residuals,
plus wall-clock. Not a result to report; it calibrates the design.
"""
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from essml.io import load_matrix  # noqa: E402

s = pd.read_csv("configs/splits.csv", index_col=0)
dev = s[s.role == "dev"]
tr, va = dev.index[dev.cv_fold != 0], dev.index[dev.cv_fold == 0]
Y = load_matrix("crispr")
X = load_matrix("expression")
Y = Y.loc[dev.index]
Y = Y.loc[:, Y.isna().mean() < 0.1]
Y = Y.fillna(Y.mean())
Ytr, Yva = Y.loc[tr].values, Y.loc[va].values
mu = Ytr.mean(0)


def pergene_r(P, T):
    P = P - P.mean(0); T = T - T.mean(0)
    den = np.sqrt((P ** 2).sum(0) * (T ** 2).sum(0))
    with np.errstate(invalid="ignore", divide="ignore"):
        r = (P * T).sum(0) / den
    return np.nan_to_num(r)  # constant prediction -> r = 0


sd = Ytr.std(0)
sel = sd > np.quantile(sd, 0.9)  # top-10% variable genes
lin = dev["OncotreeLineage"].fillna("Other")
t0 = time.time()
L = pd.DataFrame(Ytr, index=tr).groupby(lin.loc[tr].values).mean()
n = lin.loc[tr].value_counts()
k = 10  # shrinkage pseudo-count toward gene mean
Lsh = (L.mul(n.loc[L.index], axis=0) + k * mu) / (n.loc[L.index].values[:, None] + k)
Plin = np.vstack([Lsh.loc[l].values if l in Lsh.index else mu for l in lin.loc[va]])
t_lin = time.time() - t0
r_lin = pergene_r(Plin, Yva)

Xd = X.loc[dev.index]
Xd = Xd.loc[:, Xd.loc[tr].std() > 0.5]
xm, xs = Xd.loc[tr].mean(), Xd.loc[tr].std()
Xtr = ((Xd.loc[tr] - xm) / xs).values
Xva = ((Xd.loc[va] - xm) / xs).values
R = Ytr - Lsh.reindex(lin.loc[tr]).values  # residual from lineage mean
res = {}
for lam in [1e3, 1e4, 3e4, 1e5]:
    t0 = time.time()
    K = Xtr @ Xtr.T
    A = np.linalg.solve(K + lam * np.eye(len(K)), R)
    P = Plin + (Xva @ Xtr.T) @ A
    res[lam] = (pergene_r(P, Yva), time.time() - t0)
print(f"genes={Y.shape[1]} train={len(tr)} val={len(va)} expr_feats={Xd.shape[1]}")
print(f"lineage-mean: mean r all={r_lin.mean():.3f} top10%var={r_lin[sel].mean():.3f} time={t_lin:.1f}s")
for lam, (r, t) in res.items():
    print(f"ridge lam={lam:g}: mean r all={r.mean():.3f} top10%var={r[sel].mean():.3f} n(r>0.2)={(r>0.2).sum()} time={t:.1f}s")
