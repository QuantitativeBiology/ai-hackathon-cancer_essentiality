"""B2 vectorised predict (e44a4c8) vs the old per-query loop (987bdec), on random data with NaNs,
all-NaN neighbourhoods, k = 1, k = n_train and k > n_train.

    python reviews/g2_review/check_b2.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from essml.baselines import KNNExpression, _nanmean_cols  # noqa: E402


def old_predict(m, index, X, k=None):
    nb = m.neighbours(X.loc[pd.Index(index)], k)
    P = np.empty((len(nb), m.Y_.shape[1]))
    for r, idx in enumerate(nb):
        mu, n = _nanmean_cols(m.Y_[idx])
        P[r] = np.where(n > 0, mu, m.mu_)
    return pd.DataFrame(P, index=pd.Index(index), columns=m.genes_)


rng = np.random.default_rng(0)
n, nq, p, G = 120, 40, 300, 60
ids = [f"L{i}" for i in range(n + nq)]
X = pd.DataFrame(rng.normal(size=(n + nq, p)), index=ids)
Y = pd.DataFrame(rng.normal(size=(n, G)), index=ids[:n])
Y = Y.mask(rng.random(Y.shape) < 0.3)
Y.iloc[:, 0] = np.nan          # gene unobserved everywhere -> mu (NaN) both ways
Y.iloc[5:, 1] = np.nan         # gene observed in 5 lines only -> many all-NaN neighbourhoods
m = KNNExpression(k=7, n_top=100).fit(Y, X=X)
worst = 0.0
for k in [1, 7, 20, n, n + 50]:
    a = m.predict(ids[n:], X=X, k=k).to_numpy()
    b = old_predict(m, ids[n:], X, k).to_numpy()
    assert np.array_equal(np.isnan(a), np.isnan(b)), k
    worst = max(worst, float(np.nanmax(np.abs(a - b))))
print(f"max |vectorised - loop| = {worst:.2e}  ->", "PASS" if worst < 1e-12 else "FAIL")
