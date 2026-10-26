"""Baselines B0/B1/B2 (DESIGN.md §5). NaN-aware throughout: CRISPR NaNs are never filled.

Common interface (used by scripts/run.py):
    model.fit(Y, lineage=None, X=None)       Y: DataFrame [train ModelID x gene], may hold NaN
    model.predict(index, lineage=None, X=None) -> DataFrame [index x gene], finite wherever the
                                                  gene has >=1 observed training value
Training statistics use observed entries only (per-gene counts differ under the NaN patterns).

- B0 GeneMean        per-gene mean over observed training lines.
- B1 LineageMean     per-gene lineage mean shrunk to the gene mean with pseudo-count kappa:
                     (sum_lg + kappa * mu_g) / (n_lg + kappa); unseen lineage (n_lg + kappa == 0)
                     -> mu_g (B0).
- B2 KNNExpression   mean CRISPR profile of the k training lines most correlated with the query
                     in expression space (top-N variance genes on training lines, z-scored with
                     training mean/SD; Pearson across genes). Per gene, neighbours with NaN are
                     skipped; if all k are NaN -> mu_g.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def _nanmean_cols(Y: np.ndarray):
    obs = ~np.isnan(Y)
    n = obs.sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(obs, Y, 0.0).sum(0) / n, n


def _onehot(labels, levels):
    pos = {l: i for i, l in enumerate(levels)}
    L = np.zeros((len(labels), len(levels)))
    for r, l in enumerate(labels):
        L[r, pos[l]] = 1.0
    return L


class GeneMean:
    """B0."""
    name = "B0"

    def fit(self, Y: pd.DataFrame, lineage=None, X=None):
        self.genes_ = Y.columns
        self.mu_, self.n_ = _nanmean_cols(Y.to_numpy(np.float64))
        return self

    def predict(self, index, lineage=None, X=None) -> pd.DataFrame:
        P = np.broadcast_to(self.mu_, (len(index), len(self.mu_)))
        return pd.DataFrame(np.array(P), index=pd.Index(index), columns=self.genes_)


class LineageMean:
    """B1: shrunk lineage mean with pseudo-count kappa (kappa=0 -> raw lineage mean)."""
    name = "B1"

    def __init__(self, kappa: float = 10.0):
        self.kappa = float(kappa)

    def fit(self, Y: pd.DataFrame, lineage: pd.Series, X=None):
        lin = pd.Series(lineage).reindex(Y.index)
        if lin.isna().any():
            raise ValueError("lineage missing for some training lines")
        Yv = Y.to_numpy(np.float64)
        obs = ~np.isnan(Yv)
        self.genes_ = Y.columns
        self.mu_, _ = _nanmean_cols(Yv)
        self.levels_ = sorted(lin.unique())
        L = _onehot(lin.to_numpy(), self.levels_)
        self.sum_ = L.T @ np.where(obs, Yv, 0.0)  # (n_lineage, n_gene) masked sums
        self.cnt_ = L.T @ obs.astype(np.float64)   # observed counts per lineage x gene
        self.means_ = self._shrunk(self.sum_, self.cnt_, self.mu_)
        return self

    def _shrunk(self, S, N, mu):
        den = N + self.kappa
        with np.errstate(invalid="ignore", divide="ignore"):
            M = (S + self.kappa * mu) / den
        return np.where(den > 0, M, mu)

    def lineage_means(self) -> pd.DataFrame:
        """Shrunk per-lineage means (lineage x gene)."""
        return pd.DataFrame(self.means_, index=self.levels_, columns=self.genes_)

    def predict(self, index, lineage: pd.Series, X=None) -> pd.DataFrame:
        lin = pd.Series(lineage).reindex(pd.Index(index))
        pos = {l: i for i, l in enumerate(self.levels_)}
        P = np.empty((len(lin), len(self.genes_)))
        for r, l in enumerate(lin.to_numpy()):
            P[r] = self.means_[pos[l]] if l in pos else self.mu_  # unseen lineage -> B0
        return pd.DataFrame(P, index=pd.Index(index), columns=self.genes_)


def leave_line_out_lineage_residual(Y: pd.DataFrame, lineage: pd.Series, kappa: float) -> pd.DataFrame:
    """Training residuals y_ig - B1_{-i}(i, g), where B1_{-i} is B1 fitted without line i.

    Exact closed form of the brute-force loop (tested): removing line i changes the gene mean,
    the lineage sum and the lineage count only for genes observed in line i. NaN where y is NaN.
    """
    lin = pd.Series(lineage).reindex(Y.index)
    Yv = Y.to_numpy(np.float64)
    obs = ~np.isnan(Yv)
    Y0 = np.where(obs, Yv, 0.0)
    S, N = Y0.sum(0), obs.sum(0).astype(np.float64)
    levels = sorted(lin.unique())
    L = _onehot(lin.to_numpy(), levels)
    Sl, Nl = L.T @ Y0, L.T @ obs.astype(np.float64)
    o = obs.astype(np.float64)
    Si, Ni = S - Y0, N - o                       # gene totals without line i
    Sli, Nli = L @ Sl - Y0, L @ Nl - o           # own-lineage totals without line i
    with np.errstate(invalid="ignore", divide="ignore"):
        mu_i = Si / Ni                           # gene mean without line i
        den = Nli + kappa
        pred = (Sli + kappa * mu_i) / den
    pred = np.where(den > 0, pred, mu_i)         # empty lineage after removal -> B0
    R = np.where(obs, Yv - pred, np.nan)
    return pd.DataFrame(R, index=Y.index, columns=Y.columns)


class KNNExpression:
    """B2: kNN in expression space."""
    name = "B2"
    predict_params = ("k",)  # k can vary at predict time without refitting (inner CV)

    def __init__(self, k: int = 20, n_top: int = 2000):
        self.k, self.n_top = int(k), int(n_top)

    def _embed(self, X: pd.DataFrame):
        Z = (X[self.feat_].to_numpy(np.float64) - self.xm_) / self.xs_
        Z = Z - Z.mean(1, keepdims=True)
        return Z / np.linalg.norm(Z, axis=1, keepdims=True)

    def fit(self, Y: pd.DataFrame, lineage=None, X: pd.DataFrame = None):
        Xt = X.loc[Y.index]
        if Xt.isna().any().any():
            raise ValueError("expression features contain NaN")
        v = Xt.var(ddof=1)
        self.feat_ = v.sort_values(ascending=False, kind="stable").index[: self.n_top]
        self.xm_ = Xt[self.feat_].mean().to_numpy()
        self.xs_ = Xt[self.feat_].std(ddof=1).to_numpy()
        self.Zt_ = self._embed(Xt)
        self.Y_ = Y.to_numpy(np.float64)
        self.genes_ = Y.columns
        self.train_ids_ = Y.index
        self.mu_, _ = _nanmean_cols(self.Y_)
        return self

    def neighbours(self, X_new: pd.DataFrame, k=None):
        """Indices (into training lines) of the k most correlated training lines per query."""
        k = min(k or self.k, len(self.train_ids_))
        C = self._embed(X_new) @ self.Zt_.T
        return np.argsort(-C, axis=1, kind="stable")[:, :k]

    def predict(self, index, lineage=None, X: pd.DataFrame = None, k=None) -> pd.DataFrame:
        nb = self.neighbours(X.loc[pd.Index(index)], k)
        P = np.empty((len(nb), self.Y_.shape[1]))
        for r, idx in enumerate(nb):
            m, n = _nanmean_cols(self.Y_[idx])
            P[r] = np.where(n > 0, m, self.mu_)
        return pd.DataFrame(P, index=pd.Index(index), columns=self.genes_)


MODELS = {"B0": GeneMean, "B1": LineageMean, "B2": KNNExpression}
