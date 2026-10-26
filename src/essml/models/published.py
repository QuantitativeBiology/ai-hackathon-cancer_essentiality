"""Published-method re-runs R2 (Shi-style elastic net) and R4 (TrPLet-like top-k), DESIGN.md §5.

Both are per-gene models on expression only (E -> ``FoldTransformer`` fitted on the training lines,
standardised, clipped at +-8), fitted on each gene's observed training rows (NaN never filled),
raw y with an intercept. Modelled genes = training-fold variable genes (SD > 0.25), optionally
intersected with a pre-declared subset; every other gene gets B0 (training gene mean). Inner splits
are patient-grouped and lineage-stratified (``g3_common.grouped_kfold``) over the training lines and
restricted to each gene's observed rows. Per-gene fits run in ``n_jobs`` processes x 1 thread.

R2 ``ElasticNetR2``: sklearn ``ElasticNetCV`` with l1_ratio = 0.5 (Shi et al. 2024: alpha = 0.5 in
glmnet notation), penalty path of ``n_alphas`` values down to ``eps`` x the largest (glmnet default
lambda.min.ratio = 0.01 for n < p), penalty chosen by inner grouped CV (MSE).

R4 ``TopKCorrR4`` (TrPLet-like): candidates = k in ``ks`` x learner in {ridge (RidgeCV, LOO over
``ridge_alphas``), RBF-SVR (C=1, epsilon=0.1, gamma='scale', on y standardised with the fitting
rows)} on the k features with the largest |Pearson r| with y on the fitting rows.
- nested=True: the candidate is chosen per gene by mean inner-CV r (``inner_folds`` grouped
  splits; the top-k screen is recomputed inside every inner training part); refit on all observed
  training rows.
- nested=False: no inner CV; all candidates are fitted on the training fold and their held-out
  predictions are kept (``cand_pred_``) for ``oracle_select`` (run.py, leaky=1 configs only), which
  picks per gene the candidate with the best mean outer-fold r = TrPLet's best-of-N selection
  inside the reported CV (deliberately optimistic diagnostic).
"""
from __future__ import annotations

import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from .g3_common import (CellFeatures, default_groups, grouped_kfold, make_pool, masked_r,
                        top_abs_corr, variable_mask)

ROOT = Path(__file__).resolve().parents[3]
_W: dict = {}  # worker globals (set once per worker by _init)


def _init(Xtr, Xho, fold_id):
    _W["Xtr"], _W["Xho"], _W["fold_id"] = Xtr, Xho, fold_id


def _rows_splits(rows, fold_id):
    f = fold_id[rows]
    return [(np.flatnonzero(f != k), np.flatnonzero(f == k)) for k in np.unique(f)]


# ============================================================================ R2 worker
def _enet_gene(task):
    from sklearn.linear_model import ElasticNetCV
    from sklearn.exceptions import ConvergenceWarning
    y = task["y"]
    rows = np.flatnonzero(~np.isnan(y))
    Xtr = _W["Xtr"]
    X = Xtr if len(rows) == len(Xtr) else Xtr[rows]
    yr = y[rows]
    m = ElasticNetCV(l1_ratio=0.5, n_alphas=task["n_alphas"], eps=task["eps"],
                     cv=_rows_splits(rows, _W["fold_id"]), max_iter=task["max_iter"], tol=task["tol"],
                     fit_intercept=True, precompute=False, selection="cyclic")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ConvergenceWarning)
        m.fit(X, yr)
    return _W["Xho"] @ m.coef_ + m.intercept_, float(m.alpha_), int((m.coef_ != 0).sum())


# ============================================================================ R4 worker
def _cand_fit_predict(Xa, ya, Xb, ks, top, ridge_alphas):
    """Predictions on Xb of every candidate (len(ks) x 2 columns: ridge then SVR per k)."""
    from sklearn.linear_model import RidgeCV
    from sklearn.svm import SVR
    out = np.empty((len(Xb), 2 * len(ks)))
    mu, sd = ya.mean(), ya.std()
    sd = sd if sd > 1e-12 else 1.0
    for i, k in enumerate(ks):
        cols = top[:k]
        A, B = Xa[:, cols], Xb[:, cols]
        out[:, 2 * i] = RidgeCV(alphas=ridge_alphas).fit(A, ya).predict(B)
        out[:, 2 * i + 1] = SVR(kernel="rbf", C=1.0, epsilon=0.1, gamma="scale").fit(
            A, (ya - mu) / sd).predict(B) * sd + mu
    return out


def _topk_gene(task):
    y = task["y"]
    ks, ra = task["ks"], task["ridge_alphas"]
    rows = np.flatnonzero(~np.isnan(y))
    Xtr = _W["Xtr"]
    kmax = max(ks)
    n_c = 2 * len(ks)
    if task["nested"]:
        sc = np.zeros((0, n_c))
        for a, b in _rows_splits(rows, _W["fold_id"]):
            ra_, rb_ = rows[a], rows[b]
            top = top_abs_corr(Xtr[ra_], y[ra_], kmax)
            P = _cand_fit_predict(Xtr[ra_], y[ra_], Xtr[rb_], ks, top, ra)
            r = masked_r(P, np.repeat(y[rb_][:, None], n_c, 1), min_obs=3)
            sc = np.vstack([sc, np.nan_to_num(r)])
        score = sc.mean(0)
        best = int(np.argmax(score))
    else:
        score, best = np.full(n_c, np.nan), 0
    top = top_abs_corr(Xtr[rows], y[rows], kmax)
    if task["store_all"]:
        P = _cand_fit_predict(Xtr[rows], y[rows], _W["Xho"], ks, top, ra)
        return P[:, best], best, score, P.astype(np.float32)
    k = ks[best // 2]
    ks1 = [k]
    P = _cand_fit_predict(Xtr[rows], y[rows], _W["Xho"], ks1, top, ra)
    return P[:, best % 2], best, score, None


# ============================================================================ shared base
class _PerGeneE:
    uses_groups = True
    worker = None

    def __init__(self, gene_subset=None, subset_first_n=None, inner_folds=5, split_seed=20261026,
                 n_jobs=6, min_rows=50, task_chunk=100, permute_seed=None):
        self.gene_subset, self.subset_first_n = gene_subset, subset_first_n
        self.inner_folds, self.split_seed = int(inner_folds), int(split_seed)
        self.n_jobs, self.min_rows, self.task_chunk = int(n_jobs), int(min_rows), int(task_chunk)
        self.permute_seed = permute_seed

    def _subset(self, genes):
        if self.gene_subset is None:
            return np.ones(len(genes), bool)
        p = Path(self.gene_subset)
        p = p if p.is_absolute() else ROOT / p
        sub = pd.read_csv(p, dtype={"entrez": str})["entrez"]
        if self.subset_first_n:
            sub = sub.iloc[: int(self.subset_first_n)]
        return pd.Index(map(str, genes)).isin(set(sub))

    def fit(self, Y: pd.DataFrame, lineage=None, X=None, groups=None):
        """Stores the training design; the per-gene fits run in ``predict`` (they need the new
        lines' design only to evaluate the fitted models: no new-line value enters a fit)."""
        t0 = time.perf_counter()
        ids = pd.Index([str(i) for i in Y.index])
        if groups is None:
            groups = default_groups(ids)
        lin = None if lineage is None else pd.Series(lineage).reindex(Y.index).astype(str)
        self.genes_ = Y.columns
        Yv = Y.to_numpy(np.float64)
        obs = ~np.isnan(Yv)
        with np.errstate(invalid="ignore"):
            self.mu_ = np.nan_to_num(np.nanmean(Yv, axis=0))
        model = variable_mask(Yv) & self._subset(Y.columns) & (obs.sum(0) >= self.min_rows)
        self.model_genes_ = np.flatnonzero(model)
        self.cf_, Z = CellFeatures(("E",), permute_seed=self.permute_seed).fit(ids, X)
        self.Ztr_ = Z.astype(np.float64)
        self.Ytr_ = Yv[:, self.model_genes_]
        fold_id = np.full(len(ids), -1)
        for k, (_, b) in enumerate(grouped_kfold(ids, np.asarray(groups), lin, seed=self.split_seed,
                                                 n_splits=self.inner_folds)):
            fold_id[b] = k
        self.fold_id_ = fold_id
        self.prep_s_ = time.perf_counter() - t0
        return self

    def _task(self, j):
        raise NotImplementedError

    def _collect(self, j, res):
        raise NotImplementedError

    def predict(self, index, lineage=None, X=None) -> pd.DataFrame:
        t0 = time.perf_counter()
        index = pd.Index(index)
        Zho = self.cf_.transform(pd.Index([str(i) for i in index]), X).astype(np.float64)
        P = np.broadcast_to(self.mu_, (len(index), len(self.mu_))).copy()
        J = len(self.model_genes_)
        self._start(len(index), J)
        with make_pool(self.n_jobs, _init, (self.Ztr_, Zho, self.fold_id_)) as ex:
            for s in range(0, J, self.task_chunk):
                js = list(range(s, min(s + self.task_chunk, J)))
                for j, res in zip(js, ex.map(self.worker, [self._task(j) for j in js], chunksize=2)):
                    P[:, self.model_genes_[j]] = res[0]
                    self._collect(j, res)
        self.report_ = dict(self._report(), n_model_genes=int(J), prep_s=round(self.prep_s_, 1),
                            genes_s=round(time.perf_counter() - t0, 1), n_features=int(Zho.shape[1]))
        self.ho_index_ = index
        self.Ztr_ = None  # free memory
        return pd.DataFrame(P, index=index, columns=self.genes_)


class ElasticNetR2(_PerGeneE):
    name = "R2"
    worker = staticmethod(_enet_gene)

    def __init__(self, n_alphas=40, eps=1e-2, max_iter=2000, tol=1e-3, **kw):
        super().__init__(**kw)
        self.n_alphas, self.eps, self.max_iter, self.tol = int(n_alphas), float(eps), int(max_iter), float(tol)

    def _task(self, j):
        return {"y": self.Ytr_[:, j], "n_alphas": self.n_alphas, "eps": self.eps,
                "max_iter": self.max_iter, "tol": self.tol}

    def _start(self, m, J):
        self.alpha_ = np.full(J, np.nan)
        self.n_nonzero_ = np.zeros(J, int)

    def _collect(self, j, res):
        self.alpha_[j], self.n_nonzero_[j] = res[1], res[2]

    def _report(self):
        nz = self.n_nonzero_
        return {"median_nonzero": float(np.median(nz)) if len(nz) else None,
                "frac_all_zero": round(float((nz == 0).mean()), 4) if len(nz) else None}


class TopKCorrR4(_PerGeneE):
    name = "R4"
    worker = staticmethod(_topk_gene)

    def __init__(self, ks=(5, 10, 20, 50, 100, 200), ridge_alphas=(0.1, 1.0, 10.0, 100.0, 1000.0),
                 nested=True, store_all=False, inner_folds=3, **kw):
        super().__init__(inner_folds=inner_folds, **kw)
        self.ks = [int(k) for k in ks]
        self.ridge_alphas = tuple(float(a) for a in ridge_alphas)
        self.nested, self.store_all = bool(nested), bool(store_all) or not bool(nested)

    def candidate_names(self):
        return [f"{l}_k{k}" for k in self.ks for l in ("ridge", "svr")]

    def _task(self, j):
        return {"y": self.Ytr_[:, j], "ks": self.ks, "ridge_alphas": self.ridge_alphas,
                "nested": self.nested, "store_all": self.store_all}

    def _start(self, m, J):
        self.choice_ = np.zeros(J, int)
        self.inner_score_ = np.full((J, 2 * len(self.ks)), np.nan)
        self.cand_pred_ = np.zeros((m, J, 2 * len(self.ks)), np.float32) if self.store_all else None

    def _collect(self, j, res):
        self.choice_[j] = res[1]
        self.inner_score_[j] = res[2]
        if self.store_all:
            self.cand_pred_[:, j, :] = res[3]

    def _report(self):
        names = self.candidate_names()
        c = pd.Series(self.choice_).map(dict(enumerate(names))).value_counts()
        return {"nested": self.nested, "choice_counts": {str(k): int(v) for k, v in c.items()}}

    @staticmethod
    def oracle_select(res, Y, folds):
        """LEAKY (run.py, leaky=1 only): per gene, the candidate with the best mean outer-fold r.

        Uses the held-out targets of every outer fold. Returns the replaced prediction frame."""
        pred = res["pred"].copy()
        per_gene: dict = {}
        models = res["models"]
        for f, m in models.items():
            ho = m.ho_index_
            T = Y.loc[ho].to_numpy(np.float64)
            for j, g in enumerate(m.model_genes_):
                C = m.cand_pred_[:, j, :].astype(np.float64)
                t = np.repeat(T[:, g][:, None], C.shape[1], 1)
                per_gene.setdefault(m.genes_[g], []).append(masked_r(C, t, min_obs=10))
        best = {g: int(np.nanargmax(np.nanmean(np.vstack(v), 0))) for g, v in per_gene.items()
                if np.isfinite(np.vstack(v)).any()}
        for f, m in models.items():
            ho = m.ho_index_
            for j, g in enumerate(m.model_genes_):
                name = m.genes_[g]
                if name in best:
                    pred.loc[ho, name] = m.cand_pred_[:, j, best[name]].astype(np.float64)
        res.setdefault("model_report", {})["oracle_choice_counts"] = (
            pd.Series(best).value_counts().sort_index().to_dict())
        return pred
