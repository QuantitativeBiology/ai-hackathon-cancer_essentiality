"""M3: LightGBM on the out-of-fold residuals of m1_E_resB0 (DESIGN.md §5, RESEARCH_TREE D18).

Per outer training fold (training lines only reach any statistic):
1. Base. Out-of-fold M1 (B0 residual, E only = ``m1_E_resB0``) predictions for every training line
   from ``inner_folds`` patient-grouped, lineage-stratified folds of the training lines
   (``m2_local.LocalStageM2.inner_oof``: same splits and on-disk cache as the M2 track). M1 refitted
   on all training lines predicts new lines.
2. Residual target r = y - m1_oof for every training-fold variable gene (SD > 0.25 on the training
   fold; optionally intersected with a pre-declared gene subset); other genes pass M1 through.
3. Tree features per gene: the gene's own M1 prediction (OOF on training lines, full-fit on new
   lines) + the top-``k_screen`` columns by |Pearson r| with r_g from a pool of raw blocks
   (``FoldTransformer`` fitted on the training lines). Screening uses only the rows the trees are
   fitted on: inside the early-stopping split it is computed on the inner-training rows; for the
   final fit on all training rows.
4. Early stopping: ONE patient-grouped 80/20 split of the training lines; LightGBM (L2) with
   early stopping on the inner-validation rows gives ``best_iter``. If the best validation MSE is
   not below that of a zero correction, the gene gets no correction (best_iter = 0). Final model:
   ``best_iter`` rounds on all observed training rows.
5. Prediction = M1 full-fit prediction + correction (boosters stored as LightGBM model strings).

Per-gene fits run in ``n_jobs`` processes x 1 thread (LightGBM num_threads = 1, deterministic).
Negative control: ``permute_seed`` permutes the pool's rows across lines and ``use_m1_feature``
False drops the M1 feature, so the trees see no information about the line.
"""
from __future__ import annotations

import time

import numpy as np
import pandas as pd

from .g3_common import CellFeatures, default_groups, inner_split, make_pool, variable_mask

LGB_DEFAULT = {"objective": "regression", "learning_rate": 0.05, "num_leaves": 15,
               "min_data_in_leaf": 20, "feature_fraction": 1.0, "bagging_fraction": 1.0,
               "bagging_freq": 0, "lambda_l2": 1.0, "max_bin": 63, "num_threads": 1,
               "deterministic": True, "force_row_wise": True, "verbose": -1, "seed": 0}


def screen_top_k(Z, R, k, rows_list):
    """Top-k |r| pool columns per gene.

    Z: (n, p) pool; R: (n, G) residuals with NaN; rows_list: per gene, observed row positions
    (boolean masks are grouped so genes sharing a row set share one standardisation).
    Returns (G, k) int array of column indices."""
    G = R.shape[1]
    out = np.zeros((G, k), dtype=np.int64)
    keys = {}
    for g, rows in enumerate(rows_list):
        keys.setdefault(rows.tobytes(), (rows, []))[1].append(g)
    for rows, genes in keys.values():
        idx = np.flatnonzero(rows)
        Zr = Z[idx].astype(np.float64)
        Zr -= Zr.mean(0)
        zs = np.sqrt((Zr ** 2).sum(0))
        zs[zs == 0] = np.inf
        Zr /= zs
        for s in range(0, len(genes), 256):
            gs = genes[s:s + 256]
            Rr = R[np.ix_(idx, gs)]
            Rr = Rr - Rr.mean(0)
            rs = np.sqrt((Rr ** 2).sum(0))
            rs[rs == 0] = np.inf
            C = np.abs(Zr.T @ (Rr / rs))  # (p, len(gs))
            kk = min(k, C.shape[0])
            top = np.argpartition(-C, kk - 1, axis=0)[:kk]
            for j, g in enumerate(gs):
                t = top[:, j]
                out[g, :kk] = t[np.argsort(-C[t, j], kind="stable")]
    return out


def _fit_gene(task):
    """Worker: early-stopped LightGBM for one gene. Returns (model string or None, best_iter, mse)."""
    import lightgbm as lgb
    params, n_round, es = task["params"], task["n_round"], task["es"]
    Xa, ya, Xb, yb = task["Xa"], task["ya"], task["Xb"], task["yb"]
    dtr = lgb.Dataset(Xa, ya, free_raw_data=True)
    dva = lgb.Dataset(Xb, yb, reference=dtr)
    bst = lgb.train(params, dtr, num_boost_round=n_round, valid_sets=[dva],
                    callbacks=[lgb.early_stopping(es, verbose=False)])
    best = int(bst.best_iteration or 0)
    mse_best = float(np.mean((bst.predict(Xb, num_iteration=best) - yb) ** 2)) if best else np.inf
    mse_zero = float(np.mean(yb ** 2))
    if not mse_best < mse_zero:
        return None, 0, (mse_best, mse_zero)
    bst = lgb.train(params, lgb.Dataset(task["Xf"], task["yf"]), num_boost_round=best)
    return bst.model_to_string(), best, (mse_best, mse_zero)


class LGBMResidualM3:
    name = "M3"
    uses_groups = True

    def __init__(self, m1_params=None, inner_folds=5, oof_seed=20261026, split_seed=20261026,
                 pool_blocks=("E", "C", "M", "H", "meta"), k_screen=50, use_m1_feature=True,
                 lgb_params=None, num_boost_round=500, early_stopping=30, n_jobs=6,
                 gene_subset=None, subset_first_n=None, min_rows=50, permute_seed=None,
                 task_chunk=200, cache_inner=True):
        self.m1_params = dict(m1_params or {"modalities": ["E"], "residual": "B0"})
        self.inner_folds, self.oof_seed, self.split_seed = int(inner_folds), int(oof_seed), int(split_seed)
        self.pool_blocks = list(pool_blocks)
        self.k_screen = int(k_screen)
        self.use_m1_feature = bool(use_m1_feature)
        self.lgb_params = dict(LGB_DEFAULT, **(lgb_params or {}))
        self.num_boost_round, self.early_stopping = int(num_boost_round), int(early_stopping)
        self.n_jobs = int(n_jobs)
        self.gene_subset, self.subset_first_n = gene_subset, subset_first_n
        self.min_rows = int(min_rows)
        self.permute_seed = permute_seed
        self.task_chunk = int(task_chunk)
        self.cache_inner = bool(cache_inner)

    # -------------------------------------------------------------- helpers
    def _subset(self, genes):
        if self.gene_subset is None:
            return np.ones(len(genes), bool)
        from pathlib import Path
        p = Path(self.gene_subset)
        if not p.is_absolute():
            p = Path(__file__).resolve().parents[3] / p
        sub = pd.read_csv(p, dtype={"entrez": str})["entrez"]
        if self.subset_first_n:
            sub = sub.iloc[: int(self.subset_first_n)]
        return pd.Index(map(str, genes)).isin(set(sub))

    def _base(self, Y, lin, groups, X_m1):
        from .kernel_ridge import KernelRidgeM1
        from .m2_local import LocalStageM2
        st = LocalStageM2(m1_params=self.m1_params, inner_folds=self.inner_folds,
                          seed=self.oof_seed, cache_inner=self.cache_inner)
        oof, _, cached = st.inner_oof(Y, lin, groups, X_m1)
        oof = np.asarray(oof, np.float64)
        bad = ~np.isfinite(oof)
        if bad.any():
            mu = np.nan_to_num(np.nanmean(Y.to_numpy(np.float64), axis=0))
            oof = np.where(bad, mu[None, :], oof)
        m1 = KernelRidgeM1(**self.m1_params).fit(Y, lineage=lin, X=X_m1, groups=groups)
        return oof, m1, cached

    # -------------------------------------------------------------- fit / predict
    def fit(self, Y: pd.DataFrame, lineage=None, X=None, groups=None):
        t0 = time.perf_counter()
        ids = pd.Index([str(i) for i in Y.index])
        if groups is None:
            groups = default_groups(ids)
        groups = np.asarray(groups)
        lin = (pd.Series("x", index=Y.index) if lineage is None
               else pd.Series(lineage).reindex(Y.index).astype(str))
        frame = isinstance(X, pd.DataFrame)
        X_m1 = X if frame else None
        self.genes_ = Y.columns
        oof, self.m1_, cached = self._base(Y, lin, groups, X_m1)
        t_base = time.perf_counter() - t0

        Yv = Y.to_numpy(np.float64)
        model = variable_mask(Yv) & self._subset(Y.columns)
        obs = ~np.isnan(Yv)
        model &= obs.sum(0) >= self.min_rows
        self.model_genes_ = np.flatnonzero(model)
        R = Yv[:, self.model_genes_] - oof[:, self.model_genes_]
        O = obs[:, self.model_genes_]
        self.cf_, Zp = CellFeatures(self.pool_blocks if not frame else ("X",),
                                    permute_seed=self.permute_seed).fit(ids, X)
        a, b = inner_split(ids, groups, lin, seed=self.split_seed)
        in_a = np.zeros(len(ids), bool)
        in_a[a] = True
        t1 = time.perf_counter()
        k = min(self.k_screen, Zp.shape[1])
        sel_in = screen_top_k(Zp, np.where(O & in_a[:, None], R, np.nan), k,
                              [O[:, j] & in_a for j in range(O.shape[1])])
        sel_full = screen_top_k(Zp, R, k, [O[:, j] for j in range(O.shape[1])])
        t_screen = time.perf_counter() - t1
        self.sel_full_ = sel_full
        cache = (Zp, oof, R, O, a, b)
        J = len(self.model_genes_)
        self.boosters_ = [None] * J
        self.best_iter_ = np.zeros(J, int)
        mse = np.full((J, 2), np.nan)
        t2 = time.perf_counter()
        with make_pool(self.n_jobs) as ex:
            for s0 in range(0, J, self.task_chunk):
                js = range(s0, min(s0 + self.task_chunk, J))
                for j, (m, it, e) in zip(js, ex.map(_fit_gene, self._tasks(cache, sel_in, js),
                                                     chunksize=8)):
                    self.boosters_[j], self.best_iter_[j], mse[j] = m, it, e
        bi = self.best_iter_
        ok = np.isfinite(mse).all(1) & (bi > 0)
        self.report_ = {"n_model_genes": int(J), "pool_cols": int(Zp.shape[1]), "k_screen": int(k),
                        "inner_cached": bool(cached), "base_s": round(t_base, 1),
                        "screen_s": round(t_screen, 1), "trees_s": round(time.perf_counter() - t2, 1),
                        "n_inner_train": int(len(a)), "n_inner_val": int(len(b)),
                        "frac_genes_no_correction": round(float((bi == 0).mean()), 4) if J else None,
                        "median_best_iter": float(np.median(bi)) if J else None,
                        "median_val_mse_ratio": (round(float(np.median(mse[ok, 0] / mse[ok, 1])), 4)
                                                 if ok.any() else None),
                        "fit_s": round(time.perf_counter() - t0, 1)}
        return self

    @staticmethod
    def _feats(use_m1, Zrows, m1col, cols):
        F = Zrows[:, cols]
        if use_m1:
            F = np.column_stack([m1col, F])
        return np.ascontiguousarray(F, dtype=np.float32)

    def _tasks(self, cache, sel_in, js):
        Zp, oof, R, O, a, b = cache
        f = lambda Zr, m1c, cols: self._feats(self.use_m1_feature, Zr, m1c, cols)  # noqa: E731
        for j in js:
            g = self.model_genes_[j]
            oa = a[O[a, j]]
            ob = b[O[b, j]]
            of = np.flatnonzero(O[:, j])
            yield {"params": self.lgb_params, "n_round": self.num_boost_round, "es": self.early_stopping,
                   "Xa": f(Zp[oa], oof[oa, g], sel_in[j]), "ya": R[oa, j],
                   "Xb": f(Zp[ob], oof[ob, g], sel_in[j]), "yb": R[ob, j],
                   "Xf": f(Zp[of], oof[of, g], self.sel_full_[j]), "yf": R[of, j]}

    def predict(self, index, lineage=None, X=None) -> pd.DataFrame:
        import lightgbm as lgb
        index = pd.Index(index)
        sid = pd.Index([str(i) for i in index])
        frame = isinstance(X, pd.DataFrame)
        lin = None if lineage is None else pd.Series(lineage).reindex(index)
        m1new = self.m1_.predict(index, lineage=lin, X=X if frame else None).to_numpy(np.float64)
        Znew = self.cf_.transform(sid, X)
        P = m1new.copy()
        for j, m in enumerate(self.boosters_):
            if m is None:
                continue
            g = self.model_genes_[j]
            bst = lgb.Booster(model_str=m)
            P[:, g] += bst.predict(self._feats(self.use_m1_feature, Znew, m1new[:, g],
                                               self.sel_full_[j]), num_threads=1)
        return pd.DataFrame(P, index=index, columns=self.genes_)
