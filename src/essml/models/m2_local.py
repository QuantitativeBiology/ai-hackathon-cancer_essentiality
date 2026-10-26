"""M2: per-gene local stage on out-of-fold M1 predictions (DESIGN.md §5; R1 algorithms review #5).

Per outer training fold:
1. Inner grouped folds (StratifiedGroupKFold by PatientID, stratified by lineage) over the
   TRAINING lines only. For each inner fold an M1/M1b model is fitted on the other inner folds
   and predicts the inner held-out lines -> out-of-fold M1 predictions m_oof for every training
   line. (Cached on disk by a hash of the training IDs + M1 config; shared by prior ablations.)
2. The M1 model refitted on all training lines gives the predictions used for new lines.
3. For every gene g a ridge with an unpenalised intercept on its observed training rows:
       y_g ~ m_oof_g + local features of g (own gene / paralogs / partners; priors.py)
   Features are standardised per gene on the observed training rows; normal equations are
   batched over genes ((G, p, p) solves, gene chunks).
4. lambda (shared within gene groups variable/other, classes from the training fold) is chosen on
   the same inner folds: fit the per-gene ridge on the other inner folds' rows, score mean per-gene
   r on the inner held-out rows.
   ``nested_lambda=True`` (G3; review finding 15): in step 4 the training rows a_s of inner split s
   get M1 features from a NESTED out-of-fold pass inside a_s only (inner-inner grouped folds), so
   no feature used to fit the local ridge was produced by an M1 model that saw the validation rows
   b_s. The validation rows b_s keep m_oof (M1 fitted on a_s), as in the final fit.
Genes with < ``min_rows`` observed training rows pass the M1 prediction through.
``local_genes='variable'`` (G3): the local stage is fitted only for training-fold variable genes;
every other gene passes the M1 prediction through (G2 showed the per-gene stage overfits
low-variance genes).

Prior-knowledge ablation (``priors``): 'none' (own gene), 'paralogs' (+ paralog features),
'partners' (+ partner mean expression).
"""
from __future__ import annotations

import hashlib
import json
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from essml import genesets
from essml import priors as PR
from .kernel_ridge import KernelRidgeM1, KernelRidgeT, _colwise_r

ROOT = Path(__file__).resolve().parents[3]
INNER_CACHE = ROOT / "cache" / "m2_inner"
OWN = ["own_expr", "not_expressed", "expr_hinge", "own_cn", "deep_del", "own_mut"]
PARALOG = ["par1_expr", "par2_expr", "par3_expr", "par_min_expr", "par_deep_del", "par_mut"]
PARTNER = ["partner_mean_expr"]
PRIOR_SETS = {"none": OWN, "paralogs": OWN + PARALOG, "partners": OWN + PARALOG + PARTNER}
LOCAL_LAM = np.logspace(0, 4, 9)
_BUILDER = None


def split_X(X):
    """X for M2: None (cache), or {"_m1": X for the M1 model, "_local": LocalFeatureBuilder}."""
    if isinstance(X, dict) and ("_m1" in X or "_local" in X):
        return X.get("_m1"), X.get("_local") or _builder()
    return X, None


def _builder():
    global _BUILDER
    if _BUILDER is None:
        _BUILDER = PR.LocalFeatureBuilder.from_cache()
    return _BUILDER


def _stack(feat: dict, names, m1: np.ndarray):
    """(p, n, G) float64: M1 prediction first, then local features."""
    arr = [m1] + [feat[k] for k in names]
    return np.stack([np.asarray(a, np.float64) for a in arr])


def normal_equations(Xs, Y, min_rows):
    """Standardised per-gene normal equations on observed rows (all genes of a chunk at once).

    Xs: (p, n, G); Y: (n, G) with NaN. Returns a dict for ``solve_ridge``.
    """
    p, n, G = Xs.shape
    w = (~np.isnan(Y)).astype(np.float64)
    nw = w.sum(0)
    y0 = np.where(w > 0, Y, 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        my = np.nan_to_num(y0.sum(0) / nw)
        mx = np.nan_to_num((Xs * w).sum(1) / nw)          # (p, G)
    Xc = ((Xs - mx[:, None, :]) * w).transpose(2, 0, 1)  # (G, p, n), zero on unobserved rows
    yc = ((y0 - my) * w).T[:, :, None]                    # (G, n, 1)
    A = Xc @ Xc.transpose(0, 2, 1)
    b = (Xc @ yc)[..., 0]
    sd = np.sqrt(np.maximum(np.einsum("gii->gi", A), 0) / np.maximum(nw, 1)[:, None])
    sd = np.where(sd > 1e-8, sd, np.inf)                  # constant features -> coefficient 0
    return {"As": A / (sd[:, :, None] * sd[:, None, :]), "bs": b / sd, "sd": sd, "mx": mx,
            "my": my, "ok": nw >= min_rows, "p": p}


def solve_ridge(ne, lams):
    """(beta (G, p) on the original feature scale, intercept (G,), ok) for lambda per gene."""
    lams = np.broadcast_to(np.asarray(lams, np.float64), ne["my"].shape)
    As = ne["As"] + (lams[:, None, None] + 1e-10) * np.eye(ne["p"])[None]
    beta = np.linalg.solve(As, ne["bs"][..., None])[..., 0] / ne["sd"]
    icpt = ne["my"] - np.einsum("gi,ig->g", beta, ne["mx"])
    return beta, icpt, ne["ok"]


def batched_ridge(Xs, Y, lams, min_rows):
    """Per-gene ridge with an unpenalised intercept on observed rows (see normal_equations)."""
    return solve_ridge(normal_equations(Xs, Y, min_rows), lams)


def apply_ridge(Xs, beta, icpt):
    return icpt[None, :] + np.einsum("ing,gi->ng", Xs, beta)


class LocalStageM2:
    name = "M2"
    uses_groups = True

    def __init__(self, m1_params=None, base="M1", priors="partners", inner_folds=5,
                 lam_grid=None, min_rows=50, chunk=1500, seed=20261026, cache_inner=True,
                 local_genes="all", nested_lambda=False):
        self.m1_params = dict(m1_params or {"modalities": ["E"]})
        self.base = base
        if priors not in PRIOR_SETS:
            raise ValueError(f"priors must be one of {list(PRIOR_SETS)}")
        self.priors = priors
        self.inner_folds = int(inner_folds)
        self.lam_grid = np.asarray(LOCAL_LAM if lam_grid is None else lam_grid, np.float64)
        self.min_rows, self.chunk, self.seed = int(min_rows), int(chunk), int(seed)
        self.cache_inner = bool(cache_inner)
        if local_genes not in ("all", "variable"):
            raise ValueError("local_genes must be 'all' or 'variable'")
        self.local_genes = local_genes
        self.nested_lambda = bool(nested_lambda)

    def _m1(self):
        cls = KernelRidgeT if self.base == "T" else KernelRidgeM1
        return cls(**self.m1_params)

    def _inner_splits(self, ids, lineage, groups):
        from sklearn.model_selection import StratifiedGroupKFold
        sgk = StratifiedGroupKFold(n_splits=self.inner_folds, shuffle=True, random_state=self.seed)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            return list(sgk.split(np.zeros(len(ids)), lineage.loc[ids].to_numpy(), np.asarray(groups)))

    def _cache_key(self, ids):
        h = hashlib.sha256(json.dumps({"ids": sorted(map(str, ids)), "base": self.base,
                                       "m1": self.m1_params, "k": self.inner_folds,
                                       "seed": self.seed}, sort_keys=True, default=str).encode())
        return h.hexdigest()[:24]

    def inner_oof(self, Y, lineage, groups, X=None):
        """Out-of-fold M1 predictions for every training line (inner folds of the training set)."""
        ids = Y.index
        splits = self._inner_splits(ids, lineage, groups)
        cacheable = self.cache_inner and X is None
        path = INNER_CACHE / f"{self._cache_key(ids)}.npy"
        if cacheable and path.exists():
            return np.load(path), splits, True
        oof = np.full(Y.shape, np.nan)
        grp = np.asarray(groups)
        sub = (lambda ix: None) if X is None else (lambda ix: X.loc[ids[ix]])
        for a, b in splits:
            m = self._m1().fit(Y.iloc[a], lineage=lineage, X=sub(a), groups=grp[a])
            oof[b] = m.predict(ids[b], lineage=lineage, X=sub(b)).to_numpy()
        if cacheable:
            INNER_CACHE.mkdir(parents=True, exist_ok=True)
            np.save(path, oof.astype(np.float32))
            oof = np.load(path)  # identical values whether cached or fresh
        return oof, splits, False

    def _local(self, ids, builder=None, gs=slice(None)):
        return (builder or _builder()).local_features(list(map(str, ids)), list(self.genes_[gs]))

    @staticmethod
    def _fill(oof, Yv):
        """Genes unobserved in an inner training part (NaN OOF): use the training gene mean."""
        oof = np.asarray(oof, np.float64)
        bad = ~np.isfinite(oof)
        if bad.any():
            mu = np.nan_to_num(np.nanmean(Yv, axis=0))
            oof = np.where(bad, mu[None, :], oof)
        return oof, int(bad.sum())

    def lambda_cv_train_features(self, Y, lin, groups, Xm1, splits, oof):
        """M1 features of the lambda-CV training rows a_s of each inner split s.

        Plain (G2): m_oof[a_s], whose models saw the validation rows b_s. Nested (G3): a fresh
        out-of-fold M1 pass inside a_s only (``inner_oof`` on Y.iloc[a_s]), so neither the targets
        nor the features of b_s reach them. Returns (list of (|a_s|, G) arrays, n cached)."""
        ids, Yv = Y.index, Y.to_numpy(np.float64)
        out, n_cached = [], 0
        for a, b in splits:
            if self.nested_lambda:
                sub = None if Xm1 is None else Xm1.loc[ids[a]]
                o, _, c = self.inner_oof(Y.iloc[a], lin, np.asarray(groups)[a], sub)
                out.append(self._fill(o, Yv[a])[0])
                n_cached += int(c)
            else:
                out.append(oof[a])
        return out, n_cached

    def fit(self, Y: pd.DataFrame, lineage=None, X=None, groups=None):
        t0 = time.perf_counter()
        ids = Y.index
        lin = pd.Series(lineage).reindex(ids)
        if groups is None:
            sp = pd.read_csv(ROOT / "configs" / "splits.csv", index_col=0)
            groups = sp["PatientID"].reindex(ids.astype(str)).fillna(pd.Series(ids, index=ids)).to_numpy()
        groups = np.asarray(groups)
        self.genes_ = Y.columns
        Yv = Y.to_numpy(np.float64)
        Xm1, builder = split_X(X)
        oof, splits, cached = self.inner_oof(Y, lin, groups, Xm1)
        oof, self.n_oof_filled_ = self._fill(oof, Yv)
        train_feat, n_cached = self.lambda_cv_train_features(Y, lin, groups, Xm1, splits, oof)
        t_inner = time.perf_counter() - t0
        self.m1_ = self._m1().fit(Y, lineage=lin, X=Xm1, groups=groups)
        names = PRIOR_SETS[self.priors]
        var = genesets.gene_classes(Y, set())["variable"].to_numpy()
        G = Yv.shape[1]
        local = var if self.local_genes == "variable" else np.ones(G, bool)
        gidx = np.flatnonzero(local)
        chunks = [gidx[s:s + self.chunk] for s in range(0, len(gidx), self.chunk)]
        # lambda per gene group on the inner folds (local features built per gene chunk: memory)
        r_l = np.full((len(splits), len(self.lam_grid), G), np.nan)
        for gs in chunks:
            feat = self._local(ids, builder, gs)
            for si, (a, b) in enumerate(splits):
                Xa = _stack({k: v[a] for k, v in feat.items()}, names, train_feat[si][:, gs])
                Xb = _stack({k: v[b] for k, v in feat.items()}, names, oof[np.ix_(b, gs)])
                ne = normal_equations(Xa, Yv[np.ix_(a, gs)], self.min_rows // 2)
                for li, lam in enumerate(self.lam_grid):
                    beta, icpt, ok = solve_ridge(ne, lam)
                    P = np.where(ok[None, :], apply_ridge(Xb, beta, icpt), oof[np.ix_(b, gs)])
                    r_l[si, li, gs] = _masked_r(P, Yv[np.ix_(b, gs)])
        score = np.full((len(self.lam_grid), 2), np.nan)
        self.lam_by_group_ = {}
        for j, m in enumerate([~var, var]):
            m = m & local
            if m.any():
                score[:, j] = np.nanmean(np.nanmean(r_l[:, :, m], axis=2), axis=0)
                self.lam_by_group_[j] = float(self.lam_grid[int(np.nanargmax(score[:, j]))])
        self.lam_score_ = score
        lam_gene = np.array([self.lam_by_group_.get(int(v), np.nan) for v in var])
        self.beta_ = np.zeros((G, len(names) + 1))
        self.icpt_ = np.zeros(G)
        self.ok_ = np.zeros(G, bool)
        for gs in chunks:
            Xs = _stack(self._local(ids, builder, gs), names, oof[:, gs])
            beta, icpt, ok = batched_ridge(Xs, Yv[:, gs], lam_gene[gs], self.min_rows)
            self.beta_[gs], self.icpt_[gs], self.ok_[gs] = beta, icpt, ok
        self.feature_names_ = ["m1_oof"] + names
        self.report_ = {"inner_cached": cached, "inner_s": round(t_inner, 1),
                        "nested_lambda": self.nested_lambda, "nested_cached": n_cached,
                        "local_genes": self.local_genes, "n_local_genes": int(local.sum()),
                        "lam_by_group": self.lam_by_group_, "n_passthrough": int((~self.ok_).sum()),
                        "m1": getattr(self.m1_, "report_", {}),
                        "mean_beta_m1": float(np.mean(self.beta_[self.ok_, 0])) if self.ok_.any() else None,
                        "fit_s": round(time.perf_counter() - t0, 1)}
        return self

    def predict(self, index, lineage=None, X=None) -> pd.DataFrame:
        index = pd.Index(index)
        Xm1, builder = split_X(X)
        m1 = self.m1_.predict(index, lineage=lineage, X=Xm1).to_numpy()
        names = self.feature_names_[1:]
        G = len(self.genes_)
        out = np.empty((len(index), G))
        out[:] = m1  # pass-through for genes without a local fit
        gidx = np.flatnonzero(self.ok_)
        for s in range(0, len(gidx), self.chunk):
            gs = gidx[s:s + self.chunk]
            Xs = _stack(self._local(index, builder, gs), names, m1[:, gs])
            out[:, gs] = apply_ridge(Xs, self.beta_[gs], self.icpt_[gs])
        return pd.DataFrame(out, index=index, columns=self.genes_)

    def compact(self):
        self.m1_.compact()
        return self


def _masked_r(P, T):
    """Per-column Pearson r over observed (non-NaN) T rows; < 10 rows -> NaN."""
    obs = ~np.isnan(T)
    n = obs.sum(0)
    P0 = np.where(obs, P, 0.0)
    T0 = np.where(obs, T, 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        mp, mt = P0.sum(0) / n, T0.sum(0) / n
        dp, dt = np.where(obs, P - mp, 0), np.where(obs, T - mt, 0)
        r = (dp * dt).sum(0) / np.sqrt((dp ** 2).sum(0) * (dt ** 2).sum(0))
    r = np.where(np.isfinite(r), r, 0.0)
    r[n < 10] = np.nan
    return r
