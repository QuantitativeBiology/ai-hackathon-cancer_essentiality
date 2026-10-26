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
Genes with < ``min_rows`` observed training rows pass the M1 prediction through.

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
                 lam_grid=None, min_rows=50, chunk=1500, seed=20261026, cache_inner=True):
        self.m1_params = dict(m1_params or {"modalities": ["E"]})
        self.base = base
        if priors not in PRIOR_SETS:
            raise ValueError(f"priors must be one of {list(PRIOR_SETS)}")
        self.priors = priors
        self.inner_folds = int(inner_folds)
        self.lam_grid = np.asarray(LOCAL_LAM if lam_grid is None else lam_grid, np.float64)
        self.min_rows, self.chunk, self.seed = int(min_rows), int(chunk), int(seed)
        self.cache_inner = bool(cache_inner)

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

    def fit(self, Y: pd.DataFrame, lineage=None, X=None, groups=None):
        t0 = time.perf_counter()
        ids = Y.index
        lin = pd.Series(lineage).reindex(ids)
        if groups is None:
            sp = pd.read_csv(ROOT / "configs" / "splits.csv", index_col=0)
            groups = sp["PatientID"].reindex(ids.astype(str)).fillna(pd.Series(ids, index=ids)).to_numpy()
        self.genes_ = Y.columns
        Xm1, builder = split_X(X)
        oof, splits, cached = self.inner_oof(Y, lin, groups, Xm1)
        oof = np.asarray(oof, np.float64)
        bad = ~np.isfinite(oof)
        if bad.any():  # genes unobserved in an inner training part: use the training gene mean
            mu = np.nan_to_num(np.nanmean(Y.to_numpy(np.float64), axis=0))
            oof = np.where(bad, mu[None, :], oof)
        self.n_oof_filled_ = int(bad.sum())
        t_inner = time.perf_counter() - t0
        self.m1_ = self._m1().fit(Y, lineage=lin, X=Xm1, groups=groups)
        names = PRIOR_SETS[self.priors]
        Yv = Y.to_numpy(np.float64)
        var = genesets.gene_classes(Y, set())["variable"].to_numpy()
        G = Yv.shape[1]
        chunks = [slice(s, min(s + self.chunk, G)) for s in range(0, G, self.chunk)]
        # lambda per gene group on the inner folds (local features built per gene chunk: memory)
        r_l = np.full((len(splits), len(self.lam_grid), G), np.nan)
        for gs in chunks:
            Xs = _stack(self._local(ids, builder, gs), names, oof[:, gs])
            for si, (a, b) in enumerate(splits):
                ne = normal_equations(Xs[:, a], Yv[a, gs], self.min_rows // 2)
                Xb = Xs[:, b]
                for li, lam in enumerate(self.lam_grid):
                    beta, icpt, ok = solve_ridge(ne, lam)
                    P = np.where(ok[None, :], apply_ridge(Xb, beta, icpt), oof[b, gs])
                    r_l[si, li, gs] = _masked_r(P, Yv[b, gs])
        score = np.zeros((len(self.lam_grid), 2))
        for j, m in enumerate([~var, var]):
            score[:, j] = np.nanmean(np.nanmean(r_l[:, :, m], axis=2), axis=0)
        self.lam_by_group_ = {0: float(self.lam_grid[int(np.argmax(score[:, 0]))]),
                              1: float(self.lam_grid[int(np.argmax(score[:, 1]))])}
        self.lam_score_ = score
        lam_gene = np.where(var, self.lam_by_group_[1], self.lam_by_group_[0])
        self.beta_ = np.zeros((G, len(names) + 1))
        self.icpt_ = np.zeros(G)
        self.ok_ = np.zeros(G, bool)
        for gs in chunks:
            Xs = _stack(self._local(ids, builder, gs), names, oof[:, gs])
            beta, icpt, ok = batched_ridge(Xs, Yv[:, gs], lam_gene[gs], self.min_rows)
            self.beta_[gs], self.icpt_[gs], self.ok_[gs] = beta, icpt, ok
        self.feature_names_ = ["m1_oof"] + names
        self.report_ = {"inner_cached": cached, "inner_s": round(t_inner, 1),
                        "lam_by_group": self.lam_by_group_, "n_passthrough": int((~self.ok_).sum()),
                        "m1": getattr(self.m1_, "report_", {}),
                        "mean_beta_m1": float(np.mean(self.beta_[self.ok_, 0])),
                        "fit_s": round(time.perf_counter() - t0, 1)}
        return self

    def predict(self, index, lineage=None, X=None) -> pd.DataFrame:
        index = pd.Index(index)
        Xm1, builder = split_X(X)
        m1 = self.m1_.predict(index, lineage=lineage, X=Xm1).to_numpy()
        names = self.feature_names_[1:]
        G = len(self.genes_)
        out = np.empty((len(index), G))
        for s in range(0, G, self.chunk):
            gs = slice(s, min(s + self.chunk, G))
            Xs = _stack(self._local(index, builder, gs), names, m1[:, gs])
            P = apply_ridge(Xs, self.beta_[gs], self.icpt_[gs])
            out[:, gs] = np.where(self.ok_[gs][None, :], P, m1[:, gs])
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
