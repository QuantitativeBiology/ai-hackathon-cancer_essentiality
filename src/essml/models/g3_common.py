"""Shared helpers for the G3 nonlinear track (M3, M4, R2, R4; DESIGN.md §5).

Everything here is fitted on TRAINING lines only:

- ``CellFeatures``: raw blocks from cache/features -> ``features.FoldTransformer`` fitted on the
  training IDs (standardisation, SD/min-count filters) -> clipped dense matrix. Frame mode (X given
  as a DataFrame, tests) treats X as one numeric block.
- ``TrainPCA``: PCA (economy SVD) of the training rows of that matrix; held-out rows are only
  projected. Scores are divided by the training SD of PC1 (one global scale, so the relative
  variance of the PCs is kept) or, with ``whiten=True``, by each PC's own training SD.
- ``inner_split``: ONE grouped (PatientID) and lineage-stratified 80/20 split of the training lines,
  for early stopping / hyper-parameter choice.
- ``masked_r``: per-gene Pearson r over observed rows (NaN-aware), used for inner criteria only
  (the reported metrics come from the frozen ``essml.evaluate``).
- ``make_pool``: process pool for per-gene fits (1 thread per process).
"""
from __future__ import annotations

import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from essml import features as F

ROOT = Path(__file__).resolve().parents[3]
_RAW: dict[str, pd.DataFrame] = {}


def load_blocks(names, feature_dir=F.FEATURE_DIR) -> dict[str, pd.DataFrame]:
    """Raw row-wise blocks (memoised per process)."""
    out = {}
    for n in names:
        if n not in _RAW:
            _RAW[n] = F.load_raw([n], feature_dir)[n]
        out[n] = _RAW[n]
    return out


def default_groups(ids):
    """PatientID per training line from configs/splits.csv (singletons where unknown)."""
    sp = pd.read_csv(ROOT / "configs" / "splits.csv", index_col=0)
    ids = pd.Index([str(i) for i in ids])
    return sp["PatientID"].reindex(ids).fillna(pd.Series(ids, index=ids)).astype(str).to_numpy()


class CellFeatures:
    """Line-level design matrix fitted on training IDs (``FoldTransformer``), clipped at +-clip.

    blocks: raw block names (cache/features) with DESIGN §2 default specs. ``permute_seed``
    (negative control): the same random permutation of ModelIDs is applied to every raw block, so
    each line receives another line's features (train and held-out alike)."""

    def __init__(self, blocks=("E",), clip=8.0, permute_seed=None):
        self.blocks = list(blocks)
        self.clip = float(clip)
        self.permute_seed = permute_seed

    def _raw(self, X):
        if isinstance(X, pd.DataFrame):
            raw = {"X": X}
        else:
            raw = dict(load_blocks(self.blocks))
        if self.permute_seed is not None:
            allids = sorted(set().union(*[set(b.index) for b in raw.values()]))
            perm = np.random.default_rng(int(self.permute_seed)).permutation(len(allids))
            remap = dict(zip(allids, np.asarray(allids)[perm]))
            raw = {k: v.rename(index=remap) for k, v in raw.items()}
        return raw

    def fit(self, train_ids, X=None):
        self.frame_mode_ = isinstance(X, pd.DataFrame)
        if self.frame_mode_:
            specs = {"X": F.BlockSpec()}
        else:
            allspec = F.default_specs(None)
            specs = {b: allspec[b] for b in self.blocks}
        self.ft_ = F.FoldTransformer(specs, check_roles=not self.frame_mode_)
        self.ft_.fit(self._raw(X), train_ids)
        Z = self.transform(train_ids, X, _first=True)
        return self, Z

    def transform(self, ids, X=None, _first=False):
        df = self.ft_.transform(self._raw(X), ids)
        if _first:
            self.columns_ = list(df.columns)
        elif list(df.columns) != self.columns_:
            raise ValueError("feature columns changed between fit and transform")
        Z = df.to_numpy(np.float32)
        np.clip(Z, -self.clip, self.clip, out=Z)
        return Z


class TrainPCA:
    """PCA on training rows only. transform(Z) = (Z - mean_train) @ V_k / scale."""

    def __init__(self, n_comp=256, whiten=False):
        self.n_comp = int(n_comp)
        self.whiten = bool(whiten)

    def fit(self, Z):
        Z = np.asarray(Z, np.float64)
        self.mean_ = Z.mean(0)
        Zc = Z - self.mean_
        # economy SVD via the n x n Gram matrix (n << p)
        K = Zc @ Zc.T
        d, U = np.linalg.eigh(K)
        order = np.argsort(d)[::-1]
        d, U = np.clip(d[order], 0, None), U[:, order]
        k = min(self.n_comp, int((d > 1e-8 * d[0]).sum()))
        s = np.sqrt(d[:k])
        self.components_ = (Zc.T @ U[:, :k]) / s  # (p, k), orthonormal columns
        sd = s / np.sqrt(len(Z))                    # training SD of each PC score
        self.pc_sd_ = sd
        self.scale_ = sd if self.whiten else np.full(k, sd[0])
        self.explained_ = d[:k] / d.sum()
        return self

    def transform(self, Z):
        return ((np.asarray(Z, np.float64) - self.mean_) @ self.components_ / self.scale_).astype(np.float32)


def inner_split(ids, groups, lineage=None, seed=20261026, n_splits=5):
    """One grouped (patient) 1/n_splits validation split of the training lines.

    Returns (train_pos, val_pos) integer positions into ``ids``."""
    from sklearn.model_selection import StratifiedGroupKFold
    ids = pd.Index(ids)
    strata = (np.zeros(len(ids), int) if lineage is None
              else pd.Series(lineage).reindex(ids).astype(str).to_numpy())
    grp = np.asarray(groups if groups is not None else np.arange(len(ids))).astype(str)
    sgk = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=int(seed))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        a, b = next(iter(sgk.split(np.zeros(len(ids)), strata, grp)))
    return np.sort(a), np.sort(b)


def grouped_kfold(ids, groups, lineage=None, seed=20261026, n_splits=3):
    """All splits of a grouped, lineage-stratified K-fold (positions)."""
    from sklearn.model_selection import StratifiedGroupKFold
    ids = pd.Index(ids)
    strata = (np.zeros(len(ids), int) if lineage is None
              else pd.Series(lineage).reindex(ids).astype(str).to_numpy())
    grp = np.asarray(groups if groups is not None else np.arange(len(ids))).astype(str)
    sgk = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=int(seed))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        return [(np.sort(a), np.sort(b)) for a, b in sgk.split(np.zeros(len(ids)), strata, grp)]


def masked_stats(Y):
    """Per-gene mean and SD (ddof=1) over observed rows; SD < 1e-8 or n < 2 -> 1."""
    Y = np.asarray(Y, np.float64)
    obs = ~np.isnan(Y)
    n = obs.sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        mu = np.where(obs, Y, 0).sum(0) / n
        sd = np.sqrt(np.where(obs, (Y - mu) ** 2, 0).sum(0) / (n - 1))
    mu = np.nan_to_num(mu)
    sd = np.where((n > 1) & (sd > 1e-8), sd, 1.0)
    return mu, sd


def masked_r(P, T, min_obs=10):
    """Per-column Pearson r over rows where T is observed. Constant -> 0; < min_obs rows -> NaN."""
    P = np.asarray(P, np.float64)
    T = np.asarray(T, np.float64)
    obs = ~np.isnan(T)
    n = obs.sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        mp = np.where(obs, P, 0).sum(0) / n
        mt = np.where(obs, T, 0).sum(0) / n
        dp = np.where(obs, P - mp, 0)
        dt = np.where(obs, T - mt, 0)
        den = np.sqrt((dp ** 2).sum(0) * (dt ** 2).sum(0))
        r = (dp * dt).sum(0) / den
    r = np.where(den > 0, r, 0.0)
    return np.where(n >= min_obs, r, np.nan)


def variable_mask(Y, sd_thr=0.25):
    """Training-fold variable genes (SD > 0.25 over observed rows; = genesets.gene_classes)."""
    Y = np.asarray(Y, np.float64)
    obs = ~np.isnan(Y)
    n = obs.sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        mu = np.where(obs, Y, 0).sum(0) / n
        sd = np.sqrt(np.where(obs, (Y - mu) ** 2, 0).sum(0) / (n - 1))
    return (n > 1) & (sd > sd_thr)


def top_abs_corr(X, y, k):
    """Indices of the k columns of X with the largest |Pearson r| with y (rows already observed)."""
    Xc = X - X.mean(0)
    yc = y - y.mean()
    num = yc @ Xc
    den = np.sqrt((Xc ** 2).sum(0)) * np.sqrt((yc ** 2).sum())
    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.where(den > 0, num / den, 0.0)
    a = np.abs(r)
    k = min(k, X.shape[1])
    idx = np.argpartition(-a, k - 1)[:k]
    return idx[np.argsort(-a[idx], kind="stable")]


def set_single_thread():
    for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[v] = "1"
    try:
        from threadpoolctl import threadpool_limits
        threadpool_limits(1)
    except ImportError:
        pass


class _SerialPool:
    def __init__(self, initializer=None, initargs=()):
        if initializer is not None:
            initializer(*initargs)

    def map(self, func, tasks, chunksize=1):
        return map(func, tasks)

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


def make_pool(n_jobs, initializer=None, initargs=()):
    """Process pool of ``n_jobs`` single-thread workers (context manager with ``.map``).

    n_jobs <= 1 runs in-process (tests). Large read-only arrays go through ``initializer`` once
    per worker instead of with every task."""
    if n_jobs <= 1:
        return _SerialPool(initializer, initargs)
    from concurrent.futures import ProcessPoolExecutor
    return ProcessPoolExecutor(max_workers=n_jobs, initializer=_pool_init,
                               initargs=(initializer, initargs))


def _pool_init(initializer, initargs):
    set_single_thread()
    if initializer is not None:
        initializer(*initargs)
