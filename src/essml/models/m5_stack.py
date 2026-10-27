"""M5: stacking with non-negative weights per gene class (DESIGN.md §5; G3 close-out).

Per outer training fold (training lines only reach any statistic):
1. Inner folds: ``inner_folds`` patient-grouped, lineage-stratified folds of the training lines
   (the splitter and seed of ``LocalStageM2.inner_oof``, so the M1 inner OOF already cached in
   cache/m2_inner by the M2/M3 runs is reused through its exact key).
2. Inner out-of-fold predictions for every input model: for inner split s the input is fitted on
   the inner-training lines a_s and predicts the inner-held-out lines b_s. Every training line
   therefore gets a prediction from a model that never saw it. M1 inputs reuse
   ``LocalStageM2.inner_oof`` (exact-key cache cache/m2_inner); other inputs are fitted here and
   cached under an exact key in cache/m5_inner (config, row order, columns, target values,
   lineage and patient labels in row order, version). ``inner_params`` may override an input's
   params for these inner fits only (a declared simplification; see the config).
3. Weights: training-fold gene classes (``genesets.gene_classes`` on the training lines) are
   partitioned into ``CLASS_ORDER``. Per class, non-negative least squares of the training
   targets on the inner OOF predictions, per gene centred on its observed training lines and
   scaled by 1/SD of its training targets (each gene weighs equally, as in P1; the centring makes
   the fit about variation across lines, which is what P1 scores). Weights are then normalised to
   sum to 1 (r is scale-invariant; this keeps the prediction level). A class with fewer than
   ``min_class_genes`` genes uses the weights pooled over all variable genes; a zero solution
   falls back to equal weights.
4. The inputs are refitted on all training lines (lazily, at predict; exact-key cache of their
   predictions in cache/m5_outer) and combined with the class weights.

``mode="average"``: the plain-average baseline (equal weights; steps 1-3 skipped).
``permute_input``: sanity control. That input's predictions are row-permuted across lines (inner
OOF among the training lines, outer predictions among the predicted lines), so it carries no
line information and should get ~0 weight.
``oof_log_`` records (input, split, fit IDs, predicted IDs) for every inner OOF block (tests).
"""
from __future__ import annotations

import hashlib
import json
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import nnls

from essml import genesets

ROOT = Path(__file__).resolve().parents[3]
INNER_CACHE = ROOT / "cache" / "m5_inner"
OUTER_CACHE = ROOT / "cache" / "m5_outer"
CACHE_VERSION = 1
CLASS_ORDER = ["nonvariable", "var_common_essential", "var_selective", "var_other"]


# ============================================================================ weights
def class_partition(Y: pd.DataFrame, neg=None) -> np.ndarray:
    """Training-fold class label per gene (CLASS_ORDER): non-variable genes; variable genes split
    by priority common essential > selective (5-95%) > other variable."""
    gc = genesets.gene_classes(Y, set() if neg is None else neg)
    lab = np.full(Y.shape[1], "nonvariable", dtype=object)
    var = gc["variable"].to_numpy()
    lab[var] = "var_other"
    lab[var & gc["selective"].to_numpy()] = "var_selective"
    lab[var & gc["common_essential"].to_numpy()] = "var_common_essential"
    return lab


def gram(P_list, Y):
    """Per-gene centred, 1/SD(y)-scaled cross-products over observed rows.

    P_list: M arrays (n, G) (no NaN); Y: (n, G) with NaN. Returns (A (G, M, M), b (G, M))."""
    Y = np.asarray(Y, np.float64)
    obs = ~np.isnan(Y)
    n = obs.sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        my = np.where(obs, Y, 0).sum(0) / n
        yc = np.where(obs, Y - my, 0.0)
        sd = np.sqrt((yc ** 2).sum(0) / np.maximum(n - 1, 1))
    scale = np.where((n > 1) & (sd > 1e-8), 1.0 / sd, 0.0)
    yc = yc * scale
    Pc = []
    for P in P_list:
        P = np.asarray(P, np.float64)
        with np.errstate(invalid="ignore", divide="ignore"):
            mp = np.where(obs, P, 0.0).sum(0) / n  # unobserved rows never enter
        Pc.append(np.where(obs, P - np.nan_to_num(mp), 0.0) * scale)
    M = len(Pc)
    A = np.empty((Y.shape[1], M, M))
    b = np.empty((Y.shape[1], M))
    for i in range(M):
        b[:, i] = np.einsum("ij,ij->j", Pc[i], yc)
        for j in range(i, M):
            A[:, i, j] = A[:, j, i] = np.einsum("ij,ij->j", Pc[i], Pc[j])
    return A, b


def nnls_gram(A, b):
    """argmin_w>=0 w'Aw - 2b'w for a PSD (M, M) A (tiny ridge for rank-deficient A, e.g. identical
    inputs). Returns w (M,)."""
    M = len(b)
    eps = 1e-10 * max(float(np.trace(A)) / M, 1e-300)
    L = np.linalg.cholesky(A + eps * np.eye(M))
    w, _ = nnls(L.T, np.linalg.solve(L, b))
    return w


def normalise(w):
    s = float(np.sum(w))
    return w / s if s > 0 else np.full(len(w), 1.0 / len(w))


def class_weights(P_list, Y, labels, min_class_genes=20):
    """Non-negative, sum-to-1 weights per class label (dict label -> (M,) array) + report."""
    A, b = gram(P_list, Y)
    labels = np.asarray(labels)
    var = labels != "nonvariable"
    pooled = normalise(nnls_gram(A[var].sum(0), b[var].sum(0))) if var.any() else None
    out, rep = {}, {}
    for lab in CLASS_ORDER:
        m = labels == lab
        k = int(m.sum())
        if k == 0:
            continue
        if k < min_class_genes and pooled is not None:
            out[lab], how = pooled, "pooled_variable"
        else:
            raw = nnls_gram(A[m].sum(0), b[m].sum(0))
            out[lab], how = normalise(raw), ("nnls" if raw.sum() > 0 else "equal_fallback")
        rep[lab] = {"n_genes": k, "weights": np.round(out[lab], 6).tolist(), "how": how}
    return out, rep


# ============================================================================ input models
def input_class(model_key):
    from .kernel_ridge import KernelRidgeM1
    from .lgbm import LGBMResidualM3
    from .m2_local import LocalStageM2
    return {"M1": KernelRidgeM1, "M2": LocalStageM2, "M3": LGBMResidualM3}[model_key]


def _hash_data(h, Y, lineage, groups):
    h.update(Y.columns.astype(str).str.cat(sep=",").encode())
    h.update(np.ascontiguousarray(Y.to_numpy(np.float64)).tobytes())
    h.update(pd.Series(lineage).reindex(Y.index).astype(str).str.cat(sep=",").encode())
    h.update(",".join(map(str, np.asarray(groups))).encode())
    return h


def _fill_mean(P, Yv):
    P = np.asarray(P, np.float64)
    bad = ~np.isfinite(P)
    if bad.any():
        mu = np.nan_to_num(np.nanmean(Yv, axis=0))
        P = np.where(bad, mu[None, :], P)
    return P


class StackM5:
    name = "M5"
    uses_groups = True

    def __init__(self, inputs, mode="stack", inner_folds=5, seed=20261026, min_class_genes=20,
                 permute_input=None, permute_seed=0, cache=True, input_factory=None):
        """inputs: list of {"name", "model" (M1/M2/M3), "params", optional "inner_params"}.
        input_factory (tests): callable(spec, inner: bool) -> unfitted model."""
        if mode not in ("stack", "average"):
            raise ValueError("mode must be 'stack' or 'average'")
        self.inputs = [dict(s) for s in inputs]
        names = [s["name"] for s in self.inputs]
        if len(set(names)) != len(names):
            raise ValueError("input names must be unique")
        if permute_input is not None and permute_input not in names:
            raise ValueError(f"permute_input {permute_input} is not an input")
        self.mode = mode
        self.inner_folds, self.seed = int(inner_folds), int(seed)
        self.min_class_genes = int(min_class_genes)
        self.permute_input, self.permute_seed = permute_input, int(permute_seed)
        self.cache = bool(cache) and input_factory is None
        self.input_factory = input_factory

    # -------------------------------------------------------------- helpers
    def _make(self, spec, inner):
        if self.input_factory is not None:
            return self.input_factory(spec, inner)
        params = dict(spec.get("params", {}))
        if inner:
            params.update(spec.get("inner_params", {}))
        return input_class(spec["model"])(**params)

    def _splitter(self):
        from .m2_local import LocalStageM2
        return LocalStageM2(m1_params={"modalities": ["E"], "residual": "B0"},
                            inner_folds=self.inner_folds, seed=self.seed)

    def _is_cached_m1(self, spec):
        """M1 inputs use LocalStageM2.inner_oof (exact-key cache shared with M2/M3)."""
        return (self.input_factory is None and spec["model"] == "M1"
                and not spec.get("inner_params"))

    def _inner_key(self, spec, Y, lineage, groups):
        h = hashlib.sha256(json.dumps({"ids": list(map(str, Y.index)), "model": spec["model"],
                                       "params": spec.get("params", {}),
                                       "inner_params": spec.get("inner_params", {}),
                                       "k": self.inner_folds, "seed": self.seed,
                                       "v": CACHE_VERSION}, sort_keys=True, default=str).encode())
        return f"{spec['name']}_{_hash_data(h, Y, lineage, groups).hexdigest()[:24]}"

    def _outer_key(self, spec, index):
        h = hashlib.sha256(json.dumps({"ids": list(map(str, self.Y_.index)), "model": spec["model"],
                                       "params": spec.get("params", {}), "pred": list(map(str, index)),
                                       "v": CACHE_VERSION}, sort_keys=True, default=str).encode())
        return f"{spec['name']}_{_hash_data(h, self.Y_, self.lin_, self.groups_).hexdigest()[:24]}"

    def _inner_oof(self, spec, Y, lin, groups, splits):
        ids = Y.index
        if self._is_cached_m1(spec):
            st = self._splitter()
            st.m1_params = dict(spec.get("params", {}))
            oof, sp2, cached = st.inner_oof(Y, lin, groups)
            if len(sp2) != len(splits) or not all(
                    np.array_equal(a, c) and np.array_equal(b, d) for (a, b), (c, d) in zip(splits, sp2)):
                raise AssertionError("M1 inner-OOF splits differ from the stacking splits")
            for s, (a, b) in enumerate(splits):
                self.oof_log_.append({"input": spec["name"], "split": s, "fit": set(ids[a]),
                                      "pred": set(ids[b])})
            return np.asarray(oof, np.float64), bool(cached)
        path = None
        if self.cache:
            path = INNER_CACHE / f"{self._inner_key(spec, Y, lin, groups)}.npy"
            if path.exists():
                for s, (a, b) in enumerate(splits):
                    self.oof_log_.append({"input": spec["name"], "split": s, "fit": set(ids[a]),
                                          "pred": set(ids[b])})
                return np.load(path).astype(np.float64), True
        oof = np.full(Y.shape, np.nan)
        grp = np.asarray(groups)
        for s, (a, b) in enumerate(splits):
            t0 = time.perf_counter()
            m = self._make(spec, inner=True)
            m.fit(Y.iloc[a], lineage=lin.iloc[a], X=None, groups=grp[a])
            oof[b] = m.predict(ids[b], lineage=lin.iloc[b]).to_numpy(np.float64)
            self.oof_log_.append({"input": spec["name"], "split": s, "fit": set(ids[a]),
                                  "pred": set(ids[b])})
            self.timing_.setdefault(spec["name"], []).append(round(time.perf_counter() - t0, 1))
            rep = getattr(m, "report_", None)
            if rep is not None:
                self.inner_reports_.setdefault(spec["name"], []).append(
                    {k: rep.get(k) for k in ("lam_by_group", "frac_genes_no_correction",
                                             "median_best_iter", "nested_es", "nested_lambda",
                                             "inner_cached", "fit_s") if k in rep})
            del m
        if path is not None:
            INNER_CACHE.mkdir(parents=True, exist_ok=True)
            tmp = path.with_name(f"{path.stem}.{os.getpid()}.tmp.npy")
            np.save(tmp, oof.astype(np.float32))
            os.replace(tmp, path)
            oof = np.load(path).astype(np.float64)  # identical values whether cached or fresh
        return oof, False

    def _perm(self, n, salt):
        return np.random.default_rng([self.permute_seed, salt]).permutation(n)

    # -------------------------------------------------------------- fit / predict
    def fit(self, Y: pd.DataFrame, lineage=None, X=None, groups=None):
        t0 = time.perf_counter()
        if X is not None and self.input_factory is None:
            raise ValueError("M5 inputs read their features from cache/features (X must be None)")
        ids = Y.index
        lin = (pd.Series("x", index=ids) if lineage is None
               else pd.Series(lineage).reindex(ids).astype(str))
        if groups is None:
            from .g3_common import default_groups
            groups = default_groups(ids)
        groups = np.asarray(groups)
        self.Y_, self.lin_, self.groups_, self.X_ = Y, lin, groups, X
        self.genes_ = Y.columns
        self.names_ = [s["name"] for s in self.inputs]
        self.labels_ = class_partition(Y)
        self.oof_log_, self.timing_, self.inner_reports_ = [], {}, {}
        M = len(self.inputs)
        rep = {"mode": self.mode, "class_n_genes": {c: int((self.labels_ == c).sum())
                                                    for c in CLASS_ORDER}}
        if self.mode == "average":
            self.weights_ = {c: np.full(M, 1.0 / M) for c in CLASS_ORDER}
            self.report_ = dict(rep, fit_s=round(time.perf_counter() - t0, 1))
            return self
        splits = self._splitter()._inner_splits(ids, lin, groups)
        splits = [(np.asarray(a), np.asarray(b)) for a, b in splits]
        Yv = Y.to_numpy(np.float64)
        oofs, cached = [], {}
        for spec in self.inputs:
            o, c = self._inner_oof(spec, Y, lin, groups, splits)
            o = _fill_mean(o, Yv)
            if spec["name"] == self.permute_input:
                o = o[self._perm(len(ids), 0)]
            oofs.append(o)
            cached[spec["name"]] = c
        self.check_nesting()
        if self.input_factory is not None:  # tests only (real runs: ~120 MB per input and fold)
            self.oof_ = dict(zip(self.names_, oofs))
        self.weights_, wrep = class_weights(oofs, Yv, self.labels_, self.min_class_genes)
        from .g3_common import masked_r
        var = self.labels_ != "nonvariable"
        oof_r = {n: float(np.nanmean(np.nanmean(
            [masked_r(o[b][:, var], Yv[b][:, var], 20) for _, b in splits], axis=0)))
            for n, o in zip(self.names_, oofs)}
        self.report_ = dict(rep, weights=wrep, inner_cached=cached, inner_s=self.timing_,
                            inner_reports=self.inner_reports_, inner_oof_P1_variable=oof_r,
                            permute_input=self.permute_input,
                            fit_s=round(time.perf_counter() - t0, 1))
        return self

    def check_nesting(self):
        """Every training line's inner OOF prediction came from a fit that did not contain it, and
        every line is predicted exactly once per input."""
        for name in self.names_:
            seen = set()
            for e in (e for e in self.oof_log_ if e["input"] == name):
                if e["fit"] & e["pred"]:
                    raise AssertionError(f"{name}: inner fit contains predicted lines")
                if e["pred"] & seen:
                    raise AssertionError(f"{name}: a line was predicted twice")
                seen |= e["pred"]
            if seen != set(self.Y_.index):
                raise AssertionError(f"{name}: inner OOF does not cover the training lines")

    def input_predictions(self, index, lineage=None):
        """Outer predictions of every input (refitted on all training lines; exact-key cache)."""
        index = pd.Index(index)
        out = {}
        lin_new = None if lineage is None else pd.Series(lineage).reindex(index)
        for spec in self.inputs:
            path = OUTER_CACHE / f"{self._outer_key(spec, index)}.npy" if self.cache else None
            if path is not None and path.exists():
                P = np.load(path)
            else:
                t0 = time.perf_counter()
                m = self._make(spec, inner=False)
                m.fit(self.Y_, lineage=self.lin_, X=self.X_, groups=self.groups_)
                P = m.predict(index, lineage=lin_new).to_numpy(np.float64)
                self.timing_.setdefault("outer_" + spec["name"], []).append(
                    round(time.perf_counter() - t0, 1))
                del m
                if path is not None:
                    OUTER_CACHE.mkdir(parents=True, exist_ok=True)
                    tmp = path.with_name(f"{path.stem}.{os.getpid()}.tmp.npy")
                    np.save(tmp, P)
                    os.replace(tmp, path)
            out[spec["name"]] = P
        return out

    def predict(self, index, lineage=None, X=None) -> pd.DataFrame:
        index = pd.Index(index)
        preds = self.input_predictions(index, lineage)
        out = np.zeros((len(index), len(self.genes_)))
        for k, name in enumerate(self.names_):
            P = preds[name]
            if name == self.permute_input:
                P = P[self._perm(len(index), 1)]
            for c, w in self.weights_.items():
                m = self.labels_ == c
                if w[k] != 0:
                    out[:, m] += w[k] * P[:, m]
        return pd.DataFrame(out, index=index, columns=self.genes_)
