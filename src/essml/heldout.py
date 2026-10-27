"""Held-out evaluation helpers for gate G4 (cancer test set) and G6 (transfer): resampling statistics.

NOT frozen evaluation code. Point estimates always come from the frozen ``essml.evaluate`` /
``essml.metrics`` / ``essml.transfer``; this module only adds the DESIGN §4 uncertainty layer:

* ``patient_bootstrap_weights``: 95% bootstrap over held-out lines, stratified by lineage,
  resampling PATIENTS (DESIGN §4 "Uncertainty: Test"). Replicate b gives every line an integer
  weight = number of times its patient was drawn. A patient's stratum is the modal lineage of its
  lines (ties: first in sorted order); within each stratum n_patients(stratum) patients are drawn
  with replacement. Weights are shared by every model, so model-vs-baseline differences are paired.
* ``weighted_pearson_cols``: per-gene Pearson r across lines with integer line weights, equal to
  ``metrics.per_gene_pearson`` on the row-duplicated matrix (tested). Constant-prediction /
  constant-truth genes (exact, over the observed rows of the full held-out set) give r = 0, as in
  the frozen code.
* ``weighted_ap_cols``: per-gene step-wise AP (sklearn semantics, ties grouped) with line weights,
  vectorised over replicates; equals ``transfer.APPrep`` / sklearn on duplicated rows (tested).
* ``bootstrap_test_metrics``: CIs of P1 (every class), dP1 vs each baseline (paired), P2a,
  dP2a vs each baseline (paired), P2b nAP (selective, prevalence <= 0.5).

Conventions fixed here (documented in the G4 output): the evaluable gene set, class membership and
the P2b gene subset are those of the point estimate (fixed; only lines are resampled); per-line AP
(P2a) is a property of the line, so its replicate value is the weighted mean of the per-line APs.
CI = percentile (numpy 'linear'), SE = SD (ddof=1) of the replicates; every summary also reports the
bootstrap mean, its bias (mean - point) and flags a point outside the percentile interval.

P2b exception (normal interval): per-gene AP counts each item in its own precision, so a positive line
duplicated by resampling credits its own copy; resampled AP is biased upward (dev fold-0 dry run: nAP
+0.017, the percentile interval [0.200, 0.221] excluded the point 0.194, in every positive-count
stratum). That shift is an artefact of drawing duplicates, not sampling variability, so the P2b CI is
point +/- z * bootstrap SE (the percentile interval is reported alongside). P1 (Pearson) is not affected.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import metrics as M

CI_LEVEL = 0.95


# ============================================================================ resampling plan
def patient_strata(lineage: pd.Series, patient: pd.Series) -> pd.Series:
    """Stratum (lineage) per patient: modal lineage of its lines, ties -> first sorted."""
    df = pd.DataFrame({"lin": lineage.astype(str), "pat": patient.astype(str)})
    out = {}
    for p, d in df.groupby("pat", sort=True):
        vc = d["lin"].value_counts()
        out[p] = sorted(vc.index[vc == vc.max()])[0]
    return pd.Series(out, name="stratum")


def patient_bootstrap_weights(lineage: pd.Series, patient: pd.Series, n_boot: int, seed: int) -> np.ndarray:
    """Integer weights (n_boot, n_lines) for lines in ``lineage.index`` order (see module doc)."""
    ids = lineage.index
    pat = patient.reindex(ids).astype(str)
    if pat.isna().any():
        raise ValueError("patient missing for some lines")
    strata = patient_strata(lineage, pat)
    pats = strata.index.to_numpy()
    line_pat = pd.Index(pats).get_indexer(pat.to_numpy())
    rng = np.random.default_rng(np.random.SeedSequence(int(seed)))
    groups = [np.flatnonzero(strata.to_numpy() == s) for s in sorted(strata.unique())]
    W = np.empty((n_boot, len(ids)), dtype=np.float64)
    for b in range(n_boot):
        cnt = np.zeros(len(pats))
        for g in groups:  # fixed stratum order -> reproducible stream use
            cnt += np.bincount(g[rng.integers(0, len(g), len(g))], minlength=len(pats))
        W[b] = cnt[line_pat]
    return W


# ============================================================================ weighted metrics
def _const_cols(A, obs):
    big = np.finfo(np.float64).max
    return np.where(obs, A, -big).max(0) == np.where(obs, A, big).min(0)


def _pearson_prep(P, T, min_obs):
    P = np.asarray(P, np.float64)
    T = np.asarray(T, np.float64)
    obs = ~np.isnan(T)
    if np.any(~np.isfinite(P[obs])):
        raise ValueError("pred is NaN/inf at observed entries")
    n_obs = obs.sum(0)
    const = _const_cols(P, obs) | _const_cols(T, obs)
    with np.errstate(invalid="ignore", divide="ignore"):
        mp = np.where(obs, P, 0).sum(0) / n_obs
        mt = np.where(obs, T, 0).sum(0) / n_obs
    Pc = np.where(obs & ~const[None, :], P - mp, 0.0)  # centring improves precision only
    Tc = np.where(obs, T - mt, 0.0)
    return Pc, Tc, obs.astype(np.float64), const, n_obs < max(min_obs, 2)


def _pearson_chunks(P, T, W, min_obs, chunk):
    W = np.atleast_2d(np.asarray(W, np.float64))
    Pc, Tc, O, const, few = _pearson_prep(P, T, min_obs)
    PP, TT, PT = Pc * Pc, Tc * Tc, Pc * Tc
    for s in range(0, W.shape[0], chunk):
        w = W[s:s + chunk]
        sw = w @ O
        sp, st = w @ Pc, w @ Tc
        with np.errstate(invalid="ignore", divide="ignore"):
            cpp = w @ PP - sp * sp / sw
            ctt = w @ TT - st * st / sw
            cpt = w @ PT - sp * st / sw
            r = cpt / np.sqrt(cpp * ctt)
        r = np.where((cpp <= 0) | (ctt <= 0) | ~np.isfinite(r), 0.0, r)
        r = np.clip(r, -1.0, 1.0)
        r[:, const] = 0.0
        r[:, few] = np.nan
        yield s, r


def weighted_pearson_cols(P, T, W, min_obs=2, chunk=128):
    """Per-gene weighted Pearson r. P, T: (n, G) (T NaN = unobserved); W: (B, n) weights.

    Returns (B, G). Genes with < min_obs observed rows (unweighted, full set) are NaN; genes whose
    prediction or truth is exactly constant over the observed rows are 0 (frozen convention).
    """
    W = np.atleast_2d(np.asarray(W, np.float64))
    out = np.empty((W.shape[0], np.asarray(P).shape[1]))
    for s, r in _pearson_chunks(P, T, W, min_obs, chunk):
        out[s:s + len(r)] = r
    return out


def weighted_pearson_class_means(P, T, W, masks: dict, min_obs=2, chunk=128) -> dict:
    """{name: (B,) mean over the masked genes of the weighted per-gene r} without storing (B, G)."""
    W = np.atleast_2d(np.asarray(W, np.float64))
    out = {k: np.empty(W.shape[0]) for k in masks}
    for s, r in _pearson_chunks(P, T, W, min_obs, chunk):
        for k, m in masks.items():
            out[k][s:s + len(r)] = _nanmean_rows(r, np.asarray(m, bool))
    return out


def weighted_ap_cols(P, T, W, genes_mask=None, thr=-0.5):
    """Per-gene AP across lines (score -pred, label truth < thr) and prevalence, weighted.

    Returns (ap, prev), each (B, n_selected) for the genes in ``genes_mask`` (column order kept);
    AP is NaN in a replicate without a weighted positive.
    """
    P = np.asarray(P, np.float64)
    T = np.asarray(T, np.float64)
    W = np.atleast_2d(np.asarray(W, np.float64))
    B, G = W.shape[0], P.shape[1]
    gm = np.ones(G, bool) if genes_mask is None else np.asarray(genes_mask, bool)
    sel = np.flatnonzero(gm)
    ap = np.full((B, len(sel)), np.nan)
    prev = np.full((B, len(sel)), np.nan)
    for j, g in enumerate(sel):
        o = ~np.isnan(T[:, g])
        s = -P[o, g]
        y = (T[o, g] < thr).astype(np.float64)
        order = np.argsort(-s, kind="mergesort")
        ss, yy = s[order], y[order]
        ends = np.r_[np.flatnonzero(np.diff(ss) != 0), len(ss) - 1]
        w = W[:, o][:, order]
        wy, wn = w * yy[None, :], w * (1.0 - yy[None, :])
        tp = np.cumsum(wy, axis=1)[:, ends]
        tot = tp + np.cumsum(wn, axis=1)[:, ends]
        last = tp[:, -1:]
        with np.errstate(invalid="ignore", divide="ignore"):
            rec = tp / last
            prec = tp / tot
            d = np.diff(np.concatenate([np.zeros((B, 1)), rec], axis=1), axis=1)
            a = np.where(d > 0, d * prec, 0.0).sum(1)
            pv = last[:, 0] / w.sum(1)
        ap[:, j] = np.where(last[:, 0] > 0, a, np.nan)
        prev[:, j] = pv
    return ap, prev


# ============================================================================ summaries
def summarise(point, boot, level=CI_LEVEL, method="percentile") -> dict:
    """method 'percentile' (default) or 'normal' (point +/- z * bootstrap SE; used for P2b, module doc)."""
    from scipy.stats import norm
    b = np.asarray(boot, np.float64)
    ok = np.isfinite(b)
    a = (1 - level) / 2 * 100
    if ok.sum() < 2:
        return {"point": _f(point), "ci_lo": None, "ci_hi": None, "se": None, "n_boot_valid": int(ok.sum()),
                "ci_method": method}
    lo, hi = np.percentile(b[ok], [a, 100 - a])
    se = float(np.std(b[ok], ddof=1))
    mean = float(b[ok].mean())
    out = {"point": _f(point), "ci_lo": float(lo), "ci_hi": float(hi), "se": se, "n_boot_valid": int(ok.sum()),
           "ci_method": method, "boot_mean": mean,
           "bias": None if point is None or not np.isfinite(point) else mean - float(point),
           "point_outside_percentile_ci": bool(point is not None and np.isfinite(point)
                                               and not (lo - 1e-12 <= point <= hi + 1e-12))}
    if method == "normal" and point is not None and np.isfinite(point):
        z = float(norm.ppf(1 - (1 - level) / 2))
        out.update(ci_lo=float(point) - z * se, ci_hi=float(point) + z * se, ci_percentile=[float(lo), float(hi)])
    return out


def _f(x):
    return None if x is None or not np.isfinite(x) else float(x)


def _nanmean_rows(A, mask):
    with np.errstate(invalid="ignore"):
        sub = A[:, mask]
        ok = np.isfinite(sub)
        n = ok.sum(1)
        return np.where(n > 0, np.where(ok, sub, 0).sum(1) / np.maximum(n, 1), np.nan)


def bootstrap_test_metrics(pred: pd.DataFrame, truth: pd.DataFrame, gene_classes: pd.DataFrame, W: np.ndarray,
                           ev: dict, baselines: dict | None = None, min_obs=50, thr=-0.5, prev_max=0.5,
                           cache: dict | None = None) -> dict:
    """Patient-bootstrap CIs around the frozen ``evaluate`` point estimates ``ev`` (single fold).

    pred/truth/baselines: [line x gene] on the evaluated genes (truth's columns), lines in W order.
    cache: optional dict reused across models for baseline replicate arrays (memory only).
    Returns {metric_key: summarise(...)} with keys P1_<class>, dP1_vs_<B>_<class>, P2a,
    dP2a_vs_<B>, P2b_nap_prev05 (selective).
    """
    baselines = baselines or {}
    cache = {} if cache is None else cache
    genes = truth.columns
    lines = truth.index
    Tm = truth.to_numpy(np.float64)
    n_obs = (~np.isnan(Tm)).sum(0)
    evalg = n_obs >= min_obs
    cls = gene_classes.reindex(genes).fillna(False).astype(bool)
    summ = ev["summary"]

    masks = {c: cls[c].to_numpy() & evalg for c in cls.columns}

    def class_means(name, df, use_cache):
        key = ("p1", name)
        if use_cache and key in cache:
            return cache[key]
        v = weighted_pearson_class_means(df.loc[lines, genes].to_numpy(np.float64), Tm, W, masks, min_obs)
        if use_cache:
            cache[key] = v
        return v

    out = {}
    R = class_means("model", pred, False)
    RB = {b: class_means(b, d, True) for b, d in baselines.items()}
    # dP1 = mean over genes of (r - r_B) = difference of the class means: both are NaN exactly where the
    # gene has < min_obs observed lines (truth-only property), so the gene sets coincide.
    for c in cls.columns:
        out[f"P1_{c}"] = summarise(summ[c]["P1"], R[c])
        for b, rb in RB.items():
            out[f"dP1_vs_{b}_{c}"] = summarise(summ[c].get(f"dP1_vs_{b}", np.nan), R[c] - rb[c])

    # P2a: per-line AP over the 'all' class genes (fixed per line) -> weighted mean
    pl = ev["per_line"].set_index("ModelID").reindex(lines)
    ap = pl["ap"].to_numpy(np.float64)

    def wmean(v):
        ok = np.isfinite(v)
        return (W[:, ok] @ v[ok]) / W[:, ok].sum(1)

    out["P2a"] = summarise(summ["all"]["P2a"], wmean(ap))
    for b in baselines:
        if f"ap_{b}" in pl:
            out[f"dP2a_vs_{b}"] = summarise(summ["all"].get(f"dP2a_vs_{b}", np.nan),
                                            wmean(ap - pl[f"ap_{b}"].to_numpy(np.float64)))

    # P2b headline: selective genes, >= 3 pos / >= 3 neg, prevalence <= prev_max (fixed at the point estimate)
    if "selective" in cls:
        pg = ev["per_gene"].set_index("gene").reindex(genes)
        e05 = (cls["selective"].to_numpy() & evalg & np.isfinite(pg["nap"].to_numpy(np.float64))
               & (pg["prev"].to_numpy(np.float64) <= prev_max))
        apb, pv = weighted_ap_cols(pred.loc[lines, genes].to_numpy(np.float64), Tm, W, e05, thr)
        with np.errstate(invalid="ignore", divide="ignore"):
            nap = (apb - pv) / (1.0 - pv)
        out["P2b_nap_prev05"] = summarise(summ["selective"]["P2b_nap_prev05"],
                                          _nanmean_rows(nap, np.ones(nap.shape[1], bool)), method="normal")
        out["_P2b_n_genes"] = int(e05.sum())
    return out
