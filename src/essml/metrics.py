"""Evaluation metrics (FROZEN: hashed in configs/frozen_hashes.json; see DESIGN.md §4, §9).

Conventions
-----------
All functions take numpy arrays shaped (n_lines, n_genes) unless stated otherwise.
``true`` may contain NaN (library-structured CRISPR NaNs); a NaN in ``true`` means "not
observed" and the entry is excluded. NaNs are never filled. ``pred`` must be finite wherever
``true`` is observed; a NaN prediction on an observed entry raises ``ValueError`` (a model that
silently returns NaN must not be scored on fewer entries). Computation is in float64.

Degenerate cases (DESIGN §4):
- Pearson/Spearman: a constant prediction or a constant truth over the observed entries gives
  r = 0 (not NaN). "Constant" is exact: max - min == 0 over the observed entries.
- A gene (or line) with fewer than ``min_obs`` observed entries gives NaN = excluded.
- AP is sklearn's step-wise ``average_precision_score``; a constant score gives AP = π exactly.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import rankdata
from sklearn.metrics import average_precision_score, roc_auc_score


def _prep(pred, true):
    """float64 copies + observed mask; raise on NaN predictions at observed entries."""
    pred = np.asarray(pred, dtype=np.float64)
    true = np.asarray(true, dtype=np.float64)
    if pred.shape != true.shape:
        raise ValueError(f"shape mismatch: pred {pred.shape} vs true {true.shape}")
    obs = ~np.isnan(true)
    if np.any(~np.isfinite(pred[obs])):
        raise ValueError("pred is NaN/inf at observed entries")
    return pred, true, obs


def _masked_pearson_cols(P, T, obs, min_obs):
    """Column-wise Pearson r over observed rows. Constant column -> 0; n_obs < min_obs -> NaN."""
    n = obs.sum(0)
    P0 = np.where(obs, P, 0.0)
    T0 = np.where(obs, T, 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        mp = P0.sum(0) / n
        mt = T0.sum(0) / n
    dp = np.where(obs, P - mp, 0.0)
    dt = np.where(obs, T - mt, 0.0)
    num = (dp * dt).sum(0)
    den = np.sqrt((dp ** 2).sum(0) * (dt ** 2).sum(0))
    # exact constancy test on observed entries (robust to rounding in the mean)
    big = np.finfo(np.float64).max
    const_p = np.where(obs, P, -big).max(0) == np.where(obs, P, big).min(0)
    const_t = np.where(obs, T, -big).max(0) == np.where(obs, T, big).min(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        r = num / den
    r = np.where(const_p | const_t | (den == 0), 0.0, r)
    r = np.clip(r, -1.0, 1.0)
    r[n < max(min_obs, 2)] = np.nan
    return r


def per_gene_pearson(pred, true, min_obs=40):
    """Pearson r per gene (column) across observed lines -> array (n_genes,)."""
    P, T, obs = _prep(pred, true)
    return _masked_pearson_cols(P, T, obs, min_obs)


def _masked_ranks(A, obs):
    """Average ranks per column among observed rows (unobserved rows get rank NaN)."""
    A = np.where(obs, A, np.inf)  # unobserved sort last, so observed ranks are 1..n_obs
    R = rankdata(A, axis=0, method="average")
    return np.where(obs, R, np.nan)


def per_gene_spearman(pred, true, min_obs=40):
    """Spearman rho per gene across observed lines (Pearson on masked average ranks)."""
    P, T, obs = _prep(pred, true)
    return _masked_pearson_cols(_masked_ranks(P, obs), _masked_ranks(T, obs), obs, min_obs)


def per_line_pearson(pred, true, min_obs=2):
    """Pearson r per line (row) across observed genes -> array (n_lines,)."""
    P, T, obs = _prep(pred, true)
    return _masked_pearson_cols(P.T, T.T, obs.T, min_obs)


def per_line_ap(pred, true, thr=-0.5):
    """Per-line AUPRC across observed genes (P2a). Score = -pred, label = true < thr.

    Lines with no positive (or no negative) observed gene get NaN.
    """
    P, T, obs = _prep(pred, true)
    out = np.full(P.shape[0], np.nan)
    for i in range(P.shape[0]):
        o = obs[i]
        y = T[i, o] < thr
        if y.any() and not y.all():
            out[i] = average_precision_score(y, -P[i, o])
    return out


def per_gene_ap_norm(pred, true, thr=-0.5, min_pos=3, min_neg=3, min_obs=0):
    """Per-gene essentiality ranking across lines (P2b).

    Returns dict of arrays (n_genes,): ap, nap = (AP - π)/(1 - π), auroc, prev (π), n_pos, n_obs.
    ap/nap/auroc are NaN unless the gene has ≥min_pos positives, ≥min_neg negatives and
    ≥min_obs observed lines. prev is the observed-line prevalence (NaN if no observations).
    """
    P, T, obs = _prep(pred, true)
    G = P.shape[1]
    lab = (T < thr) & obs
    n_obs = obs.sum(0)
    n_pos = lab.sum(0)
    n_neg = n_obs - n_pos
    with np.errstate(invalid="ignore", divide="ignore"):
        prev = n_pos / n_obs
    ap = np.full(G, np.nan)
    auroc = np.full(G, np.nan)
    ok = (n_pos >= min_pos) & (n_neg >= min_neg) & (n_obs >= min_obs)
    for g in np.flatnonzero(ok):
        o = obs[:, g]
        y = lab[o, g]
        s = -P[o, g]
        ap[g] = average_precision_score(y, s)
        auroc[g] = roc_auc_score(y, s)
    with np.errstate(invalid="ignore", divide="ignore"):
        nap = (ap - prev) / (1.0 - prev)
    return {"ap": ap, "nap": nap, "auroc": auroc, "prev": np.where(n_obs > 0, prev, np.nan),
            "n_pos": n_pos, "n_obs": n_obs}


def recall_at_fpr(scores, positives, negatives, fpr=0.05):
    """Recall of ``positives`` at a fixed false-positive rate on ``negatives`` (DESIGN §7.5).

    scores: 1-D, higher = more essential. positives/negatives: boolean masks (same length),
    NaN scores are excluded from both. The threshold is the (k+1)-th highest negative score with
    k = floor(fpr * n_neg); a gene is called iff score > threshold (strict), so at most k
    negatives are called (ties at the threshold are not called). Returns NaN if either set is empty.
    """
    s = np.asarray(scores, dtype=np.float64)
    pos = np.asarray(positives, bool) & np.isfinite(s)
    neg = np.asarray(negatives, bool) & np.isfinite(s)
    if pos.sum() == 0 or neg.sum() == 0:
        return np.nan
    ns = np.sort(s[neg])[::-1]
    k = int(np.floor(fpr * len(ns)))
    if k >= len(ns):
        return 1.0
    thr = ns[k]
    return float(np.mean(s[pos] > thr))


def deviation_corr(pred, true, gene_mean):
    """corr(pred - gene_mean, true - gene_mean) across observed genes, per line (DESIGN §7.5).

    pred/true: (n_genes,) or (n_lines, n_genes); gene_mean: (n_genes,). Returns a float for 1-D
    input, else an array (n_lines,). A constant deviation (e.g. B0, deviation ≡ 0) gives 0.
    """
    P = np.asarray(pred, dtype=np.float64)
    T = np.asarray(true, dtype=np.float64)
    m = np.asarray(gene_mean, dtype=np.float64)
    one = P.ndim == 1
    P2, T2 = np.atleast_2d(P), np.atleast_2d(T)
    r = per_line_pearson(P2 - m, T2 - m, min_obs=2)
    return float(r[0]) if one else r
