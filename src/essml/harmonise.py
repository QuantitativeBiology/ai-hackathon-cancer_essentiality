"""Screen harmonisation for external CRISPR screens (DESIGN v3 §7.1-7.2; frozen at tag freeze-transfer).

Pure functions only: no file I/O and no target access (the pipeline driver is
``scripts/harmonise_screens.py``, which routes every DepMap and set-B read through ``essml.guard``).

Steps
-----
(i)   counts -> guide LFC -> gene score (``guide_lfc``, ``gene_scores``):
      RPM-normalise each sample (within library pool if pools were sequenced separately), add a
      pseudocount in read units, drop guides whose mean T0/plasmid reads are below ``min_t0_reads``
      (or ``min_t0_rpm`` for tables deposited already normalised), LFC per end replicate vs the mean
      T0, replicates averaged per guide (mean of replicate LFCs), gene score = MEDIAN guide LFC
      (DESIGN v3 default) over genes with >= ``min_guides`` passing guides.
(ii)  per-screen affine scaling (``affine_params`` / ``affine_apply``): median of the scaling set
      -> t_S and median of the NEGv1 scaling half -> t_N. Per PI decision D10 the targets are DepMap-matched
      (t_S ~ -1.468, t_N ~ 0.007: medians of the same sets in DepMap dev lines, computed by the driver);
      DESIGN's literal -1 / 0 is the function default and is reported as the alternative.
(iii) gene-wise empirical-Bayes location (and, with >= 8 anchors, scale) correction per lab/library
      (``fit_genewise`` / ``apply_genewise``), fitted on TRANSFORMED anchor lines against their DepMap
      dev profiles, after step (ii). p53/DDR genes are excluded from the fit and never corrected.
      Model (location): d_ag = x_ag - y_ag = delta_g + e_ag, delta_g ~ N(0, tau2) (shrinkage toward
      NO correction), e_ag ~ N(0, sigma2_a). tau2 = mean over anchor pairs of the uncentred
      cross-anchor covariance of d (needs >= 2 anchors; with one anchor tau2 is not identifiable and
      the correction is the identity), sigma2_a = mean(d_a^2) - tau2. Posterior mean with the noise
      variance inflated by ``kappa`` (strong shrinkage; kappa = 1 is plain EB):
          delta_hat_g = tau2 * sum_a d_ag / (kappa sigma2_a)  /  (1 + tau2 * sum_a 1 / (kappa sigma2_a)).
      Scale (only if n_anchors >= ``min_anchors_scale``): log SD ratio across anchors
      l_g = log sd_a(x_ag - delta_hat_g) - log sd_a(y_ag), shrunk toward 0 with the same EB recipe
      (per-gene sampling variance by leave-one-anchor-out jackknife, since x and y share anchors),
      and x' = c_g + (x - delta_hat_g - c_g) / exp(l_hat_g), with
      c_g = mean_a y_ag.
      ``loao`` evaluates the step leaving each anchor out (the left-out anchor never enters the fit).
(iv)  QC (``qc_metrics``) on the held-out EVALUATION halves of CEGv2 / NEGv1: medians, NNMD
      (median difference / MAD of NEG), SSMD (mean difference / sqrt(var sum)), AUROC.

Agreement metrics (anchors, replicate ceilings): ``deviation_corr`` = Pearson correlation of
deviations from the DepMap dev gene mean; ``nn_rank`` = rank of the true line among all DepMap dev
lines by deviation correlation.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

# ---------------------------------------------------------------- defaults (draft; frozen in configs/harmonisation_frozen.json)
DEFAULTS = {
    "pseudocount_reads": 0.5,
    "min_t0_reads": 30,          # mean raw T0/plasmid reads per guide (Hart-lab/BAGEL convention)
    "min_t0_rpm": 5.0,           # same filter for tables deposited as RPM / normalised counts
    "min_guides": 2,             # passing guides per gene for a gene score
    "gene_score": "median",      # DESIGN v3 default
    "replicates": "mean_of_replicate_lfc",
    "kappa": 2.0,                # strong shrinkage: noise variance inflated 2x in the EB posterior
    "min_anchors_location": 2,   # tau2 identifiable only from >= 2 anchors
    "min_anchors_scale": 8,      # DESIGN v3 §7.2: scale correction only with >= 8 anchors
    "improvement_margin": 0.01,  # LOAO mean gain in deviation correlation required to apply step (iii)
    "qc_nnmd_max": -1.0,         # DESIGN v3 §7.6 (Avana-type libraries)
    "qc_ssmd_max": -1.0,         # control-gap SSMD for non-Avana libraries (AUROC ~0.84 under normality)
    "qc_min_controls": 10,       # QC is 'not assessable' with < 10 evaluation-half CEG or NEG genes scored
    "affine_target_mode": "depmap_matched",  # PI decision D10; the driver sets the two targets below from dev lines
    "affine_target_scaling": -1.0,  # function default; DepMap-matched value ~ -1.468 is set by the driver
    "affine_target_neg": 0.0,       # function default; DepMap-matched value ~ 0.007
}

# DESIGN v3 §7.2: p53 / DNA-damage-response genes excluded from fitting and application of step (iii)
P53_DDR_CORE = ["TP53", "CDKN1A", "TP53BP1", "USP28", "CHEK2", "ATM", "MDM2", "MDM4", "PPM1D"]


def p53_ddr_symbols():
    """DESIGN §7.2 list: core p53/DDR genes + the p53-target score genes (features.P53_TARGETS)."""
    from .features import P53_TARGETS
    return sorted(set(P53_DDR_CORE) | set(P53_TARGETS))


# ---------------------------------------------------------------- (i) counts -> guide LFC -> gene score
def rpm(counts: pd.DataFrame, pseudocount=0.0, groups=None) -> pd.DataFrame:
    """Reads per million per column (+pseudocount in read units); ``groups`` = per-row pool labels,
    normalised within pool (pools sequenced separately)."""
    c = counts.astype(float)
    if groups is None:
        return (c + pseudocount).div(c.sum(axis=0), axis=1) * 1e6
    out = pd.DataFrame(index=c.index, columns=c.columns, dtype=float)
    g = pd.Series(np.asarray(groups), index=c.index)
    for k in pd.unique(g):
        m = (g == k).to_numpy()
        sub = c.loc[m]
        out.loc[m] = ((sub + pseudocount).div(sub.sum(axis=0), axis=1) * 1e6).to_numpy()
    return out


def guide_lfc(t0: pd.DataFrame, end: pd.DataFrame, normalised=False, pseudocount=None,
              min_t0=None, groups=None):
    """Guide LFC of the end replicates vs the mean T0 (or plasmid).

    t0, end: guides x samples (raw reads, or RPM if ``normalised``). Returns (lfc Series over kept
    guides, report dict). Replicates: LFC per end column, then the mean across end columns.
    """
    pc = DEFAULTS["pseudocount_reads"] if pseudocount is None else pseudocount
    t0 = t0.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    end = end.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    if min_t0 is None:
        min_t0 = DEFAULTS["min_t0_rpm"] if normalised else DEFAULTS["min_t0_reads"]
    keep = t0.mean(axis=1) >= min_t0
    if normalised:  # already per-million: pseudocount in RPM units, no re-normalisation
        r0, re = t0 + pc, end + pc
    else:
        r0, re = rpm(t0, pc, groups), rpm(end, pc, groups)
    l0 = np.log2(r0.mean(axis=1))
    reps = np.log2(re).sub(l0, axis=0)
    lfc = reps.mean(axis=1)[keep]
    rep_r = np.nan
    if reps.shape[1] >= 2:
        cc = reps[keep].corr().to_numpy()
        rep_r = float(np.nanmean(cc[np.triu_indices_from(cc, 1)]))
    rep = {"n_guides": int(len(t0)), "n_guides_kept": int(keep.sum()),
           "t0_filter": f"mean T0 {'RPM' if normalised else 'reads'} >= {min_t0}",
           "pseudocount": f"{pc} {'RPM' if normalised else 'reads'}",
           "n_end_reps": int(end.shape[1]), "replicate_guide_lfc_r": rep_r}
    return lfc, rep


def gene_scores(values: pd.Series, genes: pd.Series, method=None, min_guides=None) -> pd.Series:
    """Gene score from guide values (index-aligned with ``genes``; NaN genes dropped)."""
    method = method or DEFAULTS["gene_score"]
    min_guides = DEFAULTS["min_guides"] if min_guides is None else min_guides
    g = genes.reindex(values.index)
    v = values[g.notna() & values.notna()]
    grp = v.groupby(g[v.index])
    agg = grp.median() if method == "median" else grp.mean()
    return agg[grp.size() >= min_guides].astype(float)


def rank_scores(scores: pd.Series, essential_is="low") -> pd.Series:
    """Rank-only treatment (gene-level BF / MAGeCK / ORCS): percentile in (0, 1], low = essential."""
    s = scores.dropna().astype(float)
    if essential_is == "high":
        s = -s
    elif essential_is != "low":
        raise ValueError("direction unknown")
    return s.rank(method="average") / len(s)


# ---------------------------------------------------------------- control split
def control_split(ceg: set, neg: set, common_essentials: set, seed: int) -> dict:
    """Fixed-seed random halves of CEGv2 and NEGv1 (sorted Entrez strings -> permutation).

    scaling set = DepMap dev common essentials minus the CEGv2 evaluation half (and minus NEGv1).
    """
    rng = np.random.default_rng(seed)
    out = {}
    for name, s in (("ceg", ceg), ("neg", neg)):
        genes = np.array(sorted(map(str, s)))
        perm = rng.permutation(len(genes))
        half = len(genes) // 2
        out[f"{name}_scaling"] = sorted(genes[perm[:half]].tolist())
        out[f"{name}_evaluation"] = sorted(genes[perm[half:]].tolist())
    sc = set(map(str, common_essentials)) - set(out["ceg_evaluation"]) - set(map(str, neg))
    out["scaling_set"] = sorted(sc)
    return out


# ---------------------------------------------------------------- (ii) affine scaling
def affine_params(scores: pd.Series, scaling_set, neg_scaling):
    """(m_S, m_N): medians of the scaling set and the NEG scaling half over scored genes."""
    s = scores.dropna()
    m_s = float(s[s.index.isin(list(scaling_set))].median())
    m_n = float(s[s.index.isin(list(neg_scaling))].median())
    if not np.isfinite(m_s) or not np.isfinite(m_n) or m_n == m_s:
        raise ValueError("affine scaling undefined (missing controls or zero gap)")
    return m_s, m_n


def affine_apply(scores, m_s, m_n, target_s=None, target_n=None):
    """Affine map with m_S -> target_s (default -1) and m_N -> target_n (default 0)."""
    t_s = DEFAULTS["affine_target_scaling"] if target_s is None else target_s
    t_n = DEFAULTS["affine_target_neg"] if target_n is None else target_n
    return t_n + (scores - m_n) * (t_n - t_s) / (m_n - m_s)


# ---------------------------------------------------------------- (iv) QC
def qc_metrics(scores: pd.Series, ceg_eval, neg_eval) -> dict:
    s = scores.dropna().astype(float)
    e = s[s.index.isin(list(ceg_eval))].to_numpy()
    n = s[s.index.isin(list(neg_eval))].to_numpy()
    out = {"n_genes": int(len(s)), "n_ceg_eval": int(len(e)), "n_neg_eval": int(len(n)),
           "median_ceg_eval": np.nan, "median_neg_eval": np.nan, "nnmd": np.nan, "ssmd": np.nan,
           "auroc": np.nan}
    if len(e) < 2 or len(n) < 2:
        return out
    mad = 1.4826 * np.median(np.abs(n - np.median(n)))
    out.update(median_ceg_eval=float(np.median(e)), median_neg_eval=float(np.median(n)),
               nnmd=float((np.median(e) - np.median(n)) / mad) if mad > 0 else np.nan,
               ssmd=float((e.mean() - n.mean()) / np.sqrt(e.var(ddof=1) + n.var(ddof=1))),
               auroc=_auroc_low(e, n))
    return out


def _auroc_low(e, n):
    """P(essential score < non-essential score), ties 1/2."""
    r = pd.Series(np.concatenate([e, n])).rank().to_numpy()
    u = r[: len(e)].sum() - len(e) * (len(e) + 1) / 2
    return float(1 - u / (len(e) * len(n)))


def qc_pass(qc: dict, avana_like=False, nnmd_max=None, ssmd_max=None) -> bool:
    """DESIGN v3 §7.6: NNMD < -1.0 for Avana-type libraries, control-gap SSMD otherwise."""
    if min(qc.get("n_ceg_eval", 0), qc.get("n_neg_eval", 0)) < DEFAULTS["qc_min_controls"]:
        return False  # not assessable -> does not pass
    nnmd_max = DEFAULTS["qc_nnmd_max"] if nnmd_max is None else nnmd_max
    ssmd_max = DEFAULTS["qc_ssmd_max"] if ssmd_max is None else ssmd_max
    v = qc["nnmd"] if avana_like else qc["ssmd"]
    thr = nnmd_max if avana_like else ssmd_max
    return bool(np.isfinite(v) and v < thr)


# ---------------------------------------------------------------- (iii) gene-wise EB correction
@dataclass
class GenewiseCorrection:
    delta: pd.Series                      # location offsets (affine-scaled units), 0 = none
    log_scale: pd.Series                  # log SD ratio (0 = none)
    center: pd.Series                     # c_g = mean DepMap value over anchors (scale pivot)
    anchors: list = field(default_factory=list)
    tau2: float = 0.0
    sigma2: dict = field(default_factory=dict)
    tau2_scale: float = 0.0
    estimable: bool = False
    scale: bool = False
    excluded: list = field(default_factory=list)
    kappa: float = 2.0


def _eb_location(D: np.ndarray, kappa: float):
    """D: anchors x genes (NaN allowed). Returns (delta_hat, tau2, sigma2_a array, estimable)."""
    A = D.shape[0]
    G = D.shape[1]
    if A < DEFAULTS["min_anchors_location"]:
        return np.zeros(G), 0.0, np.full(A, np.nan), False
    cov = []
    for i in range(A):
        for j in range(i + 1, A):
            m = np.isfinite(D[i]) & np.isfinite(D[j])
            if m.sum() > 10:
                cov.append(np.mean(D[i, m] * D[j, m]))
    tau2 = max(0.0, float(np.mean(cov))) if cov else 0.0
    sig2 = np.array([max(np.nanmean(D[a] ** 2) - tau2, 1e-12) for a in range(A)])
    w = np.where(np.isfinite(D), 1.0 / (kappa * sig2[:, None]), 0.0)
    num = tau2 * np.nansum(np.where(np.isfinite(D), D, 0.0) * w, axis=0)
    den = 1.0 + tau2 * w.sum(axis=0)
    return num / den, tau2, sig2, True


def fit_genewise(X: pd.DataFrame, Y: pd.DataFrame, exclude=(), kappa=None, scale=False,
                 min_anchors_scale=None) -> GenewiseCorrection:
    """Fit the step-(iii) correction on anchors.

    X: anchors x genes, external profiles AFTER affine scaling; Y: same anchors x genes, DepMap dev
    profiles (rows matched by position). ``exclude`` genes (p53/DDR) are dropped from the fit and get
    no correction. ``scale`` requests the gene-wise scale term; it is used only with
    >= ``min_anchors_scale`` anchors.
    """
    kappa = DEFAULTS["kappa"] if kappa is None else kappa
    min_anchors_scale = DEFAULTS["min_anchors_scale"] if min_anchors_scale is None else min_anchors_scale
    if list(X.index) != list(Y.index):
        raise ValueError("X and Y must list the same anchors in the same order")
    genes = X.columns.intersection(Y.columns)
    excl = sorted(set(map(str, exclude)) & set(genes))
    fit_genes = genes.difference(excl)
    Xg = X[fit_genes].to_numpy(float)
    Yg = Y[fit_genes].to_numpy(float)
    D = Xg - Yg
    delta, tau2, sig2, ok = _eb_location(D, kappa)
    delta = np.where(np.isfinite(D).any(axis=0), delta, 0.0)
    A = X.shape[0]
    corr = GenewiseCorrection(
        delta=pd.Series(0.0, index=genes), log_scale=pd.Series(0.0, index=genes),
        center=pd.Series(np.nanmean(Y[genes].to_numpy(float), axis=0) if A else 0.0, index=genes),
        anchors=list(X.index), tau2=tau2, sigma2=dict(zip(X.index, map(float, sig2))),
        estimable=ok, excluded=excl, kappa=kappa)
    corr.delta.loc[fit_genes] = delta
    if scale and ok and A >= min_anchors_scale:
        Xc = Xg - delta[None, :]
        ell = _log_sd_ratio(Xc, Yg)
        # paired data: sd_x and sd_y share the anchors, so the sampling variance of l_g is estimated
        # per gene by the leave-one-anchor-out jackknife, not from the independent-SD formula
        jk = np.array([_log_sd_ratio(np.delete(Xc, i, 0), np.delete(Yg, i, 0)) for i in range(A)])
        with np.errstate(invalid="ignore"):
            v = (A - 1) / A * np.nansum((jk - np.nanmean(jk, axis=0)) ** 2, axis=0)
        good = np.isfinite(ell) & np.isfinite(v) & (v > 0)
        t2 = max(0.0, float(np.mean(ell[good] ** 2) - np.mean(v[good]))) if good.any() else 0.0
        ell_hat = np.zeros_like(ell)
        if t2 > 0:
            ell_hat[good] = t2 * ell[good] / (t2 + kappa * v[good])
        corr.log_scale.loc[fit_genes] = ell_hat
        corr.tau2_scale, corr.scale = t2, True
    return corr


def _log_sd_ratio(Xc, Yg):
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.log(np.nanstd(Xc, axis=0, ddof=1)) - np.log(np.nanstd(Yg, axis=0, ddof=1))


def apply_genewise(x, corr: GenewiseCorrection):
    """Apply a fitted correction to a profile (Series) or profiles (DataFrame rows); genes without a
    fitted term (incl. excluded p53/DDR genes) are returned unchanged."""
    if isinstance(x, pd.DataFrame):
        return x.apply(lambda r: apply_genewise(r, corr), axis=1)
    out = x.astype(float).copy()
    g = out.index.intersection(corr.delta.index)
    d, ls, c = corr.delta[g], corr.log_scale[g], corr.center[g]
    v = out[g] - d
    if corr.scale:
        v = c + (v - c) / np.exp(ls)
    out.loc[g] = v
    return out


def loao(X: pd.DataFrame, Y: pd.DataFrame, metric, exclude=(), fit_kwargs=None, targets=None):
    """Leave-one-anchor-out: for each anchor a (or each in ``targets``), fit on X.drop(a) and apply
    to a. ``metric(corrected_profile, anchor_name) -> dict``. Returns (DataFrame of metrics,
    {anchor: corrected profile}, {anchor: fitted correction}). The left-out anchor never enters the
    fit (each correction's ``anchors`` attribute records the fit set)."""
    fit_kwargs = fit_kwargs or {}
    rows, prof, fits = [], {}, {}
    for a in (targets or list(X.index)):
        keep = [i for i in X.index if i != a]
        c = fit_genewise(X.loc[keep], Y.loc[keep], exclude=exclude, **fit_kwargs)
        assert a not in c.anchors
        xc = apply_genewise(X.loc[a], c)
        prof[a], fits[a] = xc, c
        rows.append({"anchor": a, "n_fit_anchors": len(keep), "estimable": c.estimable,
                     "tau2": c.tau2, **metric(xc, a)})
    return pd.DataFrame(rows), prof, fits


# ---------------------------------------------------------------- agreement metrics
def deviation_corr(x: pd.Series, y: pd.Series, mean: pd.Series, genes=None) -> float:
    """corr(x - mean, y - mean) over genes observed in all three (optionally restricted)."""
    idx = x.dropna().index.intersection(y.dropna().index).intersection(mean.dropna().index)
    if genes is not None:
        idx = idx.intersection(pd.Index(list(genes)))
    if len(idx) < 10:
        return np.nan
    a = x[idx].to_numpy(float) - mean[idx].to_numpy(float)
    b = y[idx].to_numpy(float) - mean[idx].to_numpy(float)
    return float(np.corrcoef(a, b)[0, 1])


def deviation_corr_all(x: pd.Series, Ydev: pd.DataFrame, mean: pd.Series) -> pd.Series:
    """Pairwise-complete corr(x - mean, Y_b - mean) for every row b of Ydev (vectorised)."""
    genes = x.dropna().index.intersection(Ydev.columns).intersection(mean.dropna().index)
    u = x[genes].to_numpy(float) - mean[genes].to_numpy(float)
    V = Ydev[genes].to_numpy(float) - mean[genes].to_numpy(float)[None, :]
    M = np.isfinite(V)
    V0 = np.where(M, V, 0.0)
    U0 = np.where(M, u[None, :], 0.0)
    n = M.sum(1)
    su, sv = U0.sum(1), V0.sum(1)
    suv, suu, svv = (U0 * V0).sum(1), (U0 ** 2).sum(1), (V0 ** 2).sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        cov = suv - su * sv / n
        r = cov / np.sqrt((suu - su ** 2 / n) * (svv - sv ** 2 / n))
    return pd.Series(r, index=Ydev.index)


def nn_rank(x: pd.Series, Ydev: pd.DataFrame, mean: pd.Series, own_id: str) -> dict:
    """Rank (1 = best) of ``own_id`` among all DepMap dev lines by deviation correlation."""
    r = deviation_corr_all(x, Ydev, mean)
    rank = int((r > r[own_id]).sum() + 1)
    n = int(r.notna().sum())
    return {"nn_rank": rank, "nn_recall1": int(rank == 1), "nn_recall5": int(rank <= 5),
            "nn_nauc": 1 - (rank - 1) / max(n - 1, 1), "n_ref_lines": n}


# ---------------------------------------------------------------- replicate-ceiling helpers
def auprc(score_essential_high: np.ndarray, label: np.ndarray) -> float:
    """Average precision (step-wise), label 1 = essential, higher score = more essential."""
    s = np.asarray(score_essential_high, float)
    y = np.asarray(label, bool)
    m = np.isfinite(s)
    s, y = s[m], y[m]
    if y.sum() == 0:
        return np.nan
    o = np.argsort(-s, kind="mergesort")
    y = y[o]
    tp = np.cumsum(y)
    prec = tp / np.arange(1, len(y) + 1)
    return float(prec[y].sum() / y.sum())


def recall_at_fpr(score_essential_high, positives, negatives, fpr=0.05) -> float:
    """Recall of ``positives`` at the threshold giving ``fpr`` on ``negatives`` (higher = essential)."""
    p = np.asarray(score_essential_high[positives], float)
    n = np.asarray(score_essential_high[negatives], float)
    p, n = p[np.isfinite(p)], n[np.isfinite(n)]
    if len(p) == 0 or len(n) == 0:
        return np.nan
    thr = np.quantile(n, 1 - fpr)
    return float(np.mean(p > thr))
