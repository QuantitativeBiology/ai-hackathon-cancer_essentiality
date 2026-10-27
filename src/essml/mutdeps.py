"""Mutation-driven dependency diagnostic (DESIGN.md §4 "Mutation-driven dependencies", D24).

Diagnostic only: it does not change P1, the selection rule, delta, gene classes, labels or the
frozen metrics. Everything here works on dev cross-validation folds.

Biomarker status (configs/oncogene_pairs.csv, configs/mutation_dependencies.json)
--------------------------------------------------------------------------------
Status per line is 1 (biomarker-positive, "mutant"), 0 (WT) or NaN (unknown), from the RAW
feature matrices with the thresholds of ``essml.features``:

* ``hotspot``  : hotspot matrix value > 0 (``features.hotspot_binary``);
* ``damaging`` : damaging-mutation dosage > 0 (``features`` block M);
* ``amp``      : relative CN > ``features.AMP_THR`` (2.0);
* ``deep_del`` : relative CN < ``features.DEEP_DEL_THR`` (0.25, = ``priors.DEEP_DEL_THR``);
* ``TP53-WT``  : "mutant" means TP53 wild type = no TP53 damaging and no TP53 hotspot call
  (``features.build_composites`` TP53_mut == 0).

Combined types ("a|b") are positive if any component is positive. A NaN CN of one gene in a line
that has a CN profile is diploid (project convention). A line without the modality a component
needs (e.g. no CN profile) has that component unknown: the status is 1 if another known
component is positive, else NaN. Unknown-status lines are left out of the mutant/WT strata and
of the biomarker AUROC, but stay in the all-line metrics.

Metrics (per pair, model and outer fold; held-out lines of the fold, observed targets only)
------------------------------------------------------------------------------------------
r (frozen ``metrics.per_gene_pearson``), AP and normalised AP for effect < -0.5 (frozen
``metrics.per_gene_ap_norm``), over all lines and within each stratum, and the biomarker
AUROC: among dependent lines (effect < -0.5) with known status, AUROC of -prediction for
mutant vs WT.

Biomarker-only reference: training-fold mean effect of the gene per status (observed training
lines with known status), applied to the held-out lines by their status; unknown status or an
empty training stratum -> the training-fold gene mean.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

from essml import features as F
from essml import metrics

THR = -0.5
TYPES = ("hotspot", "damaging", "amp", "deep_del", "TP53-WT")

# DESIGN §4 pair set: (pair_id, group, dependency symbol, biomarker symbol, biomarker type)
PAIRS = [
    ("BRAF_hotspot", "oncogene", "BRAF", "BRAF", "hotspot"),
    ("KRAS_hotspot", "oncogene", "KRAS", "KRAS", "hotspot"),
    ("NRAS_hotspot", "oncogene", "NRAS", "NRAS", "hotspot"),
    ("PIK3CA_hotspot", "oncogene", "PIK3CA", "PIK3CA", "hotspot"),
    ("EGFR_hotspot_amp", "oncogene", "EGFR", "EGFR", "hotspot|amp"),
    ("ERBB2_amp", "oncogene", "ERBB2", "ERBB2", "amp"),
    ("CTNNB1_hotspot", "oncogene", "CTNNB1", "CTNNB1", "hotspot"),
    ("CCND1_amp", "oncogene", "CCND1", "CCND1", "amp"),
    ("MDM2_TP53WT", "oncogene", "MDM2", "TP53", "TP53-WT"),
    ("MYC_amp", "oncogene", "MYC", "MYC", "amp"),
    ("ARID1A_ARID1B", "paralog", "ARID1B", "ARID1A", "damaging|deep_del"),
    ("SMARCA4_SMARCA2", "paralog", "SMARCA2", "SMARCA4", "damaging|deep_del"),
    ("STAG2_STAG1", "paralog", "STAG1", "STAG2", "damaging|deep_del"),
    ("ENO1_ENO2", "paralog", "ENO2", "ENO1", "damaging|deep_del"),
    ("VPS4B_VPS4A", "paralog", "VPS4A", "VPS4B", "damaging|deep_del"),
    ("MTAP_PRMT5", "paralog", "PRMT5", "MTAP", "damaging|deep_del"),
    ("MTAP_MAT2A", "paralog", "MAT2A", "MTAP", "damaging|deep_del"),
    ("EP300_CREBBP", "paralog", "CREBBP", "EP300", "damaging|deep_del"),
    ("CREBBP_EP300", "paralog", "EP300", "CREBBP", "damaging|deep_del"),
]


# ============================================================================ status
def _component(kind, entrez, ids, hot, dam, cn):
    """(positive 0/1 array, known bool array) for one status component over ``ids``."""
    ids = pd.Index(ids)
    if kind in ("hotspot", "damaging"):
        mat = hot if kind == "hotspot" else dam
        known = ids.isin(mat.index)
        v = (mat[entrez].reindex(ids).fillna(0.0).to_numpy() if entrez in mat.columns
             else np.zeros(len(ids)))
        return (v > 0) & known, known
    if kind in ("amp", "deep_del"):
        known = ids.isin(cn.index)
        v = (cn[entrez].reindex(ids).to_numpy(np.float64) if entrez in cn.columns
             else np.ones(len(ids)))
        v = np.where(np.isnan(v), 1.0, v)  # NaN CN of a gene in a profiled line -> diploid
        pos = v > F.AMP_THR if kind == "amp" else v < F.DEEP_DEL_THR
        return pos & known, known
    raise ValueError(f"unknown component {kind}")


def biomarker_status(btype, entrez, ids, hot, dam, cn) -> pd.Series:
    """1 = biomarker-positive ("mutant"), 0 = WT, NaN = unknown (module docstring)."""
    ids = pd.Index([str(i) for i in ids])
    if btype == "TP53-WT":
        h, kh = _component("hotspot", entrez, ids, hot, dam, cn)
        d, kd = _component("damaging", entrez, ids, hot, dam, cn)
        known = kh & kd
        out = np.where(known, (~h & ~d).astype(float), np.nan)
        return pd.Series(out, index=ids, name="status")
    parts = btype.split("|")
    if any(p not in TYPES for p in parts):
        raise ValueError(f"unknown biomarker type {btype}")
    pos = np.zeros(len(ids), bool)
    all_known = np.ones(len(ids), bool)
    for p in parts:
        v, k = _component(p, entrez, ids, hot, dam, cn)
        pos |= v
        all_known &= k
    out = np.where(pos, 1.0, np.where(all_known, 0.0, np.nan))
    return pd.Series(out, index=ids, name="status")


# ============================================================================ metrics
def _r(p, t, min_obs):
    if len(t) == 0:
        return np.nan
    return float(metrics.per_gene_pearson(p[:, None], t[:, None], min_obs=min_obs)[0])


def _ap(p, t, min_pos, min_neg):
    if len(t) == 0:
        return np.nan, np.nan
    d = metrics.per_gene_ap_norm(p[:, None], t[:, None], thr=THR, min_pos=min_pos, min_neg=min_neg)
    return float(d["ap"][0]), float(d["nap"][0])


def biomarker_auroc(pred, true, status, min_each=3, thr=THR):
    """AUROC of -pred for mutant (1) vs WT (0) among dependent (true < thr) known-status lines."""
    pred, true, status = (np.asarray(a, np.float64) for a in (pred, true, status))
    m = ~np.isnan(true) & (true < thr) & ~np.isnan(status)
    y = status[m] == 1
    if y.sum() < min_each or (~y).sum() < min_each:
        return np.nan
    return float(roc_auc_score(y, -pred[m]))


def fold_metrics(pred, true, status, min_obs_all=40, min_obs_stratum=5, min_pos=3, min_neg=3,
                 min_auroc=3) -> dict:
    """All pair metrics for one model in one fold (1-D arrays over the held-out lines)."""
    pred, true, status = (np.asarray(a, np.float64) for a in (pred, true, status))
    obs = ~np.isnan(true)
    out = {"n_obs": int(obs.sum()), "n_dep": int((true[obs] < THR).sum())}
    p, t, s = pred[obs], true[obs], status[obs]
    out["r_all"] = _r(p, t, min_obs_all)
    out["ap_all"], out["nap_all"] = _ap(p, t, min_pos, min_neg)
    for name, val in (("mut", 1.0), ("wt", 0.0)):
        m = s == val
        out[f"n_{name}"] = int(m.sum())
        out[f"n_dep_{name}"] = int((t[m] < THR).sum())
        out[f"r_{name}"] = _r(p[m], t[m], min_obs_stratum)
        out[f"ap_{name}"], out[f"nap_{name}"] = _ap(p[m], t[m], min_pos, min_neg)
    out["auroc_bm"] = biomarker_auroc(p, t, s, min_auroc)
    return out


def biomarker_only(true_train, status_train, status_test) -> np.ndarray:
    """Training-fold mean effect by status, applied to held-out status (module docstring)."""
    t = np.asarray(true_train, np.float64)
    s = np.asarray(status_train, np.float64)
    st = np.asarray(status_test, np.float64)
    obs = ~np.isnan(t)
    mu = float(t[obs].mean())
    out = np.full(len(st), mu)
    for val in (0.0, 1.0):
        m = obs & (s == val)
        if m.any():
            out[st == val] = float(t[m].mean())
    return out


# ============================================================================ rule
def paired_ci(d) -> dict:
    """Mean of fold-wise paired differences with a t-based 95% CI (as scripts/g2_report.paired)."""
    d = np.asarray(d, np.float64)
    d = d[np.isfinite(d)]
    n = len(d)
    m = float(d.mean()) if n else np.nan
    if n < 2:
        return {"mean": m, "se": np.nan, "ci95": [np.nan, np.nan], "n_folds": n,
                "n_folds_pos": int((d > 0).sum())}
    se = float(d.std(ddof=1) / np.sqrt(n))
    h = float(stats.t.ppf(0.975, n - 1) * se)
    return {"mean": m, "se": se, "ci95": [m - h, m + h], "n_folds": n,
            "n_folds_pos": int((d > 0).sum())}


def clearly_helps(delta_ap: dict, eligible, gain=0.10, min_folds=4, min_pairs=3, n_folds=5) -> dict:
    """DESIGN §4 "genomics clearly helps" for one comparator.

    delta_ap: pair -> array (n_folds,) of AP(comparator) - AP(m1_E_resB0), NaN where undefined.
    Condition 1: mean over the eligible pairs whose AP is defined in every fold, fold by fold, then
    a t-based 95% CI over the folds with lower bound > 0. Condition 2: >= min_pairs eligible pairs
    with delta >= gain in >= min_folds folds (an undefined fold counts as no gain).
    """
    elig = [p for p in eligible if p in delta_ap]
    full = [p for p in elig if np.isfinite(np.asarray(delta_ap[p], float)).all()
            and len(delta_ap[p]) == n_folds]
    if full:
        pooled = np.mean(np.stack([np.asarray(delta_ap[p], float) for p in full]), axis=0)
        ci = paired_ci(pooled)
        ci["folds"] = pooled.round(5).tolist()
    else:
        ci = {"mean": np.nan, "se": np.nan, "ci95": [np.nan, np.nan], "folds": []}
    cond1 = bool(np.isfinite(ci["ci95"][0]) and ci["ci95"][0] > 0)
    n_gain = {p: int(np.nansum(np.asarray(delta_ap[p], float) >= gain)) for p in elig}
    winners = sorted(p for p, k in n_gain.items() if k >= min_folds)
    cond2 = len(winners) >= min_pairs
    return {"pooled": ci, "pooled_pairs": full, "cond1_ci_excludes_0": cond1,
            "n_folds_gain_ge": n_gain, "pairs_gain": winners, "cond2_ge3_pairs": cond2,
            "met": cond1 and cond2}
