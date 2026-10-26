"""Transfer verdicts for non-transformed screens (DESIGN v3 §7.4-7.6 as amended by PI decision D13).

Pure computation: no model code and no file I/O. Every function takes OBSERVED profiles and PREDICTION
vectors as inputs (pandas Series indexed by Entrez strings). Thresholds and populations come from
``configs/transfer_criteria.json`` (written by ``scripts/transfer_criteria.py``); the defaults below
are the frozen values and are asserted equal to the file in tests/test_transfer.py.

Definitions (D13)
-----------------
- Unit of evidence = STUDY (publication). A study has one or more ENTRIES: one consensus profile per
  TP53 status = per-gene mean (NaN-skipping) of its QC-passing genome-wide fitness arms of that status.
  Study value = mean over its entries; background value = mean over its verdict-eligible studies.
- Gene universe U_e = eval genes (configs/eval_genes.csv) observed in the entry consensus, minus the
  scaling-half control genes (ceg_scaling ∪ neg_scaling ∪ scaling_set of configs/control_split.json,
  pinned by the affine step). Evaluation-half CEG/NEG genes stay in.
- Q3a (primary): label = rank label (genes at or below the 10% quantile of the consensus over U_e);
  statistic = AP(model) - AP(B0); AP = sklearn step-wise average precision of score -pred.
- Q3b (primary): positives = consensus < -0.5 AND DepMap-dev dependency fraction <= 0.5 (not NEGv1
  evaluation genes); negatives = NEGv1 evaluation half; statistic = recall at 5% FPR (threshold on the
  negatives, ``metrics.recall_at_fpr`` convention) of the model minus that of B1-nearest.
- D7: with fold predictions {fold: Series}, skill is computed per fold model and the skills averaged.
- Replicate ceiling of entry e (study s): the same statistic with the predictor = leave-one-study-out
  consensus = mean of the OTHER studies' study profiles of the background (descriptive studies may
  serve as predictors); genes unobserved in every other study are filled with B0.
- Hierarchical bootstrap (``BootstrapPlan``): per replicate, (1) resample the verdict-eligible studies
  with replacement, (2) resample paralog-family gene blocks with replacement (gene weights =
  multiplicities, shared by every entry), (3) ceiling only: resample each target's predictor pool with
  replacement. Independent streams from SeedSequence(seed).spawn(3), so the model-skill bootstrap
  and the ceiling bootstrap use identical study and gene draws. CI = percentile 95%;
  MDE = 2 x bootstrap SE of the background ceiling skill.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

# ---------------------------------------------------------------- frozen constants (mirrored in the criteria file)
RANK_Q = 0.10            # Q3a primary rank label: 10% most depleted
EFFECT_THR = -0.5        # effect label (Q3a secondary; Q3b positives)
FRAC_MAX = 0.5           # Q3b screen-specific: DepMap-dev dependency fraction <= 0.5
FRAC_MAX_RARE = 0.05     # Q3b secondary definition
FPR = 0.05               # Q3b threshold on the NEGv1 evaluation half
N_BOOT = 2000
CI_LEVEL = 0.95
MDE_MULT = 2.0
CEILING_FRACTION = 0.5   # transfers: mean skill >= 0.5 x ceiling skill
PARALOG_MIN_IDENTITY = 30.0
MIN_UNIVERSE = 1000      # entry statistic is NaN (entry dropped) below these counts
MIN_Q3B_POS = 10
MIN_Q3B_NEG = 100
H2_MIN_GENES = 10
H2_MIN_SPEARMAN = 0.3
H2_MIN_SIGN = 0.70
H3_MIN_OBSERVED = 3
H3_PASS_FRACTION = 0.5

TRANSFERS, FAILS, INCONCLUSIVE, MIXED = "transfers", "fails", "inconclusive", "mixed"


# ================================================================ gene blocks
def paralog_blocks(genes, pairs: pd.DataFrame, min_identity=PARALOG_MIN_IDENTITY) -> pd.Series:
    """Block id per gene: connected components of paralog pairs with min(identity, rev_identity) >=
    ``min_identity`` (columns a, b, identity, rev_identity); unpaired genes are singleton blocks."""
    g = pd.Index(sorted(set(map(str, genes))))
    p = pairs[np.minimum(pairs["identity"], pairs["rev_identity"]) >= min_identity]
    a = g.get_indexer(p["a"].astype(str))
    b = g.get_indexer(p["b"].astype(str))
    m = (a >= 0) & (b >= 0)
    adj = coo_matrix((np.ones(int(m.sum())), (a[m], b[m])), shape=(len(g), len(g)))
    _, lab = connected_components(adj, directed=False)
    return pd.Series(lab, index=g, name="block")


# ================================================================ weighted metrics (== frozen metrics on duplicated genes)
class APPrep:
    """Step-wise average precision (sklearn ``average_precision_score`` semantics, ties grouped) with
    optional per-gene weights; integer weights equal duplicating genes. Higher score = more essential."""

    def __init__(self, score_high, y):
        s = np.asarray(score_high, dtype=np.float64)
        if not np.all(np.isfinite(s)):
            raise ValueError("non-finite scores")
        self.o = np.argsort(-s, kind="mergesort")
        ss = s[self.o]
        self.y = np.asarray(y, bool)[self.o].astype(np.float64)
        self.ends = np.r_[np.flatnonzero(np.diff(ss) != 0), len(ss) - 1] if len(ss) else np.array([], int)

    def __call__(self, w=None) -> float:
        if len(self.y) == 0:
            return np.nan
        if w is None:
            wy, wn = self.y, 1.0 - self.y
        else:
            ww = np.asarray(w, dtype=np.float64)[self.o]
            wy, wn = ww * self.y, ww * (1.0 - self.y)
        tp = np.cumsum(wy)[self.ends]
        tot = tp + np.cumsum(wn)[self.ends]
        if tp[-1] <= 0:
            return np.nan
        keep = np.r_[tot[0] > 0, np.diff(tot) > 0]  # drop tie groups with zero total weight
        tp, tot = tp[keep], tot[keep]
        rec = tp / tp[-1]
        return float(np.sum(np.diff(np.r_[0.0, rec]) * (tp / tot)))


def average_precision(score_high, y, w=None) -> float:
    return APPrep(score_high, y)(w)


class RecallPrep:
    """Recall of positives at ``fpr`` on negatives, ``metrics.recall_at_fpr`` convention (threshold =
    (k+1)-th highest negative score, k = floor(fpr * n_neg); called iff score > threshold), with
    optional integer gene weights (== duplicating genes)."""

    def __init__(self, score_high, pos, neg, fpr=FPR):
        s = np.asarray(score_high, dtype=np.float64)
        pos, neg = np.asarray(pos, bool), np.asarray(neg, bool)
        if np.any(pos & neg):
            raise ValueError("a gene cannot be both positive and negative")
        if not np.all(np.isfinite(s[pos | neg])):
            raise ValueError("non-finite scores")
        ni = np.flatnonzero(neg)
        o = np.argsort(-s[ni], kind="mergesort")
        self.ni, self.ns = ni[o], s[ni][o]
        self.pi = np.flatnonzero(pos)
        self.ps = s[self.pi]
        self.fpr = fpr

    def __call__(self, w=None) -> float:
        if w is None:
            wn, wp = np.ones(len(self.ni)), np.ones(len(self.pi))
        else:
            w = np.asarray(w, dtype=np.float64)
            wn, wp = w[self.ni], w[self.pi]
        n_neg, n_pos = wn.sum(), wp.sum()
        if n_neg <= 0 or n_pos <= 0:
            return np.nan
        k = int(np.floor(self.fpr * n_neg))
        if k >= n_neg:
            return 1.0
        j = int(np.searchsorted(np.cumsum(wn), k, side="right"))
        thr = self.ns[j]
        return float(wp[self.ps > thr].sum() / n_pos)


# ================================================================ universes, labels, consensus
def scaling_halves(control_split: dict) -> set:
    """Genes pinned by the affine step: ceg_scaling ∪ neg_scaling ∪ scaling_set (D13.2)."""
    return set(control_split["ceg_scaling"]) | set(control_split["neg_scaling"]) | set(control_split["scaling_set"])


def gene_universe(observed: pd.Series, eval_genes, exclude) -> pd.Index:
    """U_e = eval genes observed (finite) in the consensus, minus ``exclude`` (scaling halves); sorted."""
    obs = observed[np.isfinite(observed.astype(float))]
    U = pd.Index(sorted(set(obs.index.astype(str)) & set(map(str, eval_genes)) - set(map(str, exclude))))
    return U


def rank_label(observed: pd.Series, U: pd.Index, q=RANK_Q) -> np.ndarray:
    """True for genes at or below the q-quantile (numpy 'linear') of the consensus over U."""
    x = observed.reindex(U).to_numpy(np.float64)
    return x <= np.quantile(x, q)


def effect_label(observed: pd.Series, U: pd.Index, thr=EFFECT_THR) -> np.ndarray:
    return observed.reindex(U).to_numpy(np.float64) < thr


def q3b_sets(observed: pd.Series, U: pd.Index, frac_below: pd.Series, neg_eval, thr=EFFECT_THR,
             frac_max=FRAC_MAX):
    """(positives, negatives) boolean arrays over U. Genes without a dependency fraction are not positives."""
    x = observed.reindex(U).to_numpy(np.float64)
    neg = U.isin(list(map(str, neg_eval)))
    fr = frac_below.reindex(U).to_numpy(np.float64)
    pos = (x < thr) & np.isfinite(fr) & (fr <= frac_max) & ~neg
    return pos, neg


def consensus(arms) -> pd.Series:
    """Per-gene mean over arms (NaN-skipping; union of genes)."""
    df = pd.concat([a.astype(float) for a in arms], axis=1)
    return df.mean(axis=1, skipna=True).dropna()


# ================================================================ entry evaluator
@dataclass
class Entry:
    name: str
    study: str
    observed: pd.Series
    tp53: str = ""


class EntryEval:
    """One entry x one question: fixes U, labels and the gene positions inside a plan's gene index."""

    def __init__(self, entry: Entry, question: str, eval_genes, exclude, frac_below=None, neg_eval=None,
                 label="rank10", plan_genes: pd.Index | None = None, universe: pd.Index | None = None,
                 frac_max=FRAC_MAX):
        self.entry, self.question, self.label = entry, question, label
        U = gene_universe(entry.observed, eval_genes, exclude) if universe is None else pd.Index(universe)
        if plan_genes is not None:
            U = U[U.isin(plan_genes)]
        self.U = U
        self.idx = None if plan_genes is None else plan_genes.get_indexer(U)
        self.valid = len(U) >= MIN_UNIVERSE
        if question == "Q3a":
            self.y = rank_label(entry.observed, U) if label == "rank10" else effect_label(entry.observed, U)
            self.n_pos = int(self.y.sum())
            self.valid &= self.n_pos > 0
        elif question == "Q3b":
            self.pos, self.neg = q3b_sets(entry.observed, U, frac_below, neg_eval, frac_max=frac_max)
            self.n_pos, self.n_neg = int(self.pos.sum()), int(self.neg.sum())
            self.valid &= (self.n_pos >= MIN_Q3B_POS) & (self.n_neg >= MIN_Q3B_NEG)
        else:
            raise ValueError(question)

    def prep(self, pred: pd.Series):
        p = pred.reindex(self.U).to_numpy(np.float64)
        if not np.all(np.isfinite(p)):
            raise ValueError(f"prediction not finite on the universe of {self.entry.name}")
        s = -p
        return APPrep(s, self.y) if self.question == "Q3a" else RecallPrep(s, self.pos, self.neg)

    def weights(self, w_plan):
        return None if w_plan is None else np.asarray(w_plan)[self.idx]

    def summary(self) -> dict:
        d = {"entry": self.entry.name, "study": self.entry.study, "n_universe": int(len(self.U)),
             "n_pos": self.n_pos, "valid": bool(self.valid)}
        if self.question == "Q3b":
            d["n_neg"] = self.n_neg
        return d


# ================================================================ bootstrap plan
class BootstrapPlan:
    """Deterministic hierarchical resampling plan (studies -> gene blocks; + predictor pools).

    studies: verdict-eligible study names (sorted internally); genes: plan gene index (sorted);
    blocks: Series gene -> block id. ``replicates()`` yields (r, study_index_array, gene_weights).
    """

    def __init__(self, studies, genes, blocks: pd.Series, n_boot=N_BOOT, seed=20261026):
        self.studies = sorted(studies)
        self.genes = pd.Index(sorted(map(str, genes)))
        b = blocks.reindex(self.genes)
        if b.isna().any():  # genes absent from the block table are singletons
            start = int(np.nanmax(b.to_numpy())) + 1 if b.notna().any() else 0
            b[b.isna()] = np.arange(start, start + int(b.isna().sum()))
        codes, uniq = pd.factorize(b.astype(np.int64))
        self.block_of = codes
        self.n_blocks = len(uniq)
        self.n_boot, self.seed = int(n_boot), int(seed)

    def _streams(self):
        return [np.random.default_rng(s) for s in np.random.SeedSequence(self.seed).spawn(3)]

    def replicates(self, with_pool=False):
        rs, rg, rp = self._streams()
        n_s = len(self.studies)
        for r in range(self.n_boot):
            si = rs.integers(0, n_s, n_s)
            bc = np.bincount(rg.integers(0, self.n_blocks, self.n_blocks), minlength=self.n_blocks)
            w = bc[self.block_of].astype(np.float64)
            yield (r, si, w, rp) if with_pool else (r, si, w)


def _summ(point, boot, level=CI_LEVEL) -> dict:
    b = np.asarray(boot, dtype=np.float64)
    ok = np.isfinite(b)
    a = (1 - level) / 2 * 100
    if ok.sum() < 2:
        return {"point": point, "ci_lo": np.nan, "ci_hi": np.nan, "se": np.nan, "n_boot_valid": int(ok.sum())}
    lo, hi = np.percentile(b[ok], [a, 100 - a])
    return {"point": float(point), "ci_lo": float(lo), "ci_hi": float(hi), "se": float(np.std(b[ok], ddof=1)),
            "n_boot_valid": int(ok.sum())}


def _bg_mean(study_vals: dict, studies, draw=None) -> float:
    names = studies if draw is None else [studies[i] for i in draw]
    v = np.array([study_vals.get(s, np.nan) for s in names], dtype=np.float64)
    v = v[np.isfinite(v)]
    return float(v.mean()) if len(v) else np.nan


def _study_vals(entry_vals: dict, evals: dict) -> dict:
    by = {}
    for name, v in entry_vals.items():
        if np.isfinite(v):
            by.setdefault(evals[name].entry.study, []).append(v)
    return {s: float(np.mean(v)) for s, v in by.items()}


def _as_folds(p):
    return p if isinstance(p, dict) else {0: p}


# ================================================================ model skill (point + hierarchical bootstrap)
def skill_bootstrap(evals: dict, preds: dict, comparator: dict, plan: BootstrapPlan, keep_boot=False) -> dict:
    """Background skill of a model over a comparator.

    evals: {entry_name: EntryEval} (built with plan.genes); preds: {entry_name: Series or {fold: Series}}
    (D7: per-fold skill, averaged over folds); comparator: {entry_name: Series}.
    Returns point/CI/SE of the background skill plus per-entry and per-study point values.
    """
    live = {n: e for n, e in evals.items() if e.valid}
    mp = {n: [e.prep(p) for p in _as_folds(preds[n]).values()] for n, e in live.items()}
    cp = {n: e.prep(comparator[n]) for n, e in live.items()}

    def entry_vals(w):
        out = {}
        for n, e in live.items():
            ww = e.weights(w)
            out[n] = float(np.mean([m(ww) for m in mp[n]])) - cp[n](ww)
        return out

    def abs_vals(w, which):
        out = {}
        for n, e in live.items():
            ww = e.weights(w)
            out[n] = float(np.mean([m(ww) for m in mp[n]])) if which == "model" else cp[n](ww)
        return out

    ev = entry_vals(None)
    sv = _study_vals(ev, live)
    point = _bg_mean(sv, plan.studies)
    boot = np.empty(plan.n_boot)
    for r, si, w in plan.replicates():
        boot[r] = _bg_mean(_study_vals(entry_vals(w), live), plan.studies, si)
    out = _summ(point, boot)
    out.update({"per_entry": ev, "per_study": sv,
                "model_abs_per_entry": abs_vals(None, "model"), "comparator_abs_per_entry": abs_vals(None, "comp"),
                "n_entries_valid": len(live), "entries_invalid": sorted(set(evals) - set(live))})
    if keep_boot:
        out["boot"] = boot
    return out


# ================================================================ replicate ceiling (LOSO consensus of other studies)
def study_profiles(entries) -> dict:
    """Study profile = per-gene mean over its entries (used only as a ceiling predictor)."""
    by = {}
    for e in entries:
        by.setdefault(e.study, []).append(e.observed)
    return {s: consensus(v) for s, v in by.items()}


def _pool_predictor(pool_names, counts, profiles, b0: pd.Series, U: pd.Index) -> pd.Series:
    num = np.zeros(len(U))
    den = np.zeros(len(U))
    for nm, c in zip(pool_names, counts):
        if c == 0:
            continue
        x = profiles[nm].reindex(U).to_numpy(np.float64)
        m = np.isfinite(x)
        num[m] += c * x[m]
        den[m] += c
    with np.errstate(invalid="ignore", divide="ignore"):
        p = np.where(den > 0, num / den, b0.reindex(U).to_numpy(np.float64))
    return pd.Series(p, index=U)


def ceiling_bootstrap(evals: dict, pools: dict, profiles: dict, comparator: dict, b0: pd.Series,
                      plan: BootstrapPlan, resample_pool=True, keep_boot=False) -> dict:
    """Replicate-ceiling skill: the predictor of entry e is the mean of the study profiles in
    pools[study(e)] (the OTHER studies of the background), B0 where none observes the gene.

    With ``resample_pool`` each eligible study's pool is resampled with replacement per replicate
    (third stream of the plan); the study and gene draws are those of ``skill_bootstrap``.
    """
    live = {n: e for n, e in evals.items() if e.valid}
    cp = {n: e.prep(comparator[n]) for n, e in live.items()}
    cache = {}

    def pred_prep(n, counts):
        key = (n, tuple(int(c) for c in counts))
        if key not in cache:
            e = live[n]
            pool = pools.get(e.entry.study, [])
            cache[key] = e.prep(_pool_predictor(pool, counts, profiles, b0, e.U))
        return cache[key]

    def entry_vals(w, counts_of):
        out = {}
        for n, e in live.items():
            ww = e.weights(w)
            out[n] = pred_prep(n, counts_of[e.entry.study])(ww) - cp[n](ww)
        return out

    full = {s: np.ones(len(pools.get(s, [])), int) for s in plan.studies}
    ev = entry_vals(None, full)
    sv = _study_vals(ev, live)
    point = _bg_mean(sv, plan.studies)
    boot = np.empty(plan.n_boot)
    for r, si, w, rp in plan.replicates(with_pool=True):
        counts = {}
        for s in plan.studies:  # fixed order -> reproducible stream use
            k = len(pools.get(s, []))
            counts[s] = np.bincount(rp.integers(0, k, k), minlength=k) if (resample_pool and k) else np.ones(k, int)
        boot[r] = _bg_mean(_study_vals(entry_vals(w, counts), live), plan.studies, si)
    out = _summ(point, boot)
    out.update({"per_entry": ev, "per_study": sv,
                "predictor_abs_per_entry": {n: pred_prep(n, full[live[n].entry.study])(None) for n in live},
                "comparator_abs_per_entry": {n: cp[n](None) for n in live},
                "n_entries_valid": len(live), "entries_invalid": sorted(set(evals) - set(live))})
    out["mde"] = MDE_MULT * out["se"] if np.isfinite(out["se"]) else np.nan
    if keep_boot:
        out["boot"] = boot
    return out


# ================================================================ anchor gate (D13.7)
def anchor_gate(anchor_results: dict) -> dict:
    """anchor_results: {anchor_id: skill_bootstrap result of the same statistic on that anchor (single
    OOF fold model, D7)}. Informative and positive iff >= 1 anchor has lower CI > 0."""
    pos = sorted(a for a, r in anchor_results.items() if np.isfinite(r.get("ci_lo", np.nan)) and r["ci_lo"] > 0)
    return {"n_anchors": len(anchor_results), "positive_anchors": pos, "informative_positive": bool(pos)}


# ================================================================ verdicts
def background_verdict(skill: dict, ceiling_skill: float, mde: float, harmonisation_ok: bool,
                       anchor_gate_positive: bool, swap_results=()) -> dict:
    """D13.7 verdict for one background x question.

    transfers: lower CI(skill) > 0 AND point skill >= 0.5 x ceiling skill (D13 literal: with a ceiling
               skill <= 0 the ratio condition is vacuous; this is flagged in the output, not altered);
    fails:     upper CI(skill) < MDE (strict) AND harmonisation validation passed AND the anchor gate is
               informative and positive;
    both conditions met -> inconclusive (rule cannot separate); neither -> inconclusive.
    Expression-swap cap: any swapped-profile skill with lower CI > 0 caps 'transfers' at 'inconclusive'.
    """
    lo, hi, pt = skill.get("ci_lo", np.nan), skill.get("ci_hi", np.nan), skill.get("point", np.nan)
    reasons = []
    if not all(np.isfinite([lo, hi, pt])) or not np.isfinite(ceiling_skill) or not np.isfinite(mde):
        return {"verdict": INCONCLUSIVE, "reasons": ["statistic, ceiling or MDE not finite"],
                "transfers_condition": False, "fails_condition": False, "swap_cap": False}
    t = bool(lo > 0 and pt >= CEILING_FRACTION * ceiling_skill)
    f = bool(hi < mde and harmonisation_ok and anchor_gate_positive)
    cap = any(np.isfinite(s.get("ci_lo", np.nan)) and s["ci_lo"] > 0 for s in swap_results)
    if ceiling_skill <= 0:
        reasons.append("ceiling skill <= 0: the 50%-of-ceiling condition is vacuous")
    if t and f:
        v = INCONCLUSIVE
        reasons.append("transfers and fails conditions both met")
    elif t:
        v = INCONCLUSIVE if cap else TRANSFERS
        if cap:
            reasons.append("expression-swap null has lower CI > 0: capped")
    elif f:
        v = FAILS
    else:
        v = INCONCLUSIVE
        if hi < mde and not harmonisation_ok:
            reasons.append("fails blocked: harmonisation validation")
        if hi < mde and not anchor_gate_positive:
            reasons.append("fails blocked: anchor gate not informative/positive")
    return {"verdict": v, "reasons": reasons, "transfers_condition": t, "fails_condition": f, "swap_cap": bool(cap),
            "ratio_condition_vacuous": bool(ceiling_skill <= 0)}


def overall_verdict(per_background: dict) -> dict:
    """All verdict-eligible backgrounds agree on transfers/fails -> that verdict; both 'transfers' and
    'fails' present -> 'mixed'; otherwise 'inconclusive'. The per-background vector is always returned."""
    vs = list(per_background.values())
    if not vs:
        return {"verdict": INCONCLUSIVE, "per_background": {}}
    if all(v == TRANSFERS for v in vs):
        o = TRANSFERS
    elif all(v == FAILS for v in vs):
        o = FAILS
    elif TRANSFERS in vs and FAILS in vs:
        o = MIXED
    else:
        o = INCONCLUSIVE
    return {"verdict": o, "per_background": dict(per_background)}


def harmonisation_ok(arm_qc_pass, consensus_qc_pass) -> bool:
    """All arms of every contributing study passed QC and every study consensus passes QC."""
    return bool(all(arm_qc_pass) and all(consensus_qc_pass))


# ================================================================ set A (tier A1, descriptive)
def percentile_in(distribution, value) -> float:
    """100 x fraction of the reference (cancer-test per-line) distribution <= value."""
    d = np.asarray(distribution, dtype=np.float64)
    d = d[np.isfinite(d)]
    return float(100.0 * np.mean(d <= value)) if len(d) and np.isfinite(value) else np.nan


def set_a_descriptive(line_results: dict, test_distribution) -> dict:
    """line_results: {line: skill_bootstrap result (gene-block bootstrap)}. Descriptive label from the
    median within-test percentile: >= 25 'comparable to cancer test lines'; 5-25 'lower than most cancer
    test lines'; < 5 'outside the cancer range'. Not verdict-bearing."""
    rows = {}
    for ln, r in line_results.items():
        cat = ("skill>0" if r["ci_lo"] > 0 else "skill<0" if r["ci_hi"] < 0 else "CI covers 0") \
            if np.isfinite(r.get("ci_lo", np.nan)) else "n/a"
        rows[ln] = {"point": r["point"], "ci_lo": r.get("ci_lo"), "ci_hi": r.get("ci_hi"), "category": cat,
                    "test_percentile": percentile_in(test_distribution, r["point"])}
    pc = np.array([v["test_percentile"] for v in rows.values()], dtype=np.float64)
    med = float(np.nanmedian(pc)) if np.isfinite(pc).any() else np.nan
    label = ("n/a" if not np.isfinite(med) else "comparable to cancer test lines" if med >= 25
             else "lower than most cancer test lines" if med >= 5 else "outside the cancer range")
    return {"per_line": rows, "median_test_percentile": med, "descriptive_verdict": label,
            "n_lines_skill_positive": int(sum(v["category"] == "skill>0" for v in rows.values()))}


# ================================================================ hypotheses (secondary, operational)
def h1(model_delta: dict, observed_delta: dict, genes=("MDM2", "PPM1D")) -> dict:
    """model_delta[g] = T(flag WT) - T(flag null), same RPE1 expression (fold-averaged);
    observed_delta[g] = Brown 2019 WT consensus - Brown 2019 null consensus. Pass iff all four < 0."""
    md = {g: model_delta.get(g, np.nan) for g in genes}
    od = {g: observed_delta.get(g, np.nan) for g in genes}
    if not all(np.isfinite(list(md.values()) + list(od.values()))):
        return {"result": INCONCLUSIVE, "reason": "missing gene", "model_delta": md, "observed_delta": od}
    model_ok = all(v < 0 for v in md.values())
    obs_ok = all(v < 0 for v in od.values())
    agree = int(sum(np.sign(md[g]) == np.sign(od[g]) for g in genes))
    res = "pass" if (model_ok and obs_ok) else "fail"
    return {"result": res, "model_condition": model_ok, "observed_condition": obs_ok,
            "sign_agreement": f"{agree}/{len(genes)}", "model_delta": md, "observed_delta": od,
            "note": "feature sanity check (single two-point counterfactual), not independent transfer evidence"}


def h2(model_delta: pd.Series, observed_deltas: dict) -> dict:
    """model_delta: per gene (p53-pathway set) T(flag WT) - T(flag null); observed_deltas: {study: Series
    WT - null consensus delta} (Brown 2019, Bowden 2020). Observed delta per gene = mean over the studies
    observing it (study = unit). Pass iff >= 10 genes, Spearman > 0.3 and sign agreement >= 70%."""
    from scipy.stats import spearmanr
    od = pd.concat(observed_deltas, axis=1).mean(axis=1, skipna=True) if observed_deltas else pd.Series(dtype=float)
    g = od.dropna().index.intersection(model_delta.dropna().index)
    n = len(g)
    per_study = {s: int(d.dropna().index.intersection(model_delta.dropna().index).size) for s, d in observed_deltas.items()}
    if n < H2_MIN_GENES:
        return {"result": INCONCLUSIVE, "n_genes": n, "reason": f"< {H2_MIN_GENES} observed genes", "n_per_study": per_study}
    m, o = model_delta[g].to_numpy(float), od[g].to_numpy(float)
    rho = float(spearmanr(m, o).statistic) if (np.ptp(m) > 0 and np.ptp(o) > 0) else 0.0
    sign = float(np.mean(np.sign(m) == np.sign(o)))
    ok = rho > H2_MIN_SPEARMAN and sign >= H2_MIN_SIGN
    return {"result": "pass" if ok else "fail", "n_genes": n, "spearman": rho, "sign_agreement": sign,
            "n_per_study": per_study, "genes": list(g)}


def h3(observed_mair: pd.Series, pred_T: pd.Series, pred_noEC: pd.Series, genes) -> dict:
    """Among the pluripotency genes observed < -0.5 in the Mair 2019 consensus (>= 3 needed), pass if
    model T predicts < -0.5 for >= 50%; the no-EC variant's count is reported; 'narrow EC transfer' if
    only the EC-trained model passes."""
    obs = observed_mair.reindex(list(genes))
    dep = [g for g in genes if np.isfinite(obs.get(g, np.nan)) and obs[g] < EFFECT_THR
           and np.isfinite(pred_T.get(g, np.nan))]
    if len(dep) < H3_MIN_OBSERVED:
        return {"result": INCONCLUSIVE, "n_observed_depleted": len(dep), "observed_depleted": dep,
                "reason": f"< {H3_MIN_OBSERVED} genes observed < -0.5 in Mair 2019"}
    k_t = int(sum(pred_T[g] < EFFECT_THR for g in dep))
    k_n = int(sum(pred_noEC.get(g, np.nan) < EFFECT_THR for g in dep))
    pass_t = k_t >= H3_PASS_FRACTION * len(dep)
    pass_n = k_n >= H3_PASS_FRACTION * len(dep)
    res = ("pass" if pass_n else "narrow EC transfer") if pass_t else "fail"
    return {"result": res, "n_observed_depleted": len(dep), "observed_depleted": dep,
            "n_predicted_T": k_t, "n_predicted_noEC": k_n, "T_passes": bool(pass_t), "noEC_passes": bool(pass_n)}


def h4() -> dict:
    return {"result": INCONCLUSIVE, "reason": "declared: primary T-cell screens fail QC (Shifrut 2018 SSMD -0.81)"}


# ================================================================ populations from the criteria file
def build_entries(bg_spec: dict, load_arm) -> tuple[list, list]:
    """Entries of one background from ``transfer_criteria.json -> backgrounds.<bg>.studies``.

    bg_spec["studies"] = {study: {"role": "verdict" | "predictor_only", "entries": {tp53: {"arms": [...]}}}}.
    load_arm(arm_id) -> Series (harmonised 'final' profile). Returns (verdict entries, all entries);
    entry name = "<study>:<tp53>", observed = consensus of its arms.
    """
    eligible, everyone = [], []
    for study, s in sorted(bg_spec["studies"].items()):
        for tp, e in sorted(s["entries"].items()):
            ent = Entry(name=f"{study}:{tp}", study=study, tp53=tp,
                        observed=consensus([load_arm(a) for a in e["arms"]]))
            everyone.append(ent)
            if s["role"] == "verdict":
                eligible.append(ent)
    return eligible, everyone


def loso_pools(bg_spec: dict) -> dict:
    """{verdict study: sorted list of the OTHER studies of the background (any role)}."""
    studies = sorted(bg_spec["studies"])
    return {s: [o for o in studies if o != s] for s in studies if bg_spec["studies"][s]["role"] == "verdict"}
