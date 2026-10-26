"""Known-answer tests for essml.transfer (D13 verdict computations; toy data only, no real screens)."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import average_precision_score

from essml import metrics as M
from essml import transfer as TR

ROOT = Path(__file__).resolve().parents[1]
G = 3000
GENES = pd.Index([f"{i}" for i in range(G)])
EXCL = set(GENES[:100])            # toy "scaling halves"
NEG = list(GENES[100:400])          # toy NEGv1 evaluation half


def toy_world(seed=0, n_arms=(1, 2, 1, 3), noise=0.15, study_noise=0.05):
    """B0 = population mean; background truth = B0 + sparse deviations; arms = truth + noise."""
    rng = np.random.default_rng(seed)
    m = rng.normal(-0.1, 0.25, G)
    m[400:700] = rng.normal(-1.5, 0.3, 300)        # common essentials (high dependency fraction)
    m[100:400] = rng.normal(0.0, 0.05, 300)        # non-essentials
    d = np.zeros(G)
    spec = rng.choice(np.arange(700, G), 250, replace=False)
    d[spec] = rng.normal(-1.2, 0.3, 250)           # background-specific dependencies
    truth = pd.Series(m + d, index=GENES)
    b0 = pd.Series(m, index=GENES)
    frac = pd.Series(np.where(m < -1.0, 0.9, 0.05), index=GENES)
    entries = []
    for k, na in enumerate(n_arms):
        shift = rng.normal(0, study_noise, G)
        arms = [truth + shift + rng.normal(0, noise, G) for _ in range(na)]
        entries.append(TR.Entry(name=f"s{k}:WT", study=f"s{k}", observed=TR.consensus(arms), tp53="WT"))
    return truth, b0, frac, entries


def evals_for(entries, q, frac, plan, label="rank10"):
    return {e.name: TR.EntryEval(e, q, GENES, EXCL, frac_below=frac, neg_eval=NEG, label=label, plan_genes=plan.genes)
            for e in entries}


def plan_for(entries, n_boot=200, seed=20261026, pairs=None):
    pairs = pd.DataFrame(columns=["a", "b", "identity", "rev_identity"]) if pairs is None else pairs
    U = sorted(set().union(*[set(TR.gene_universe(e.observed, GENES, EXCL)) for e in entries]))
    return TR.BootstrapPlan(sorted({e.study for e in entries}), U, TR.paralog_blocks(U, pairs), n_boot, seed)


# ---------------------------------------------------------------- metric equivalence with the frozen code
def test_weighted_ap_and_recall_match_frozen_metrics():
    rng = np.random.default_rng(3)
    for _ in range(100):
        n = int(rng.integers(20, 200))
        s = np.round(rng.normal(size=n), 1)            # ties
        y = rng.random(n) < 0.3
        y[0] = True
        w = rng.integers(0, 3, n).astype(float)
        w[0] = 1
        idx = np.repeat(np.arange(n), w.astype(int))
        assert TR.average_precision(s, y, w) == pytest.approx(average_precision_score(y[idx], s[idx]), abs=1e-12)
        assert TR.average_precision(s, y) == pytest.approx(average_precision_score(y, s), abs=1e-12)
        neg = (~y) & (rng.random(n) < 0.8)
        neg[1] = not y[1]
        r = TR.RecallPrep(s, y, neg)(w)
        ref = M.recall_at_fpr(s[idx], y[idx], neg[idx])
        assert (np.isnan(r) and np.isnan(ref)) or r == pytest.approx(ref, abs=1e-12)


def test_constant_score_ap_is_prevalence():
    y = np.r_[np.ones(10, bool), np.zeros(90, bool)]
    assert TR.average_precision(np.zeros(100), y) == pytest.approx(0.1)


# ---------------------------------------------------------------- universes and labels
def test_universe_excludes_scaling_halves_and_rank_label_is_ten_percent():
    _, _, _, entries = toy_world()
    U = TR.gene_universe(entries[0].observed, GENES, EXCL)
    assert not set(U) & EXCL and len(U) == G - 100
    assert set(NEG) <= set(U)                      # evaluation-half controls retained
    lab = TR.rank_label(entries[0].observed, U)
    assert lab.sum() == int(np.ceil(0.1 * len(U))) or abs(lab.mean() - 0.1) < 1e-3
    cs = {"ceg_scaling": ["1"], "neg_scaling": ["2"], "scaling_set": ["3", "1"], "ceg_evaluation": ["4"]}
    assert TR.scaling_halves(cs) == {"1", "2", "3"}


def test_q3b_positives_exclude_negatives_and_common_essentials():
    truth, b0, frac, entries = toy_world()
    U = TR.gene_universe(entries[0].observed, GENES, EXCL)
    pos, neg = TR.q3b_sets(entries[0].observed, U, frac, NEG)
    assert not np.any(pos & neg)
    assert np.all(frac.reindex(U).to_numpy()[pos] <= 0.5)
    assert neg.sum() == len(NEG)


# ---------------------------------------------------------------- verdicts end to end
def _skills(entries, preds, b0, frac, q="Q3a", seed=20261026, n_boot=200):
    plan = plan_for(entries, n_boot, seed)
    ev = evals_for(entries, q, frac, plan)
    return TR.skill_bootstrap(ev, {n: preds[n] for n in ev}, {n: b0 for n in ev}, plan, keep_boot=True), plan, ev


def _ceiling(entries, b0, frac, plan, q="Q3a"):
    ev = evals_for(entries, q, frac, plan)
    spec = {"studies": {e.study: {"role": "verdict", "entries": {}} for e in entries}}
    return TR.ceiling_bootstrap(ev, TR.loso_pools(spec), TR.study_profiles(entries), {n: b0 for n in ev}, b0, plan)


def test_perfect_predictor_transfers():
    truth, b0, frac, entries = toy_world()
    preds = {e.name: e.observed for e in entries}      # oracle: predicts each entry's own consensus
    sk, plan, _ = _skills(entries, preds, b0, frac)
    ce = _ceiling(entries, b0, frac, plan)
    assert sk["point"] > 0 and sk["ci_lo"] > 0
    assert ce["point"] > 0 and sk["point"] > ce["point"]   # oracle beats independent replicates
    v = TR.background_verdict(sk, ce["point"], ce["mde"], True, False)
    assert v["verdict"] == TR.TRANSFERS


def test_b0_equal_predictor_has_zero_skill_and_never_transfers():
    truth, b0, frac, entries = toy_world()
    sk, plan, _ = _skills(entries, {e.name: b0 for e in entries}, b0, frac)
    assert sk["point"] == 0.0 and sk["ci_lo"] == 0.0 and sk["ci_hi"] == 0.0
    ce = _ceiling(entries, b0, frac, plan)
    v = TR.background_verdict(sk, ce["point"], ce["mde"], harmonisation_ok=True, anchor_gate_positive=True)
    assert v["verdict"] == TR.FAILS                                  # 0 < MDE, anchors positive
    v2 = TR.background_verdict(sk, ce["point"], ce["mde"], harmonisation_ok=True, anchor_gate_positive=False)
    assert v2["verdict"] == TR.INCONCLUSIVE                          # fails needs the anchor gate
    v3 = TR.background_verdict(sk, ce["point"], ce["mde"], harmonisation_ok=False, anchor_gate_positive=True)
    assert v3["verdict"] == TR.INCONCLUSIVE                          # and harmonisation validation


def test_noise_predictor_is_inconclusive_without_anchor_evidence():
    truth, b0, frac, entries = toy_world()
    rng = np.random.default_rng(7)
    noise = {e.name: pd.Series(rng.normal(-0.3, 0.6, G), index=GENES) for e in entries}
    sk, plan, _ = _skills(entries, noise, b0, frac)
    assert sk["ci_lo"] < 0                                           # cannot transfer
    ce = _ceiling(entries, b0, frac, plan)
    v = TR.background_verdict(sk, ce["point"], ce["mde"], True, anchor_gate_positive=False)
    assert v["verdict"] == TR.INCONCLUSIVE


def test_ci_covering_zero_above_mde_is_inconclusive():
    sk = {"point": 0.02, "ci_lo": -0.05, "ci_hi": 0.09}
    v = TR.background_verdict(sk, ceiling_skill=0.03, mde=0.08, harmonisation_ok=True, anchor_gate_positive=True)
    assert v["verdict"] == TR.INCONCLUSIVE and not v["transfers_condition"] and not v["fails_condition"]


def test_ratio_condition_and_boundaries():
    sk = {"point": 0.014, "ci_lo": 0.001, "ci_hi": 0.03}
    assert TR.background_verdict(sk, 0.03, 0.1, True, False)["verdict"] == TR.INCONCLUSIVE   # 0.014 < 0.015
    sk["point"] = 0.015
    assert TR.background_verdict(sk, 0.03, 0.1, True, False)["verdict"] == TR.TRANSFERS     # >= is inclusive
    # upper CI == MDE is not < MDE (strict)
    z = {"point": 0.0, "ci_lo": -0.1, "ci_hi": 0.1}
    assert TR.background_verdict(z, 0.03, 0.1, True, True)["verdict"] == TR.INCONCLUSIVE
    # both conditions met -> inconclusive
    b = {"point": 0.02, "ci_lo": 0.01, "ci_hi": 0.03}
    out = TR.background_verdict(b, 0.03, 0.1, True, True)
    assert out["transfers_condition"] and out["fails_condition"] and out["verdict"] == TR.INCONCLUSIVE
    # non-positive ceiling: ratio condition vacuous, flagged
    out = TR.background_verdict({"point": 0.02, "ci_lo": 0.01, "ci_hi": 0.03}, -0.1, 0.005, True, False)
    assert out["verdict"] == TR.TRANSFERS and out["ratio_condition_vacuous"]


def test_expression_swap_null_caps_transfers():
    truth, b0, frac, entries = toy_world()
    preds = {e.name: e.observed for e in entries}
    sk, plan, ev = _skills(entries, preds, b0, frac)
    ce = _ceiling(entries, b0, frac, plan)
    # swapped profile that still carries skill (here: the same oracle) -> cap
    swap = TR.skill_bootstrap(ev, preds, {n: b0 for n in ev}, plan)
    v = TR.background_verdict(sk, ce["point"], ce["mde"], True, False, swap_results=[swap])
    assert v["swap_cap"] and v["verdict"] == TR.INCONCLUSIVE
    # swapped profile without skill (B0) -> no cap
    flat = TR.skill_bootstrap(ev, {n: b0 for n in ev}, {n: b0 for n in ev}, plan)
    assert TR.background_verdict(sk, ce["point"], ce["mde"], True, False, swap_results=[flat])["verdict"] == TR.TRANSFERS


def test_q3b_end_to_end_perfect_vs_comparator():
    truth, b0, frac, entries = toy_world()
    sk, plan, ev = _skills(entries, {e.name: e.observed for e in entries}, b0, frac, q="Q3b")
    assert all(e.valid for e in ev.values())
    assert sk["point"] > 0 and sk["ci_lo"] > 0


# ---------------------------------------------------------------- bootstrap reproducibility and pairing
def test_hierarchical_bootstrap_reproducible_with_seed():
    truth, b0, frac, entries = toy_world()
    rng = np.random.default_rng(1)
    preds = {e.name: b0 + pd.Series(rng.normal(0, 0.2, G), index=GENES) for e in entries}
    a, _, _ = _skills(entries, preds, b0, frac, seed=20261026)
    b, _, _ = _skills(entries, preds, b0, frac, seed=20261026)
    c, _, _ = _skills(entries, preds, b0, frac, seed=1)
    np.testing.assert_array_equal(a["boot"], b["boot"])
    assert not np.array_equal(a["boot"], c["boot"])
    assert (a["ci_lo"], a["ci_hi"]) == (b["ci_lo"], b["ci_hi"])


def test_skill_and_ceiling_share_study_and_gene_draws():
    _, _, _, entries = toy_world()
    plan = plan_for(entries, n_boot=5)
    draws1 = [(si.copy(), w.copy()) for _, si, w in plan.replicates()]
    draws2 = [(si.copy(), w.copy()) for _, si, w, _ in plan.replicates(with_pool=True)]
    for (s1, w1), (s2, w2) in zip(draws1, draws2):
        np.testing.assert_array_equal(s1, s2)
        np.testing.assert_array_equal(w1, w2)


def test_gene_blocks_resampled_together():
    _, _, _, entries = toy_world()
    pairs = pd.DataFrame({"a": ["500", "501", "900"], "b": ["501", "502", "901"],
                          "identity": [50, 50, 10], "rev_identity": [60, 40, 10]})
    blocks = TR.paralog_blocks(GENES, pairs)
    assert blocks["500"] == blocks["501"] == blocks["502"]
    assert blocks["900"] != blocks["901"]                            # below 30% identity
    plan = plan_for(entries, n_boot=50, pairs=pairs)
    i = plan.genes.get_indexer(["500", "501", "502"])
    for _, _, w in plan.replicates():
        assert w[i[0]] == w[i[1]] == w[i[2]]


# ---------------------------------------------------------------- study is the unit
def test_study_is_the_unit():
    truth, b0, frac, entries = toy_world(n_arms=(1, 1, 1, 1))
    # (a) replicated arms do not add weight: a study made of 5 copies of one arm == that arm
    e0 = entries[0]
    dup = TR.Entry(name=e0.name, study=e0.study, observed=TR.consensus([e0.observed] * 5))
    pd.testing.assert_series_equal(dup.observed, e0.observed.astype(float), check_names=False)
    # (b) a study with two entries counts once: background = mean over studies of entry means
    extra = TR.Entry(name="s0:null", study="s0", observed=entries[0].observed * 0.5)
    ents = entries + [extra]
    preds = {e.name: e.observed for e in ents}
    sk, plan, ev = _skills(ents, preds, b0, frac)
    pe = sk["per_entry"]
    s0 = np.mean([pe["s0:WT"], pe["s0:null"]])
    expect = np.mean([s0, pe["s1:WT"], pe["s2:WT"], pe["s3:WT"]])
    assert sk["point"] == pytest.approx(expect)
    assert plan.studies == ["s0", "s1", "s2", "s3"]                 # 4 resampling units, not 5 entries
    # (c) study draws are with replacement over studies
    sis = np.array([si for _, si, _ in plan.replicates()])
    assert sis.shape == (plan.n_boot, 4) and sis.max() < 4
    assert any(len(set(r)) < 4 for r in sis)


def test_between_study_variance_enters_the_ci():
    # same predictor, same genes: when half of the studies are unrelated to it, the study level of the
    # bootstrap must widen the CI far beyond the gene-sampling spread of four concordant studies
    truth, b0, frac, _ = toy_world()
    rng = np.random.default_rng(11)
    same = [TR.Entry(f"s{k}:WT", f"s{k}", truth + pd.Series(rng.normal(0, 0.05, G), index=GENES)) for k in range(4)]
    mixed = same[:2] + [TR.Entry(f"s{k}:WT", f"s{k}", pd.Series(rng.normal(-0.3, 0.6, G), index=GENES))
                        for k in (2, 3)]
    a, _, _ = _skills(same, {e.name: truth for e in same}, b0, frac)
    b, _, _ = _skills(mixed, {e.name: truth for e in mixed}, b0, frac)
    assert b["se"] > 2 * a["se"]
    assert b["ci_hi"] - b["ci_lo"] > 0.1


def test_loso_pools_exclude_own_study_and_include_predictor_only():
    spec = {"studies": {"mair": {"role": "verdict"}, "ihry": {"role": "predictor_only"},
                        "yilmaz": {"role": "predictor_only"}}}
    assert TR.loso_pools(spec) == {"mair": ["ihry", "yilmaz"]}


def test_single_verdict_study_ceiling_uses_predictor_only_studies():
    truth, b0, frac, entries = toy_world(n_arms=(2, 1, 1))
    target, others = entries[0], entries[1:]
    plan = plan_for([target], n_boot=100)
    ev = evals_for([target], "Q3a", frac, plan)
    pools = {target.study: [e.study for e in others]}
    ce = TR.ceiling_bootstrap(ev, pools, TR.study_profiles(entries), {target.name: b0}, b0, plan)
    assert ce["point"] > 0 and np.isfinite(ce["mde"]) and ce["mde"] == pytest.approx(2 * ce["se"])


def test_d7_fold_skills_are_averaged():
    truth, b0, frac, entries = toy_world()
    rng = np.random.default_rng(5)
    p1 = {e.name: e.observed for e in entries}
    p2 = {e.name: b0 + pd.Series(rng.normal(0, 0.3, G), index=GENES) for e in entries}
    s1, _, _ = _skills(entries, p1, b0, frac)
    s2, _, _ = _skills(entries, p2, b0, frac)
    s12, _, _ = _skills(entries, {n: {0: p1[n], 1: p2[n]} for n in p1}, b0, frac)
    assert s12["point"] == pytest.approx((s1["point"] + s2["point"]) / 2)


def test_nan_prediction_on_universe_raises():
    truth, b0, frac, entries = toy_world()
    bad = b0.copy()
    bad.iloc[500] = np.nan
    with pytest.raises(ValueError):
        _skills(entries, {e.name: bad for e in entries}, b0, frac)


# ---------------------------------------------------------------- anchor gate, overall, set A
def test_anchor_gate_and_overall():
    assert TR.anchor_gate({"HAP1": {"ci_lo": -0.01}, "A375": {"ci_lo": 0.02}})["informative_positive"]
    assert not TR.anchor_gate({"HAP1": {"ci_lo": -0.01}})["informative_positive"]
    assert not TR.anchor_gate({})["informative_positive"]
    ov = TR.overall_verdict
    assert ov({"RPE1": "transfers", "hPSC": "transfers"})["verdict"] == "transfers"
    assert ov({"RPE1": "fails", "hPSC": "fails"})["verdict"] == "fails"
    assert ov({"RPE1": "transfers", "hPSC": "fails"})["verdict"] == "mixed"
    assert ov({"RPE1": "transfers", "hPSC": "inconclusive"})["verdict"] == "inconclusive"


def test_set_a_descriptive():
    res = {"l1": {"point": 0.05, "ci_lo": 0.01, "ci_hi": 0.09}, "l2": {"point": -0.2, "ci_lo": -0.3, "ci_hi": -0.1}}
    dist = np.linspace(-0.1, 0.1, 101)
    out = TR.set_a_descriptive(res, dist)
    assert out["per_line"]["l1"]["category"] == "skill>0" and out["per_line"]["l2"]["category"] == "skill<0"
    assert out["per_line"]["l1"]["test_percentile"] == pytest.approx(100 * 76 / 101)
    assert out["per_line"]["l2"]["test_percentile"] == 0.0


# ---------------------------------------------------------------- hypotheses
def test_h1():
    assert TR.h1({"MDM2": -0.3, "PPM1D": -0.1}, {"MDM2": -0.8, "PPM1D": -0.2})["result"] == "pass"
    r = TR.h1({"MDM2": -0.3, "PPM1D": 0.1}, {"MDM2": -0.8, "PPM1D": -0.2})
    assert r["result"] == "fail" and r["sign_agreement"] == "1/2"
    assert TR.h1({"MDM2": -0.3}, {"MDM2": -0.8, "PPM1D": -0.2})["result"] == "inconclusive"


def test_h2():
    g = [str(i) for i in range(15)]
    obs = pd.Series(np.linspace(-1, 1, 15), index=g)
    assert TR.h2(obs * 0.3, {"brown": obs})["result"] == "pass"
    assert TR.h2(-obs, {"brown": obs})["result"] == "fail"
    assert TR.h2(obs, {"brown": obs.iloc[:9]})["result"] == "inconclusive"     # < 10 genes
    # study = unit: the observed delta is the mean over studies observing the gene
    r = TR.h2(obs, {"brown": obs, "bowden": obs.iloc[:5] * 3})
    assert r["n_genes"] == 15 and r["result"] == "pass"


def test_h3_narrow_ec():
    genes = ["POU5F1", "LIN28A", "SOX2", "NANOG"]
    obs = pd.Series({"POU5F1": -1.5, "LIN28A": -1.0, "SOX2": -0.8, "NANOG": 0.1})
    t = pd.Series({"POU5F1": -0.9, "LIN28A": -0.6, "SOX2": -0.1, "NANOG": -0.7})
    noec = pd.Series({"POU5F1": -0.6, "LIN28A": -0.1, "SOX2": -0.1, "NANOG": -0.7})
    r = TR.h3(obs, t, noec, genes)
    assert r["n_observed_depleted"] == 3 and r["n_predicted_T"] == 2 and r["result"] == "narrow EC transfer"
    assert TR.h3(obs, t, t, genes)["result"] == "pass"
    assert TR.h3(obs.clip(lower=-0.6) * 0 + 0.0, t, t, genes)["result"] == "inconclusive"
    assert TR.h4()["result"] == "inconclusive"


# ---------------------------------------------------------------- frozen file consistency
def test_criteria_file_matches_module_constants():
    p = ROOT / "configs" / "transfer_criteria.json"
    if not p.exists():
        pytest.skip("criteria file not written")
    c = json.loads(p.read_text(encoding="utf-8"))
    seeds = json.loads((ROOT / "configs" / "seeds.json").read_text())
    assert c["bootstrap"]["n_resamples"] == TR.N_BOOT == 2000
    assert c["bootstrap"]["seed"] == seeds["split_seed"]
    for bg in ("RPE1", "hPSC"):
        b = c["backgrounds"][bg]
        for k in ("Q3a_primary", "Q3b_primary"):
            r = b["replicate_ceiling"][k]
            assert r["verdict_bearing"] and r["MDE"] == pytest.approx(2 * r["bootstrap_se"], abs=2e-4)
        assert b["nearest_depmap_lineage"]["lineage"]
    assert c["backgrounds"]["hPSC"]["studies"]["ihry2019"]["role"] == "predictor_only"
    assert not any(s.startswith("bowden") for s in c["backgrounds"]["RPE1"]["studies"])
