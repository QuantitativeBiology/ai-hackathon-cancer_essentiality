"""Per-gene predictability tiers, confidence and the null-calibration escalation (DESIGN §4, D22); toy data."""
import numpy as np
import pandas as pd
import pytest

from essml import predictability as PR


def stats_from(rows):
    """rows: {gene: [r per fold]} -> fold_stats output."""
    r = pd.DataFrame(rows).T
    r.columns = list(range(r.shape[1]))
    return PR.fold_stats(r)


def test_tier_boundaries_inclusive_and_strict():
    s = stats_from({
        "hi_exact": [0.30, 0.30, 0.30, 0.20, 0.40],     # mean 0.30 (>=), min 0.20 (>=) -> high
        "hi_min_below": [0.40, 0.40, 0.40, 0.19, 0.11],  # mean 0.30, min 0.11 < 0.20 -> medium (min > 0)
        "med_exact": [0.15, 0.15, 0.15, 0.15, 0.15],     # mean 0.15 (>=), min 0.15 > 0 -> medium
        "med_min_zero": [0.30, 0.30, 0.15, 0.0, 0.0],    # mean 0.15, min 0 (not > 0) -> low
        "low_mean": [0.149, 0.149, 0.149, 0.149, 0.149],
        "negative": [-0.2, -0.1, 0.0, 0.1, -0.3],
    })
    t = PR.assign_tiers(s)
    assert t["hi_exact"] == "high"
    assert t["hi_min_below"] == "medium"
    assert t["med_exact"] == "medium"
    assert t["med_min_zero"] == "low"
    assert t["low_mean"] == "low"
    assert t["negative"] == "low"


def test_not_evaluable_if_any_fold_nan():
    s = stats_from({"g_nan": [0.9, 0.9, np.nan, 0.9, 0.9], "g_ok": [0.9] * 5,
                    "g_all_nan": [np.nan] * 5})
    t = PR.assign_tiers(s)
    assert t["g_nan"] == "not_evaluable" and t["g_all_nan"] == "not_evaluable" and t["g_ok"] == "high"
    c = PR.confidence(s)
    assert np.isnan(c["g_nan"]) and np.isnan(c["g_all_nan"])
    assert s.loc["g_nan", "n_folds_evaluable"] == 4


def test_fold_stats_and_confidence_values():
    s = stats_from({"g": [0.1, 0.2, 0.3, 0.4, 0.5], "neg": [-0.9, -0.95, -1.0, -0.9, -1.0]})
    assert s.loc["g", "mean_r"] == pytest.approx(0.3)
    assert s.loc["g", "min_r"] == pytest.approx(0.1)
    assert s.loc["g", "sd_r"] == pytest.approx(np.std([0.1, 0.2, 0.3, 0.4, 0.5], ddof=1))
    c = PR.confidence(s)
    assert c["g"] == pytest.approx(0.3 - np.std([0.1, 0.2, 0.3, 0.4, 0.5], ddof=1))
    assert c["neg"] == -1.0  # clipped


def test_exclude_folds_uses_only_listed_folds():
    r = pd.DataFrame({"g": [np.nan, 0.5, 0.5, 0.5, 0.5]}).T
    r.columns = range(5)
    assert PR.assign_tiers(PR.fold_stats(r))["g"] == "not_evaluable"
    assert PR.assign_tiers(PR.fold_stats(r, folds=[1, 2, 3, 4]))["g"] == "high"


def _n_genes(n, mean, spread=0.0, prefix="g"):
    return {f"{prefix}{i}": [mean + spread] * 5 for i in range(n)}


def test_escalation_not_triggered_at_exactly_5_percent():
    obs = stats_from(_n_genes(100, 0.6))
    null = stats_from(_n_genes(5, 0.32, prefix="n") | _n_genes(50, 0.0, prefix="z"))
    res = PR.null_calibrate(obs, null)
    assert res["status"] == "not_triggered" and res["threshold"] == 0.30
    assert res["steps"][0]["n_high_null"] == 5 and res["steps"][0]["n_high_obs"] == 100


def test_escalation_picks_smallest_passing_ladder_value():
    obs = stats_from(_n_genes(100, 0.6))
    # 6 null genes at 0.32 (fail at 0.30), 6 at 0.37 (fail at 0.35 too), none above 0.40 -> 0.40
    null = stats_from(_n_genes(6, 0.32, prefix="a") | _n_genes(6, 0.37, prefix="b"))
    res = PR.null_calibrate(obs, null)
    assert res["status"] == "escalated" and res["threshold"] == 0.40
    assert [s["pass"] for s in res["steps"]] == [False, False, True, True]


def test_escalation_ladder_exhausted():
    obs = stats_from(_n_genes(10, 0.9))
    null = stats_from(_n_genes(5, 0.9, prefix="n"))
    res = PR.null_calibrate(obs, null)
    assert res["status"] == "ladder_exhausted" and res["threshold"] == 0.50


def test_escalation_counts_observed_high_at_each_threshold():
    # observed highs sit at 0.36: escalating to 0.40 would remove them all -> 0.35 must be chosen
    obs = stats_from(_n_genes(100, 0.36))
    null = stats_from(_n_genes(10, 0.31, prefix="n"))
    res = PR.null_calibrate(obs, null)
    assert res["threshold"] == 0.35 and res["steps"][1]["n_high_obs"] == 100


def test_build_table_covers_all_genes_and_records_threshold():
    r = pd.DataFrame({"1": [0.5] * 5, "2": [0.2] * 5, "3": [np.nan] + [0.5] * 4}).T
    r.columns = range(5)
    sym = pd.Series({"1": "A", "2": "B", "3": "C", "4": "D"})
    t = PR.build_table(r, range(5), 0.35, sym, "run_x", all_genes=["1", "2", "3", "4"])
    assert list(t["entrez"]) == ["1", "2", "3", "4"]
    assert list(t["tier"]) == ["high", "medium", "not_evaluable", "not_evaluable"]
    assert (t["threshold_high_mean_r"] == 0.35).all() and (t["source_run_id"] == "run_x").all()
    assert {"r_fold0", "r_fold4", "mean_r", "min_r", "sd_r", "confidence"} <= set(t.columns)


def test_fold_r_table_matches_frozen_min_obs_rule():
    rng = np.random.default_rng(0)
    ids = [f"L{i}" for i in range(100)]
    folds = pd.Series(np.arange(100) % 2, index=ids)
    T = pd.DataFrame(rng.normal(size=(100, 3)), index=ids, columns=["a", "b", "c"])
    T.loc[folds.index[folds == 1][:20], "c"] = np.nan          # fold 1: 30 observed lines of c
    P = T.fillna(0) + rng.normal(0, 0.5, size=T.shape)
    gc = pd.DataFrame({"all": True, "variable": True}, index=["a", "b", "c"])
    long = PR.fold_r_table(P, T, folds, gc, min_obs=40)
    w = PR.wide_r(long)
    assert np.isfinite(w.loc["c", 0]) and np.isnan(w.loc["c", 1])
    assert PR.assign_tiers(PR.fold_stats(w))["c"] == "not_evaluable"


def test_annotate_long_carries_confidence_and_tier(tmp_path):
    tab = pd.DataFrame({"entrez": ["1", "2"], "symbol": ["A", "B"], "mean_r": [0.5, 0.1],
                        "confidence": [0.4, 0.0], "tier": ["high", "low"], "source_run_id": "x",
                        "threshold_high_mean_r": 0.3})
    p = tmp_path / "pg.csv"
    tab.to_csv(p, index=False)
    t = PR.load_table(p)
    pred = pd.DataFrame([[1.0, 2.0, 3.0]], index=pd.Index(["L0"], name="ModelID"), columns=["1", "2", "9"])
    long = PR.annotate_long(pred, t)
    assert list(long.columns) == ["ModelID", "entrez", "symbol", "predicted_effect", "confidence", "tier"]
    assert long.set_index("entrez").loc["9", "tier"] == "not_evaluable"
    assert long.set_index("entrez").loc["1", "confidence"] == 0.4
