"""Toy tests for the mutation-driven dependency diagnostic (essml.mutdeps)."""
import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import average_precision_score, roc_auc_score

from essml import mutdeps as MD

IDS = [f"L{i}" for i in range(6)]


def _mats():
    hot = pd.DataFrame({"G": [1, 0, 0, 0, 0, 0], "TP53": [0, 1, 0, 0, 0, 0]}, index=IDS, dtype=float)
    dam = pd.DataFrame({"G": [0, 0, 2, 0, 0, 0], "TP53": [0, 0, 1, 0, 0, 0]}, index=IDS, dtype=float)
    # L5 has no CN profile; L4 has NaN CN for G (diploid)
    cn = pd.DataFrame({"G": [1.0, 2.5, 1.0, 0.1, np.nan]}, index=IDS[:5])
    return hot, dam, cn


def test_status_components_and_unknown():
    hot, dam, cn = _mats()
    s = lambda t, g="G": MD.biomarker_status(t, g, IDS, hot, dam, cn).tolist()
    assert s("hotspot") == [1, 0, 0, 0, 0, 0]
    assert s("damaging") == [0, 0, 1, 0, 0, 0]
    amp = s("amp")
    assert amp[:5] == [0, 1, 0, 0, 0] and np.isnan(amp[5])          # no CN profile -> unknown
    dd = s("deep_del")
    assert dd[:5] == [0, 0, 0, 1, 0] and np.isnan(dd[5])           # NaN CN -> diploid
    hd = s("hotspot|amp")
    assert hd[:5] == [1, 1, 0, 0, 0] and np.isnan(hd[5])
    # positive via a known component even if another component is unknown
    hot2 = hot.copy()
    hot2.loc["L5", "G"] = 1
    assert MD.biomarker_status("hotspot|amp", "G", IDS, hot2, dam, cn).tolist()[5] == 1
    assert s("TP53-WT", "TP53") == [1, 0, 0, 1, 1, 1]                # mutant = wild type
    with pytest.raises(ValueError):
        s("bogus")


def test_fold_metrics_strata_match_direct_computation():
    rng = np.random.default_rng(1)
    n = 60
    status = np.r_[np.ones(20), np.zeros(35), np.full(5, np.nan)]
    true = rng.normal(size=n) - 0.8 * (status == 1)
    pred = true + rng.normal(scale=0.7, size=n)
    true[3] = np.nan                                                  # unobserved entry ignored
    out = MD.fold_metrics(pred, true, status, min_obs_all=10)
    obs = ~np.isnan(true)
    for name, m in (("mut", status == 1), ("wt", status == 0)):
        mm = m & obs
        assert out[f"n_{name}"] == mm.sum()
        assert out[f"r_{name}"] == pytest.approx(np.corrcoef(pred[mm], true[mm])[0, 1])
        y = true[mm] < -0.5
        assert out[f"ap_{name}"] == pytest.approx(average_precision_score(y, -pred[mm]))
        pi = y.mean()
        assert out[f"nap_{name}"] == pytest.approx((out[f"ap_{name}"] - pi) / (1 - pi))
    # all-line metrics include the unknown-status lines
    assert out["n_obs"] == obs.sum() == 59
    assert out["r_all"] == pytest.approx(np.corrcoef(pred[obs], true[obs])[0, 1])
    # biomarker AUROC: dependent known-status lines only
    dep = obs & (true < -0.5) & ~np.isnan(status)
    assert out["auroc_bm"] == pytest.approx(roc_auc_score(status[dep] == 1, -pred[dep]))


def test_small_strata_are_undefined():
    status = np.r_[np.ones(3), np.zeros(50)]
    true = np.r_[np.full(3, -1.0), np.linspace(-1, 1, 50)]
    pred = true.copy()
    out = MD.fold_metrics(pred, true, status, min_obs_all=10, min_obs_stratum=5)
    assert np.isnan(out["r_mut"]) and np.isnan(out["ap_mut"])       # 3 < 5 lines; no negatives
    assert np.isfinite(out["r_wt"]) and np.isfinite(out["ap_all"])
    assert out["auroc_bm"] == out["auroc_bm"]                        # 3 mutant dependent -> defined


def test_biomarker_only_uses_training_means_only():
    t_tr = np.array([-1.0, -2.0, 0.0, 1.0, np.nan, 5.0])
    s_tr = np.array([1, 1, 0, 0, 1, np.nan])
    s_te = np.array([1, 0, np.nan])
    p = MD.biomarker_only(t_tr, s_tr, s_te)
    # mutant mean -1.5 (NaN target ignored), WT 0.5, unknown -> overall training mean (0.6)
    assert p.tolist() == pytest.approx([-1.5, 0.5, 0.6])
    # empty training stratum -> training gene mean
    p2 = MD.biomarker_only(np.array([0.0, 2.0]), np.array([0, 0]), np.array([1, 0]))
    assert p2.tolist() == pytest.approx([1.0, 1.0])
    # within a stratum the reference is constant -> r = 0 by the frozen convention
    out = MD.fold_metrics(np.r_[np.full(10, -1.5), np.full(10, 0.5)],
                          np.r_[np.linspace(-2, 0, 10), np.linspace(-1, 1, 10)],
                          np.r_[np.ones(10), np.zeros(10)], min_obs_all=5)
    assert out["r_mut"] == 0.0 and out["r_wt"] == 0.0


def test_clearly_helps_rule():
    big = {f"p{i}": np.array([0.2, 0.15, 0.12, 0.11, 0.3]) for i in range(3)}
    res = MD.clearly_helps({**big, "q": np.array([0.0, 0.01, -0.01, 0.0, 0.02])},
                           eligible=["p0", "p1", "p2", "q"])
    assert res["met"] and res["cond1_ci_excludes_0"] and res["pairs_gain"] == ["p0", "p1", "p2"]
    # 3 of 5 folds >= 0.10 is not enough; a NaN fold counts as no gain
    two = {"p0": np.array([0.2, 0.2, 0.2, np.nan, 0.0]), "p1": big["p1"], "p2": big["p2"]}
    r2 = MD.clearly_helps(two, eligible=list(two))
    assert r2["n_folds_gain_ge"]["p0"] == 3 and not r2["cond2_ge3_pairs"]
    assert r2["pooled_pairs"] == ["p1", "p2"]                         # p0 undefined in a fold
    # ineligible pairs are ignored; a CI that straddles 0 fails condition 1
    noisy = {"a": np.array([0.3, -0.3, 0.2, -0.2, 0.05])}
    r3 = MD.clearly_helps({**noisy, **big}, eligible=["a"])
    assert not r3["cond1_ci_excludes_0"] and not r3["met"] and r3["pooled_pairs"] == ["a"]


def test_paired_ci_matches_t_interval():
    from scipy import stats
    d = np.array([0.01, 0.02, 0.0, 0.03, 0.015])
    ci = MD.paired_ci(d)
    se = d.std(ddof=1) / np.sqrt(5)
    h = stats.t.ppf(0.975, 4) * se
    assert ci["ci95"] == pytest.approx([d.mean() - h, d.mean() + h])
