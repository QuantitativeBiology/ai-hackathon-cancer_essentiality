"""Toy-data tests with known answers for essml.harmonise (DESIGN v3 §7.2, §8)."""
import numpy as np
import pandas as pd
import pytest

from essml import harmonise as H


def _genes(n, prefix="g"):
    return [f"{prefix}{i}" for i in range(n)]


# ---------------------------------------------------------------- (i) counts -> LFC -> gene score
def test_guide_lfc_rpm_and_filter():
    t0 = pd.DataFrame({"T0": [1000, 1000, 10, 2000]}, index=list("abcd"))
    end = pd.DataFrame({"E1": [500, 2000, 10, 2000], "E2": [500, 2000, 10, 2000]}, index=list("abcd"))
    lfc, rep = H.guide_lfc(t0, end, pseudocount=0.0, min_t0=30)
    assert "c" not in lfc.index and rep["n_guides_kept"] == 3  # guide c below 30 T0 reads
    # equal library sizes up to the dropped guide: totals 4010 vs 4510
    exp_a = np.log2(500 / 4510) - np.log2(1000 / 4010)
    assert lfc["a"] == pytest.approx(exp_a)
    assert lfc["b"] - lfc["a"] == pytest.approx(2.0)
    assert rep["t0_filter"] == "mean T0 reads >= 30"


def test_rpm_within_pool():
    c = pd.DataFrame({"s": [10, 30, 100, 300]}, index=list("abcd"))
    r = H.rpm(c, groups=["p1", "p1", "p2", "p2"])
    assert r["s"].tolist() == pytest.approx([250000, 750000, 250000, 750000])


def test_gene_score_median_and_min_guides():
    v = pd.Series([1.0, 2.0, 9.0, -1.0, 5.0], index=list("abcde"))
    g = pd.Series(["X", "X", "X", "Y", np.nan], index=list("abcde"))
    s = H.gene_scores(v, g, min_guides=2)
    assert s["X"] == 2.0 and "Y" not in s.index  # median; Y has a single guide


# ---------------------------------------------------------------- control split
def test_control_split_fixed_seed_and_scaling_set():
    ceg = {str(i) for i in range(100)}
    neg = {str(i) for i in range(1000, 1080)}
    ce = {str(i) for i in range(50, 300)}
    a = H.control_split(ceg, neg, ce, seed=1)
    b = H.control_split(ceg, neg, ce, seed=1)
    assert a == b
    assert set(a["ceg_scaling"]) | set(a["ceg_evaluation"]) == ceg
    assert not set(a["ceg_scaling"]) & set(a["ceg_evaluation"])
    assert len(a["neg_scaling"]) == 40
    assert not set(a["scaling_set"]) & set(a["ceg_evaluation"])
    assert set(a["scaling_set"]) == ce - set(a["ceg_evaluation"])


# ---------------------------------------------------------------- (ii) affine scaling maps planted medians exactly
def test_affine_maps_planted_medians_exactly():
    rng = np.random.default_rng(0)
    genes = _genes(2000)
    x = pd.Series(rng.normal(0, 1, 2000), index=genes)
    scal, neg = genes[:101], genes[200:301]  # odd sizes: medians are data points
    x[scal] = rng.normal(-3.7, 0.5, 101)
    x[neg] = rng.normal(0.4, 0.2, 101)
    ms, mn = H.affine_params(x, scal, neg)
    assert ms == pytest.approx(np.median(x[scal])) and mn == pytest.approx(np.median(x[neg]))
    y = H.affine_apply(x, ms, mn)
    assert np.median(y[scal]) == pytest.approx(-1.0, abs=1e-12)
    assert np.median(y[neg]) == pytest.approx(0.0, abs=1e-12)
    # an affine transform of the input gives the same output
    y2 = H.affine_apply(2.5 * x - 1.3, *H.affine_params(2.5 * x - 1.3, scal, neg))
    assert np.allclose(y, y2)


def test_qc_metrics_known_values_and_affine_invariance():
    e = pd.Series([-3.0, -2.0, -1.0], index=["e1", "e2", "e3"])
    n = pd.Series([-1.0, 0.0, 1.0], index=["n1", "n2", "n3"])
    s = pd.concat([e, n])
    q = H.qc_metrics(s, e.index, n.index)
    assert q["median_ceg_eval"] == -2 and q["median_neg_eval"] == 0
    assert q["nnmd"] == pytest.approx(-2 / 1.4826)
    assert q["ssmd"] == pytest.approx(-2 / np.sqrt(2))
    assert q["auroc"] == pytest.approx(1 - 0.5 / 9)  # one tie (-1 vs -1)
    q2 = H.qc_metrics(3 * s + 7, e.index, n.index)
    assert q2["nnmd"] == pytest.approx(q["nnmd"]) and q2["ssmd"] == pytest.approx(q["ssmd"])
    assert H.qc_pass(q) and not H.qc_pass({**q, "ssmd": -0.5})


# ---------------------------------------------------------------- (iii) EB correction
def _anchor_data(n_anchor=6, G=3000, delta_sd=0.4, noise=0.15, seed=0, p53=()):
    rng = np.random.default_rng(seed)
    genes = _genes(G)
    Y = pd.DataFrame(rng.normal(-0.2, 0.5, (n_anchor, G)), index=[f"a{i}" for i in range(n_anchor)], columns=genes)
    delta = pd.Series(rng.normal(0, delta_sd, G), index=genes)
    X = Y + delta + rng.normal(0, noise, (n_anchor, G))
    return X, Y, delta


def test_eb_recovers_planted_offsets():
    X, Y, delta = _anchor_data()
    c = H.fit_genewise(X, Y, kappa=1.0)
    # tau2 ~ var(delta) = 0.16; posterior weight B = tau2 / (tau2 + sigma2/A)
    assert c.tau2 == pytest.approx(0.16, rel=0.1)
    B = c.tau2 / (c.tau2 + 0.15 ** 2 / 6)
    assert np.corrcoef(c.delta, delta)[0, 1] > 0.98
    slope = np.polyfit(delta, c.delta, 1)[0]
    assert slope == pytest.approx(B, abs=0.03)
    # corrected new line with the same planted offset is closer to the truth
    rng = np.random.default_rng(5)
    y_new = pd.Series(rng.normal(0, 0.5, X.shape[1]), index=X.columns)
    x_new = y_new + delta + rng.normal(0, 0.15, X.shape[1])
    xc = H.apply_genewise(x_new, c)
    assert np.sqrt(((xc - y_new) ** 2).mean()) < 0.5 * np.sqrt(((x_new - y_new) ** 2).mean())


def test_eb_single_gene_offset_recovered_in_direction():
    """A large planted offset on one gene among many offset genes is recovered (shrunk)."""
    X, Y, delta = _anchor_data(delta_sd=0.3, seed=2)
    X = X.copy()
    X["g7"] = Y["g7"] + 1.5
    c = H.fit_genewise(X, Y, kappa=2.0)
    assert 1.0 < c.delta["g7"] < 1.5


def test_eb_strong_shrinkage_and_no_signal():
    X, Y, delta = _anchor_data(delta_sd=0.0, noise=0.3)
    c = H.fit_genewise(X, Y)
    assert np.abs(c.delta).max() < 0.1  # no planted offsets -> near-identity correction
    c1 = H.fit_genewise(*_anchor_data()[:2], kappa=1.0)
    c4 = H.fit_genewise(*_anchor_data()[:2], kappa=4.0)
    assert np.abs(c4.delta).mean() < np.abs(c1.delta).mean()  # larger kappa = stronger shrinkage


def test_single_anchor_not_estimable():
    X, Y, _ = _anchor_data(n_anchor=1)
    c = H.fit_genewise(X, Y)
    assert not c.estimable and (c.delta == 0).all()
    xc = H.apply_genewise(X.iloc[0], c)
    assert np.allclose(xc, X.iloc[0])


def test_scale_only_with_enough_anchors():
    rng = np.random.default_rng(3)
    X, Y, delta = _anchor_data(n_anchor=12, delta_sd=0.2, noise=0.05)
    lam = pd.Series(np.exp(rng.normal(0, 0.3, X.shape[1])), index=X.columns)
    c_g = Y.mean()
    Xs = c_g + (Y - c_g) * lam + delta + rng.normal(0, 0.05, X.shape)
    c = H.fit_genewise(Xs, Y, scale=True)
    assert c.scale and np.corrcoef(c.log_scale, np.log(lam))[0, 1] > 0.7
    c7 = H.fit_genewise(Xs.iloc[:7], Y.iloc[:7], scale=True)
    assert not c7.scale and (c7.log_scale == 0).all()


def test_p53_ddr_genes_untouched():
    X, Y, _ = _anchor_data()
    excl = ["g1", "g2", "g3"]
    X = X.copy()
    X[excl] = Y[excl] + 2.0  # huge offset that the fit must ignore
    c = H.fit_genewise(X, Y, exclude=excl, scale=False)
    assert (c.delta[excl] == 0).all() and set(c.excluded) == set(excl)
    x = X.iloc[0]
    xc = H.apply_genewise(x, c)
    assert (xc[excl] == x[excl]).all()
    assert not np.allclose(xc.drop(excl), x.drop(excl))
    X12, Y12, _ = _anchor_data(n_anchor=12)
    c12 = H.fit_genewise(X12, Y12, exclude=excl, scale=True)
    assert (H.apply_genewise(X12.iloc[0], c12)[excl] == X12.iloc[0][excl]).all()
    # the DESIGN §7.2 list is complete
    s = set(H.p53_ddr_symbols())
    assert {"TP53", "CDKN1A", "TP53BP1", "USP28", "CHEK2", "ATM", "MDM2", "MDM4", "PPM1D"} <= s
    assert {"RRM2B", "ZMAT3", "DDB2"} <= s


def test_loao_never_uses_left_out_anchor():
    X, Y, _ = _anchor_data(n_anchor=4)
    # make anchor a2 carry a unique offset on g0: if it leaked, a2's own correction would remove part of it
    X = X.copy()
    X.loc[:, "g0"] = Y.loc[:, "g0"]
    X.loc["a2", "g0"] = Y.loc["a2", "g0"] + 3.0
    seen = []

    def metric(xc, a):
        seen.append(a)
        return {"err_g0": float(xc["g0"] - Y.loc[a, "g0"])}

    res, prof, fits = H.loao(X, Y, metric)
    assert seen == list(X.index)
    for a, c in fits.items():
        assert a not in c.anchors and len(c.anchors) == 3
    # a2's correction comes from a0,a1,a3 only (offset 0 on g0), so its own +3 offset is untouched
    assert abs(fits["a2"].delta["g0"]) < 0.05
    assert res.set_index("anchor").loc["a2", "err_g0"] == pytest.approx(3.0, abs=0.05)
    # whereas a fit that includes a2 moves g0 (B * 3 / 4 ~ 0.7)
    full = H.fit_genewise(X, Y)
    assert full.delta["g0"] > 0.4


# ---------------------------------------------------------------- agreement metrics
def test_deviation_corr_and_nn_rank():
    rng = np.random.default_rng(1)
    genes = _genes(500)
    Ydev = pd.DataFrame(rng.normal(0, 1, (30, 500)), index=[f"L{i}" for i in range(30)], columns=genes)
    m = Ydev.mean()
    x = Ydev.loc["L7"] * 2 + 0.3 * rng.normal(0, 1, 500)
    r = H.deviation_corr_all(x, Ydev, m)
    assert r.idxmax() == "L7"
    assert r["L7"] == pytest.approx(H.deviation_corr(x, Ydev.loc["L7"], m))
    nn = H.nn_rank(x, Ydev, m, "L7")
    assert nn["nn_rank"] == 1 and nn["nn_recall1"] == 1
    # NaNs are handled pairwise
    Yn = Ydev.copy()
    Yn.iloc[3, :100] = np.nan
    r2 = H.deviation_corr_all(x, Yn, m)
    assert r2["L3"] == pytest.approx(H.deviation_corr(x, Yn.loc["L3"], m))


def test_auprc_and_recall_at_fpr():
    s = np.array([5, 4, 3, 2, 1.0])
    assert H.auprc(s, np.array([1, 1, 0, 0, 0])) == pytest.approx(1.0)
    assert H.auprc(s, np.array([0, 0, 0, 1, 1])) == pytest.approx((1 / 4 + 2 / 5) / 2)
    sc = pd.Series(np.arange(100, dtype=float), index=_genes(100))
    pos, neg = pd.Index(_genes(100)[90:]), pd.Index(_genes(100)[:80])
    assert H.recall_at_fpr(sc, pos, neg, fpr=0.05) == 1.0


def test_rank_scores_direction():
    s = pd.Series([10.0, -5.0, 0.0], index=["a", "b", "c"])
    assert H.rank_scores(s, "high").idxmin() == "a"  # BF high = essential -> lowest rank
    assert H.rank_scores(s, "low").idxmin() == "b"
    with pytest.raises(ValueError):
        H.rank_scores(s, "unknown")
