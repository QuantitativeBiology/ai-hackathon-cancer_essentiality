"""Known-answer tests for essml.metrics (DESIGN §8)."""
import numpy as np
import pytest
from sklearn.metrics import average_precision_score

from essml import metrics as M

rng = np.random.default_rng(0)


def _toy(n=120, g=30, nan_frac=0.1):
    T = rng.normal(size=(n, g)) - 0.4
    T[rng.random((n, g)) < nan_frac] = np.nan
    return T


def test_perfect_and_anticorrelated():
    T = _toy()
    P = np.where(np.isnan(T), 0.0, T)  # pred must be finite; NaN-true entries are ignored
    assert np.allclose(M.per_gene_pearson(P, T, min_obs=10), 1.0)
    assert np.allclose(M.per_gene_spearman(P, T, min_obs=10), 1.0)
    assert np.allclose(M.per_line_pearson(P, T), 1.0)
    assert np.allclose(M.per_gene_pearson(-P, T, min_obs=10), -1.0)
    assert np.allclose(M.per_gene_spearman(-P, T, min_obs=10), -1.0)
    ap = M.per_line_ap(P, T)
    assert np.allclose(ap[np.isfinite(ap)], 1.0) and np.isfinite(ap).sum() > 100
    g = M.per_gene_ap_norm(P, T)
    ok = np.isfinite(g["ap"])
    assert ok.sum() > 20 and np.allclose(g["ap"][ok], 1) and np.allclose(g["nap"][ok], 1)
    assert np.allclose(g["auroc"][ok], 1)


def test_nan_values_in_unobserved_entries_do_not_matter():
    T = _toy()
    P1 = np.where(np.isnan(T), 0.0, T + rng.normal(size=T.shape))
    P2 = np.where(np.isnan(T), 1e6, P1)
    assert np.allclose(M.per_gene_pearson(P1, T, 10), M.per_gene_pearson(P2, T, 10))


def test_constant_predictor_r0_ap_pi():
    T = _toy()
    P = np.ones_like(T) * 0.123456789  # constant per gene
    assert np.all(M.per_gene_pearson(P, T, 10) == 0.0)
    assert np.all(M.per_gene_spearman(P, T, 10) == 0.0)
    g = M.per_gene_ap_norm(P, T)
    ok = np.isfinite(g["ap"])
    assert ok.sum() > 20
    assert np.array_equal(g["ap"][ok], g["prev"][ok])  # AP == π exactly
    assert np.all(g["nap"][ok] == 0.0) and np.all(g["auroc"][ok] == 0.5)
    # per-line: constant across genes -> AP == per-line prevalence
    ap = M.per_line_ap(P, T)
    obs = ~np.isnan(T)
    pi = ((T < -0.5) & obs).sum(1) / obs.sum(1)
    assert np.allclose(ap[np.isfinite(ap)], pi[np.isfinite(ap)], rtol=0, atol=1e-15)


def test_column_mean_like_constant_is_exactly_zero():
    # gene mean repeated over lines: a rounding-prone constant still gives exactly 0
    T = _toy()
    mu = np.nanmean(T, 0) + 1e-9 / 3
    P = np.tile(mu, (T.shape[0], 1))
    assert np.all(M.per_gene_pearson(P, T, 10) == 0.0)


def test_constant_truth_gives_zero():
    T = np.zeros((50, 3)); T[:, 1] = 2.0; T[:, 2] = rng.normal(size=50)
    P = rng.normal(size=(50, 3))
    r = M.per_gene_pearson(P, T, 10)
    assert r[0] == 0 and r[1] == 0 and r[2] != 0


def test_random_predictor_ap_close_to_prevalence():
    n, g = 4000, 5
    T = rng.normal(size=(n, g)) - 0.3
    P = rng.normal(size=(n, g))
    res = M.per_gene_ap_norm(P, T)
    assert np.all(np.abs(res["ap"] - res["prev"]) < 0.03)
    assert np.all(np.abs(res["auroc"] - 0.5) < 0.04)


def test_min_obs_excludes():
    T = _toy(n=60, g=4, nan_frac=0)
    T[:30, 0] = np.nan
    P = rng.normal(size=T.shape)
    r = M.per_gene_pearson(P, T, min_obs=40)
    assert np.isnan(r[0]) and np.isfinite(r[1:]).all()


def test_nan_pred_raises():
    T = _toy(); P = np.where(np.isnan(T), 0, T); P[0, 0] = np.nan; T[0, 0] = 1.0
    with pytest.raises(ValueError):
        M.per_gene_pearson(P, T, 10)


def test_against_numpy_and_sklearn():
    from scipy.stats import spearmanr
    T = _toy(n=80, g=6)
    P = np.where(np.isnan(T), 0, T + rng.normal(size=T.shape))
    P = np.round(P, 1)  # ties for Spearman
    r = M.per_gene_pearson(P, T, 10); rho = M.per_gene_spearman(P, T, 10)
    for j in range(6):
        o = ~np.isnan(T[:, j])
        assert np.isclose(r[j], np.corrcoef(P[o, j], T[o, j])[0, 1])
        assert np.isclose(rho[j], spearmanr(P[o, j], T[o, j])[0])
    ap = M.per_line_ap(P, T)
    o = ~np.isnan(T[3])
    assert np.isclose(ap[3], average_precision_score(T[3, o] < -0.5, -P[3, o]))


def test_p2b_eligibility():
    T = np.full((20, 1), 0.0); T[:2, 0] = -1  # only 2 positives
    g = M.per_gene_ap_norm(rng.normal(size=T.shape), T)
    assert np.isnan(g["ap"][0]) and g["prev"][0] == 0.1


def test_recall_at_fpr():
    neg_scores = np.arange(100, dtype=float)          # 100 negatives, scores 0..99
    pos_scores = np.array([96.5, 95.5, 94.9, 50, 200])
    s = np.r_[neg_scores, pos_scores]
    pos = np.r_[np.zeros(100, bool), np.ones(5, bool)]
    # fpr 5% -> k = 5 allowed FPs; threshold = 6th highest negative = 94 -> called: >94
    assert M.recall_at_fpr(s, pos, ~pos, 0.05) == pytest.approx(4 / 5)
    assert M.recall_at_fpr(s, pos, ~pos, 0.0) == pytest.approx(1 / 5)   # > 99


def test_deviation_corr():
    gm = rng.normal(size=50)
    T = gm + rng.normal(size=50)
    assert M.deviation_corr(gm, T, gm) == 0.0                      # B0: deviation ≡ 0
    assert M.deviation_corr(T, T, gm) == pytest.approx(1.0)
    assert M.deviation_corr(2 * gm - T, T, gm) == pytest.approx(-1.0)
