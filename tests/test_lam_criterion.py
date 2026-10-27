"""D20c lambda-criterion diagnostic: LOO-MSE and grouped inner-CV r (toy data, frame mode)."""
import numpy as np
import pandas as pd
import pytest

from essml.models.kernel_ridge import KernelRidgeM1, LAM_GRID, fit_pattern_ridge


def _toy(n=160, p=40, G=30, signal=True, seed=0):
    rng = np.random.default_rng(seed)
    ids = pd.Index([f"L{i}" for i in range(n)])
    X = pd.DataFrame(rng.normal(size=(n, p)), index=ids, columns=[f"x{j}" for j in range(p)])
    B = rng.normal(size=(p, G)) * (0.3 if signal else 0.0)
    Y = pd.DataFrame(X.to_numpy() @ B + rng.normal(size=(n, G)), index=ids,
                     columns=[str(100 + g) for g in range(G)])
    return X, Y


def _fit(criterion, X, Y):
    m = KernelRidgeM1(modalities=["X"], residual="B0", lam_criterion=criterion, min_rows=20)
    return m.fit(Y, lineage=None, X=X)


def test_noise_loo_r_picks_smallest_lambda_but_mse_and_inner_cv_do_not():
    X, Y = _toy(signal=False)
    lam = {c: _fit(c, X, Y).lam_by_group_ for c in ("loo_r", "loo_mse", "inner_cv_r")}
    # the D3 artefact on the B0 path: LOO r -> -1 as lambda -> inf, so pure noise picks the minimum
    assert all(v == LAM_GRID[0] for v in lam["loo_r"].values())
    # LOO-MSE and inner-CV r prefer strong shrinkage on noise (MSE flat beyond ~1e3 here)
    assert all(v >= 1e3 for v in lam["loo_mse"].values())
    assert all(v >= 1e3 for v in lam["inner_cv_r"].values())


def test_signal_all_criteria_predict_well():
    X, Y = _toy(signal=True, seed=1)
    tr, te = Y.index[:120], Y.index[120:]
    for c in ("loo_r", "loo_mse", "inner_cv_r"):
        m = _fit(c, X.loc[tr], Y.loc[tr])
        P = m.predict(te, X=X.loc[te])
        r = [np.corrcoef(P[g], Y.loc[te, g])[0, 1] for g in Y.columns]
        assert np.mean(r) > 0.5, c
        if c != "loo_r":
            assert m.report_["lam_criterion"] == c and m.report_["score_by_group"]


def test_default_criterion_unchanged():
    X, Y = _toy(seed=2)
    a = KernelRidgeM1(modalities=["X"], residual="B0", min_rows=20).fit(Y, X=X)
    b = _fit("loo_r", X, Y)
    assert a.lam_by_group_ == b.lam_by_group_
    assert "lam_criterion" not in a.report_


def test_inner_cv_validation_targets_do_not_reach_their_own_fit():
    """Perturbing the regression targets of inner-validation fold b_0 leaves split 0's scores
    unchanged (its fit uses rows a_0 only) and changes the other splits (b_0 trains them)."""
    from essml.models.kernel_ridge import cluster_patterns, inner_cv_r_grid
    X, Y = _toy(seed=3)
    Z = X.to_numpy()
    K = Z @ Z.T
    Yv = Y.to_numpy()
    R = Yv - Yv.mean(0)
    idx = np.arange(len(Y))
    splits = [(idx[idx % 4 != s], idx[idx % 4 == s]) for s in range(4)]
    cl, _, _ = cluster_patterns(~np.isnan(R), 1, 20)
    g1 = inner_cv_r_grid(K, R, Yv, cl, splits, LAM_GRID)
    R2 = R.copy()
    b0 = splits[0][1]
    R2[b0] += np.random.default_rng(9).normal(size=(len(b0), R.shape[1])) * 5
    g2 = inner_cv_r_grid(K, R2, Yv, cl, splits, LAM_GRID)  # same truth, perturbed fit targets
    assert np.allclose(g1[0], g2[0], rtol=0, atol=1e-10)
    for s in (1, 2, 3):
        assert not np.allclose(g1[s], g2[s])


def test_inner_cv_requires_b0_residual_and_splits():
    with pytest.raises(ValueError):
        KernelRidgeM1(residual="B1", lam_criterion="inner_cv_r")
    X, Y = _toy(seed=4)
    with pytest.raises(ValueError):
        fit_pattern_ridge(X.to_numpy() @ X.to_numpy().T, Y.to_numpy(), Y.to_numpy(),
                          np.zeros(Y.shape[1], int), criterion="inner_cv_r")
