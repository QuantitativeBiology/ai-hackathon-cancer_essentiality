"""G3 nonlinear-track models on toy data (frame mode, in-process): shapes, NaN masking, and that
no held-out row (target or feature) reaches any fit."""
import numpy as np
import pandas as pd
import pytest

import run
from essml.models import factorised as FZ
from essml.models.g3_common import TrainPCA, inner_split, masked_r, top_abs_corr
from essml.models.lgbm import LGBMResidualM3, screen_top_k
from essml.models.published import ElasticNetR2, TopKCorrR4

torch = pytest.importorskip("torch")


def _toy(n=120, g=12, p=30, seed=0):
    rng = np.random.default_rng(seed)
    ids = pd.Index([f"L{i}" for i in range(n)])
    X = pd.DataFrame(rng.normal(size=(n, p)), index=ids, columns=[f"x{j}" for j in range(p)])
    W = rng.normal(size=(4, g))
    Y = pd.DataFrame(np.tanh(X.iloc[:, :4].to_numpy()) @ W * 0.5 + rng.normal(size=(n, g)) * 0.5,
                     index=ids, columns=[f"g{j}" for j in range(g)])
    Y = Y.mask(rng.random(Y.shape) < 0.08)
    lineage = pd.Series(rng.choice(list("abc"), size=n), index=ids)
    groups = pd.Series([f"P{i // 2}" for i in range(n)], index=ids)
    folds = pd.Series(pd.factorize(groups)[0] % 3, index=ids)
    return Y, X, lineage, groups, folds


M1P = {"modalities": ["E"], "residual": "B0", "min_rows": 20}
MODELS = {
    "RidgePC": lambda: FZ.RidgePC(n_pc=10),
    "MLP": lambda: FZ.NeuralPC(n_pc=10, hidden=16, max_epochs=8, patience=3,
                               grid={"whiten": [False], "dropout": [0.1]}),
    "M4": lambda: FZ.FactorisedM4(n_pc=10, hidden=16, emb_dim=4, max_epochs=8, patience=3,
                                  grid={"whiten": [False], "dropout": [0.1]}),
    "M3": lambda: LGBMResidualM3(n_jobs=1, k_screen=5, min_rows=20, m1_params=M1P,
                                 num_boost_round=20, early_stopping=5, cache_inner=False),
    "R2": lambda: ElasticNetR2(n_jobs=1, min_rows=20, n_alphas=5, inner_folds=3),
    "R4": lambda: TopKCorrR4(n_jobs=1, min_rows=20, ks=(2, 5)),
}


@pytest.mark.parametrize("key", list(MODELS))
def test_shapes_and_finite(key):
    Y, X, lin, grp, folds = _toy()
    tr, ho = Y.index[folds != 0], Y.index[folds == 0]
    m = MODELS[key]().fit(Y.loc[tr], lineage=lin.loc[tr], X=X, groups=grp.loc[tr].to_numpy())
    P = m.predict(ho, lineage=lin.loc[ho], X=X)
    assert P.shape == (len(ho), Y.shape[1])
    assert list(P.columns) == list(Y.columns) and list(P.index) == list(ho)
    assert np.isfinite(P.to_numpy()).all()


@pytest.mark.parametrize("key", list(MODELS))
def test_no_fit_sees_held_out_and_held_out_invariance(key):
    """Through the launcher: no fit receives a held-out ID, and perturbing held-out TARGETS and the
    held-out FEATURES of one line leaves the other held-out lines' predictions unchanged."""
    Y, X, lin, grp, folds = _toy()
    cfg = {"model": key if key in ("M3", "R2", "R4") else {"RidgePC": "M4pre_RidgePC",
                                                           "MLP": "M4pre_MLP", "M4": "M4"}[key],
           "params": {}, "grid": {}, "inner_folds": 3, "inner_min_obs": 5}
    orig = run.MODELS[cfg["model"]]
    run.MODELS[cfg["model"]] = lambda **kw: MODELS[key]()
    try:
        run.FIT_LOG.clear()
        a = run.run_cv(Y, folds, grp, lin, X, cfg, Y.columns, set(), seed=7, folds_to_run=[0],
                       log=lambda *x: None)
        for e in run.FIT_LOG:
            assert not set(folds.index[folds == e["outer_fold"]]).intersection(e["ids"])
        ho = folds.index[folds == 0]
        Y2, X2 = Y.copy(), X.copy()
        Y2.loc[ho] = np.random.default_rng(5).normal(size=(len(ho), Y.shape[1])) * 10
        Y2.loc[ho[:5], "g0"] = np.nan
        X2.loc[ho[0]] = X2.loc[ho[0]] * 50 + 3  # one held-out line's features
        b = run.run_cv(Y2, folds, grp, lin, X2, cfg, Y.columns, set(), seed=7, folds_to_run=[0],
                       log=lambda *x: None)
    finally:
        run.MODELS[cfg["model"]] = orig
    pd.testing.assert_frame_equal(a["pred"].loc[ho[1:]], b["pred"].loc[ho[1:]])


def test_ridge_patterns_equals_ridge_on_observed_rows():
    rng = np.random.default_rng(1)
    P = rng.normal(size=(40, 5))
    Y = rng.normal(size=(40, 3))
    Y[[1, 7, 9], 1] = np.nan
    coef, icpt = FZ.ridge_patterns(P, Y, [2.0])
    rows = ~np.isnan(Y[:, 1])
    Pr, yr = P[rows], Y[rows, 1]
    Pc = Pr - Pr.mean(0)
    c = np.linalg.solve(Pc.T @ Pc + 2.0 * np.eye(5), Pc.T @ (yr - yr.mean()))
    np.testing.assert_allclose(coef[0][:, 1], c, atol=1e-10)
    np.testing.assert_allclose(icpt[0, 1], yr.mean() - Pr.mean(0) @ c, atol=1e-10)


def test_masked_mse_ignores_masked_cells():
    out = torch.tensor([[1.0, 2.0], [3.0, 4.0]])
    y = torch.tensor([[1.5, 1e6], [2.0, -1e6]])
    m = torch.tensor([[1.0, 0.0], [1.0, 0.0]])
    assert float(FZ.masked_mse(out, y, m)) == pytest.approx((0.25 + 1.0) / 2)


def test_neural_nan_targets_train_and_predict_finite():
    Y, X, lin, grp, folds = _toy()
    Y.iloc[:100, 0] = np.nan  # a gene with few observed rows
    m = MODELS["MLP"]().fit(Y, lineage=lin, X=X, groups=grp.to_numpy())
    assert np.isfinite(m.predict(Y.index[:5], X=X).to_numpy()).all()


def test_train_pca_uses_only_given_rows():
    rng = np.random.default_rng(2)
    Z = rng.normal(size=(50, 8))
    a = TrainPCA(3).fit(Z[:40])
    Z2 = Z.copy()
    Z2[40:] *= 100
    b = TrainPCA(3).fit(Z2[:40])
    np.testing.assert_allclose(a.transform(Z[:5]), b.transform(Z2[:5]), atol=1e-5)


def test_inner_split_respects_groups():
    ids = pd.Index([f"L{i}" for i in range(60)])
    groups = np.array([f"P{i // 3}" for i in range(60)])
    a, b = inner_split(ids, groups, None, seed=1)
    assert not set(groups[a]).intersection(groups[b])
    assert len(a) + len(b) == 60


def test_screen_top_k_matches_bruteforce_with_nan():
    rng = np.random.default_rng(3)
    Z = rng.normal(size=(60, 20))
    R = rng.normal(size=(60, 4)) + Z[:, [3, 5, 7, 9]]
    R[:10, 2] = np.nan
    rows = [~np.isnan(R[:, j]) for j in range(4)]
    sel = screen_top_k(Z, R, 3, rows)
    for j in range(4):
        exp = top_abs_corr(Z[rows[j]], R[rows[j], j], 3)
        assert list(sel[j]) == list(exp)


def test_oracle_select_is_optimistic():
    """oracle_select picks per gene the best outer-fold candidate, so its r >= the nested choice."""
    Y, X, lin, grp, folds = _toy(n=150)
    cfg = {"model": "R4", "params": {}, "grid": {}, "inner_folds": 3}
    orig = run.MODELS["R4"]
    run.MODELS["R4"] = lambda **kw: TopKCorrR4(n_jobs=1, min_rows=20, ks=(2, 5), nested=False)
    try:
        res = run.run_cv(Y, folds, grp, lin, X, cfg, Y.columns, set(), seed=7, log=lambda *x: None)
        pred0 = res["pred"].copy()
        pred = TopKCorrR4.oracle_select(res, Y, folds)
    finally:
        run.MODELS["R4"] = orig
    r0 = np.nanmean([np.nanmean(masked_r(pred0.loc[folds.index[folds == f]].to_numpy(),
                                         Y.loc[folds.index[folds == f]].to_numpy())) for f in range(3)])
    r1 = np.nanmean([np.nanmean(masked_r(pred.loc[folds.index[folds == f]].to_numpy(),
                                         Y.loc[folds.index[folds == f]].to_numpy())) for f in range(3)])
    assert r1 >= r0 - 1e-12
