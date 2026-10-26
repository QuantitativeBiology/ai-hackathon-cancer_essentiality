"""M1 pattern-wise kernel ridge: known answers against brute-force per-gene ridge (DESIGN §8).

The pattern fit (one centred sub-kernel of a single fold-level K per NaN pattern) must equal an
explicit primal ridge with an unpenalised intercept fitted per gene on that gene's observed rows,
including closed-form leave-one-out and leave-group-out residuals. NaN rows hold wildly different
values in the toy so that centring with full-fold means (a common bug) would fail.
"""
import numpy as np
import pandas as pd
import pytest

import run
from essml.models.kernel_ridge import (KernelRidgeM1, KernelRidgeR1, _colwise_r, cluster_patterns,
                                       fit_pattern_ridge, loo_cluster, predict_pattern_ridge)

rng = np.random.default_rng(11)


def _ridge(Z, y, lam):
    zm, ym = Z.mean(0), y.mean()
    Zc = Z - zm
    beta = np.linalg.solve(Zc.T @ Zc + lam * np.eye(Z.shape[1]), Zc.T @ (y - ym))
    return lambda Znew: ym + (Znew - zm) @ beta


def _toy(n=60, p=8, G=5):
    Z = rng.normal(size=(n, p)) * rng.uniform(0.5, 2, size=p)
    B = rng.normal(size=(p, G))
    Y = Z @ B + rng.normal(size=(n, G))
    Y[:, 1] += 50.0  # large offset
    Y[:20] += 30.0   # rows that are NaN for genes 2-3: full-fold centring would be visibly wrong
    obs = np.ones((n, G), bool)
    obs[:20, 2] = False   # pattern A (genes 2, 3)
    obs[:20, 3] = False
    obs[40:, 4] = False   # pattern B (gene 4)
    Yn = np.where(obs, Y, np.nan)
    return Z, Yn, obs


@pytest.mark.parametrize("lam", [0.5, 5.0, 50.0])
def test_pattern_fit_equals_per_gene_ridge_with_loo(lam):
    Z, Y, obs = _toy()
    n, G = Y.shape
    K = Z @ Z.T
    sol = fit_pattern_ridge(K, Y, Y, np.zeros(G, int), lam_grid=np.array([lam]), merge_frac=0.0,
                            min_rows=5)
    assert sol["report"]["n_clusters"] == 3
    Znew = rng.normal(size=(7, Z.shape[1]))
    P = predict_pattern_ridge(Znew @ Z.T, sol["clusters"], G)
    for g in range(G):
        rows = np.flatnonzero(obs[:, g])
        f = _ridge(Z[rows], Y[rows, g], lam)
        assert np.allclose(P[:, g], f(Znew), atol=1e-8)
        # brute-force LOO r of y
        loo = np.array([_ridge(np.delete(Z[rows], i, 0), np.delete(Y[rows, g], i), lam)(Z[rows[i]][None])[0]
                        for i in range(len(rows))])
        r_bf = np.corrcoef(loo, Y[rows, g])[0, 1]
        assert np.isclose(sol["r_grid"][0, g], r_bf, atol=1e-8)


def test_loo_residuals_and_group_loo_match_brute_force():
    Z, Y, _ = _toy()
    rows = np.arange(20, 60)  # fully observed rows
    y = Y[rows, :2]
    Zs = Z[rows]
    K = Zs @ Zs.T
    kbar = K.mean(0)
    Kc = K - kbar[None] - kbar[:, None] + kbar.mean()
    groups = np.array([i // 3 for i in range(len(rows))])  # patient triplets
    lam = 3.0
    d, V, U, rg = loo_cluster(Kc, y - y.mean(0), y, groups, [lam])
    pred = np.zeros_like(y)
    for gid in np.unique(groups):
        te = groups == gid
        for j in range(2):
            pred[te, j] = _ridge(Zs[~te], y[~te, j], lam)(Zs[te])
    assert np.allclose(rg[0], _colwise_r(pred, y), atol=1e-8)
    d, V, U, rg1 = loo_cluster(Kc, y - y.mean(0), y, None, [lam])
    pred1 = np.array([[_ridge(np.delete(Zs, i, 0), np.delete(y[:, j], i), lam)(Zs[i][None])[0]
                       for j in range(2)] for i in range(len(rows))])
    assert np.allclose(rg1[0], _colwise_r(pred1, y), atol=1e-8)


def test_pattern_merging_and_fallback():
    n = 200
    obs = np.ones((n, 6), bool)
    obs[0, 1] = False          # differs by 1 row (<= 1% = 2 rows) -> merged with the full pattern
    obs[:100, 2] = False       # differs by 100 rows -> own cluster
    obs[:100, 3] = False
    obs[:170, 4] = False       # 30 rows < min_rows -> fallback
    obs[5:9, 5] = False        # differs by 4 rows -> own cluster
    cl, fb, rep = cluster_patterns(obs, tol=2, min_rows=50)
    assert rep["n_patterns"] == 5 and rep["n_patterns_merged"] == 1
    assert sorted(fb.tolist()) == [4]
    byg = {int(g): c for c in cl for g in c["genes"]}
    assert byg[0] is byg[1] and len(byg[0]["rows"]) == n - 1   # merged: intersection rows
    assert len(byg[2]["rows"]) == 100 and byg[2] is byg[3]
    assert len(byg[5]["rows"]) == n - 4
    # merged gene equals per-gene ridge on the intersection rows
    Z = rng.normal(size=(n, 4))
    Y = Z @ rng.normal(size=(4, 6)) + rng.normal(size=(n, 6))
    Yn = np.where(obs, Y, np.nan)
    sol = fit_pattern_ridge(Z @ Z.T, Yn, Yn, np.zeros(6, int), lam_grid=np.array([1.0]),
                            merge_frac=0.01, min_rows=50)
    P = predict_pattern_ridge(Z[:3] @ Z.T, sol["clusters"], 6)
    rows = np.arange(1, n)
    assert np.allclose(P[:, 0], _ridge(Z[rows], Y[rows, 0], 1.0)(Z[:3]))
    assert np.allclose(P[:, 4], 0.0)  # fallback genes: residual model predicts 0


def test_grouped_lambda_choice_is_argmax_of_mean_loo_r():
    Z, Y, _ = _toy(n=80, G=6)
    gg = np.array([0, 0, 1, 1, 1, 0])
    grid = np.logspace(-1, 3, 9)
    sol = fit_pattern_ridge(Z @ Z.T, Y, Y, gg, lam_grid=grid, merge_frac=0.0, min_rows=5)
    for lab in (0, 1):
        sc = sol["r_grid"][:, gg == lab].mean(1)
        assert sol["lam_by_group"][lab] == grid[int(np.argmax(sc))]
        assert np.all(sol["lam"][gg == lab] == sol["lam_by_group"][lab])


def test_model_class_frame_mode_matches_core_and_offsets():
    n, p, G = 90, 6, 4
    ids = pd.Index([f"L{i}" for i in range(n)])
    X = pd.DataFrame(rng.normal(size=(n, p)), index=ids, columns=[f"x{j}" for j in range(p)])
    Y = pd.DataFrame(X.to_numpy() @ rng.normal(size=(p, G)) + rng.normal(size=(n, G)), index=ids,
                     columns=[f"g{j}" for j in range(G)])
    lin = pd.Series(rng.choice(list("abc"), size=n), index=ids)
    tr, te = ids[:70], ids[70:]
    m = KernelRidgeM1(modalities=["X"], lam_grid=[10.0], kappa=5.0, min_rows=10)
    m.fit(Y.loc[tr], lineage=lin, X=X.loc[tr], groups=np.arange(70))
    P = m.predict(te, lineage=lin, X=X.loc[te])
    assert P.shape == (20, G) and np.isfinite(P.to_numpy()).all()
    # R1 on raw y equals a per-gene primal ridge on standardised (clipped) features
    r1 = KernelRidgeR1(modalities=["X"], lam_grid=[10.0], min_rows=10)
    r1.fit(Y.loc[tr], X=X.loc[tr], groups=np.arange(70))
    Pr = r1.predict(te, X=X.loc[te])
    Xt = X.loc[tr].to_numpy()
    mu, sd = Xt.mean(0), Xt.std(0)
    Zt, Zte = (Xt - mu) / sd, (X.loc[te].to_numpy() - mu) / sd
    for g in range(G):
        assert np.allclose(Pr.iloc[:, g], _ridge(Zt, Y.loc[tr].iloc[:, g].to_numpy(), 10.0)(Zte))
    m.compact()
    assert np.allclose(m.predict(te, lineage=lin, X=X.loc[te]), P, atol=1e-4)


def _toy_run(n=150, g=30, p=40):
    ids = pd.Index([f"L{i}" for i in range(n)])
    lineage = pd.Series(rng.choice(list("abcde"), size=n), index=ids)
    groups = pd.Series([f"P{i // 2}" for i in range(n)], index=ids)
    folds = pd.Series(pd.factorize(groups)[0] % 5, index=ids)
    X = pd.DataFrame(rng.normal(size=(n, p)), index=ids, columns=[f"x{j}" for j in range(p)])
    Y = pd.DataFrame(rng.normal(size=(n, g)), index=ids, columns=[f"g{j}" for j in range(g)])
    Y.iloc[:, :5] += X.iloc[:, :5].to_numpy()
    Y = Y.mask(rng.random(Y.shape) < 0.03)
    return Y, folds, groups, lineage, X


@pytest.mark.parametrize("key", ["M1", "R1"])
def test_launcher_no_leakage_and_held_out_invariance(key):
    Y, folds, groups, lineage, X = _toy_run()
    cfg = {"model": key, "params": {"modalities": ["X"], "min_rows": 20, "lam_grid": [1, 10, 100]},
           "inner_folds": 3, "inner_min_obs": 5, "reference_baselines": {"B0": {}}}
    run.FIT_LOG.clear()
    a = run.run_cv(Y, folds, groups, lineage, X, cfg, Y.columns, set(), seed=7, log=lambda *x: None)
    for e in run.FIT_LOG:
        assert not set(folds.index[folds == e["outer_fold"]]).intersection(e["ids"])
    Y2 = Y.copy()
    ho = folds.index[folds == 0]
    Y2.loc[ho] = rng.normal(size=(len(ho), Y.shape[1])) * 10
    b = run.run_cv(Y2, folds, groups, lineage, X, cfg, Y.columns, set(), seed=7, log=lambda *x: None)
    pd.testing.assert_frame_equal(a["pred"].loc[ho], b["pred"].loc[ho])
    # the model learns the planted signal (positive control)
    r = [np.corrcoef(a["pred"].loc[ho, f"g{j}"], Y.loc[ho, f"g{j}"].fillna(0))[0, 1] for j in range(5)]
    assert np.mean(r) > 0.3


def test_m2_batched_ridge_equals_per_gene_ridge():
    from essml.models.m2_local import batched_ridge
    p, n, G = 4, 70, 9
    Xs = rng.normal(size=(p, n, G)) * rng.uniform(0.2, 3, size=(p, 1, 1))
    Xs[2, :, 3] = 1.7  # a constant feature for gene 3 -> coefficient 0
    Y = np.einsum("ing,i->ng", Xs, rng.normal(size=p)) + rng.normal(size=(n, G))
    Y[rng.random(Y.shape) < 0.2] = np.nan
    lams = rng.uniform(0.5, 20, size=G)
    beta, icpt, ok = batched_ridge(Xs, Y, lams, min_rows=10)
    for g in range(G):
        o = ~np.isnan(Y[:, g])
        X = Xs[:, o, g].T
        mx, sd = X.mean(0), X.std(0)
        sd = np.where(sd > 1e-8, sd, np.inf)
        Zs = (X - mx) / sd
        bs = np.linalg.solve(Zs.T @ Zs + lams[g] * np.eye(p), Zs.T @ (Y[o, g] - Y[o, g].mean()))
        assert np.allclose(beta[g], bs / sd, atol=1e-8)
        assert np.isclose(icpt[g], Y[o, g].mean() - mx @ (bs / sd))
    assert beta[3, 2] == 0.0 and ok.all()


def test_m2_recovers_planted_local_signal_and_uses_training_rows_only():
    from essml.models.m2_local import LocalStageM2
    from essml.priors import LocalFeatureBuilder
    n, G = 160, 8
    ids = pd.Index([f"L{i}" for i in range(n)])
    genes = [str(100 + j) for j in range(G)]
    E = pd.DataFrame(rng.normal(3, 1, size=(n, G)), index=ids, columns=genes)
    cn = pd.DataFrame(1.0, index=ids, columns=genes)
    dam = pd.DataFrame(0.0, index=ids, columns=genes)
    top = pd.DataFrame({"par1": [genes[(j + 1) % G] for j in range(G)], "par2": np.nan, "par3": np.nan,
                        "par1_score": 50.0, "par2_score": np.nan, "par3_score": np.nan,
                        "family_size": 1, "family_size_id20": 1, "max_identity": 50.0},
                       index=pd.Index(genes, name="entrez"))
    partners = pd.DataFrame({"a": genes, "b": genes[1:] + genes[:1], "source": "union"})
    builder = LocalFeatureBuilder(E, cn, dam, top, partners)
    Xm1 = pd.DataFrame(rng.normal(size=(n, 5)), index=ids, columns=[f"x{j}" for j in range(5)])
    Y = pd.DataFrame(-0.8 * E.to_numpy() + 0.3 * rng.normal(size=(n, G)), index=ids, columns=genes)
    Y.iloc[:5, 0] = np.nan
    lin = pd.Series(rng.choice(list("ab"), size=n), index=ids)
    tr, te = ids[:120], ids[120:]
    m = LocalStageM2(m1_params={"modalities": ["X"], "lam_grid": [10.0, 100.0], "min_rows": 20},
                     priors="none", inner_folds=3, min_rows=20, cache_inner=False)
    m.fit(Y.loc[tr], lineage=lin, X={"_m1": Xm1.loc[tr], "_local": builder},
          groups=np.arange(120))
    P = m.predict(te, lineage=lin, X={"_m1": Xm1.loc[te], "_local": builder})
    r = [np.corrcoef(P.loc[te, g], Y.loc[te, g])[0, 1] for g in genes]
    assert np.mean(r) > 0.8
    # own_expr and expr_hinge are near-collinear here (TPM >> 1): their summed slope ~ -0.8
    assert m.feature_names_[1] == "own_expr" and m.feature_names_[3] == "expr_hinge"
    assert -0.95 < (m.beta_[:, 1] + m.beta_[:, 3]).mean() < -0.6
    # held-out targets cannot matter: fit sees training rows only
    Y2 = Y.copy()
    Y2.loc[te] = 99.0
    m2 = LocalStageM2(m1_params={"modalities": ["X"], "lam_grid": [10.0, 100.0], "min_rows": 20},
                      priors="none", inner_folds=3, min_rows=20, cache_inner=False)
    m2.fit(Y2.loc[tr], lineage=lin, X={"_m1": Xm1.loc[tr], "_local": builder}, groups=np.arange(120))
    pd.testing.assert_frame_equal(P, m2.predict(te, lineage=lin, X={"_m1": Xm1.loc[te], "_local": builder}))
    # gene chunking (memory) does not change the result
    m3 = LocalStageM2(m1_params={"modalities": ["X"], "lam_grid": [10.0, 100.0], "min_rows": 20},
                      priors="none", inner_folds=3, min_rows=20, cache_inner=False, chunk=3)
    m3.fit(Y.loc[tr], lineage=lin, X={"_m1": Xm1.loc[tr], "_local": builder}, groups=np.arange(120))
    P3 = m3.predict(te, lineage=lin, X={"_m1": Xm1.loc[te], "_local": builder})
    assert np.allclose(P, P3, atol=1e-10) and m3.lam_by_group_ == m.lam_by_group_
