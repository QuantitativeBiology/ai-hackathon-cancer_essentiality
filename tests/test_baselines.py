"""B0/B1/B2 known answers; leave-line-out lineage residual vs brute force (DESIGN §8)."""
import numpy as np
import pandas as pd
import pytest

from essml.baselines import GeneMean, KNNExpression, LineageMean, leave_line_out_lineage_residual

rng = np.random.default_rng(1)


def _toy(n=40, g=6):
    idx = [f"L{i}" for i in range(n)]
    Y = pd.DataFrame(rng.normal(size=(n, g)), index=idx, columns=[f"g{j}" for j in range(g)])
    Y = Y.mask(rng.random(Y.shape) < 0.2)
    lin = pd.Series(rng.choice(["a", "b", "c"], size=n), index=idx)
    lin.iloc[-1] = "single"  # a lineage with one line
    return Y, lin


def test_b0_nanmean():
    Y, _ = _toy()
    P = GeneMean().fit(Y).predict(["new1", "new2"])
    assert np.allclose(P.iloc[0], Y.mean()) and np.allclose(P.iloc[1], Y.mean())


def test_b1_formula_and_fallback():
    Y, lin = _toy()
    k = 3.0
    m = LineageMean(k).fit(Y, lin)
    P = m.predict(["q1", "q2"], pd.Series({"q1": "a", "q2": "unseen"}))
    ya = Y[lin == "a"]
    exp = (ya.sum() + k * Y.mean()) / (ya.notna().sum() + k)
    assert np.allclose(P.loc["q1"], exp)
    assert np.allclose(P.loc["q2"], Y.mean())                 # unseen lineage -> B0
    P0 = LineageMean(0).fit(Y, lin).predict(["q"], pd.Series({"q": "b"}))
    assert np.allclose(P0.loc["q"], Y[lin == "b"].mean())     # kappa = 0: raw lineage nanmean
    Pinf = LineageMean(1e12).fit(Y, lin).predict(["q"], pd.Series({"q": "b"}))
    assert np.allclose(Pinf.loc["q"], Y.mean())               # huge kappa -> B0


@pytest.mark.parametrize("kappa", [0.0, 0.5, 5.0])
def test_leave_line_out_residual_matches_brute_force(kappa):
    Y, lin = _toy()
    R = leave_line_out_lineage_residual(Y, lin, kappa)
    for i in Y.index:
        rest = Y.index.drop(i)
        m = LineageMean(kappa).fit(Y.loc[rest], lin.loc[rest])
        exp = Y.loc[i] - m.predict([i], lin.loc[[i]]).loc[i]
        assert np.allclose(R.loc[i].to_numpy(), exp.to_numpy(), equal_nan=True), i
    assert R.isna().equals(Y.isna())


def test_b2_knn():
    n = 30
    idx = [f"L{i}" for i in range(n)]
    X = pd.DataFrame(rng.normal(size=(n, 50)), index=idx)
    Y = pd.DataFrame(rng.normal(size=(n, 4)), index=idx, columns=list("abcd"))
    Y.iloc[:, 3] = np.nan
    Y.iloc[0, 3] = 7.0
    tr = idx[:25]
    m = KNNExpression(k=1, n_top=50).fit(Y.loc[tr], X=X.loc[tr])
    P = m.predict(tr[:5], X=X)
    assert np.allclose(P.iloc[:, :3], Y.loc[tr[:5]].iloc[:, :3])  # nearest neighbour is itself
    assert P.loc["L1", "d"] == 7.0 and P.loc["L0", "d"] == 7.0      # all-NaN neighbours -> gene mean
    Pk = m.predict(idx[25:], X=X, k=25)                               # k = all lines -> B0
    assert np.allclose(Pk.to_numpy(), np.tile(Y.loc[tr].mean().to_numpy(), (5, 1)))
    # feature scaling uses training lines only: held-out X cannot change the fit
    X2 = X.copy()
    X2.loc[idx[25:]] *= 100
    a = KNNExpression(k=3, n_top=50).fit(Y.loc[tr], X=X2.loc[tr]).predict(tr, X=X2)
    b = KNNExpression(k=3, n_top=50).fit(Y.loc[tr], X=X.loc[tr]).predict(tr, X=X)
    assert np.allclose(a, b)
