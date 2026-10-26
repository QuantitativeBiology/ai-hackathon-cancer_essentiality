"""G3 report columns (D16/D17): within-group centring and the Nadeau-Bengio SE."""
import numpy as np
import pandas as pd

import g3_report as R

rng = np.random.default_rng(5)


def test_group_centre_matches_lineage_partialling_and_kills_group_means():
    n, G = 60, 7
    lab = rng.choice(list("abc"), size=n)
    T = rng.normal(size=(n, G)) + np.array([{"a": 1, "b": -1, "c": 3}[l] for l in lab])[:, None]
    T[rng.random(T.shape) < 0.1] = np.nan
    P = rng.normal(size=(n, G))
    Pc, Tc = R.group_centre(P, T, lab)
    # = OLS residual on [1, lineage dummies] over each gene's observed rows (neg_genes.partial_out)
    D = pd.get_dummies(lab).to_numpy(float)
    for g in range(G):
        o = ~np.isnan(T[:, g])
        for A, Ac in [(T, Tc), (P, Pc)]:
            beta, *_ = np.linalg.lstsq(D[o], A[o, g], rcond=None)
            assert np.allclose(Ac[o, g], A[o, g] - D[o] @ beta)
    # a pure group-mean predictor (B1-like) has zero within-group skill
    B = np.array([{"a": 1.0, "b": -1.0, "c": 3.0}[l] for l in lab])[:, None].repeat(G, 1)
    Bc, _ = R.group_centre(B, T, lab)
    assert np.allclose(Bc, 0)


def test_nadeau_bengio_se():
    d = np.array([0.01, 0.02, 0.0, 0.015, 0.005])
    se = d.std(ddof=1) / np.sqrt(5)
    assert np.isclose(R.nb_se(d, 0.25), se * np.sqrt(0.45 / 0.2))
    out = R.paired_nb(d + 1, np.ones(5), 0.25)
    assert np.isclose(out["se"], se) and out["se_nb"] > out["se"] and out["p_nb"] > out["p"]
    folds = pd.Series([0] * 20 + [1] * 20 + [2] * 20 + [3] * 20 + [4] * 20)
    assert np.isclose(R.test_train_ratio(folds), 0.25)
