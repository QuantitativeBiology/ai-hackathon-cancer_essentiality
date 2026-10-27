"""run.py exclude_from_training and leave-one-lineage-out folds (G3; ported from the G3 leakage
review's toy check): excluded lines reach no fit, are still predicted, and no prediction (model or
B0/B1 references) depends on their targets. Also the exact M2 inner-OOF cache key."""
import numpy as np
import pandas as pd
import pytest

import run
from essml.models.m2_local import LocalStageM2

rng = np.random.default_rng(3)
n, g, p = 150, 20, 30
ids = pd.Index([f"L{i}" for i in range(n)])
lineage = pd.Series(rng.choice(list("abcde"), size=n), index=ids)
groups = pd.Series([f"P{i // 2}" for i in range(n)], index=ids)
X = pd.DataFrame(rng.normal(size=(n, p)), index=ids, columns=[f"x{j}" for j in range(p)])
Y = pd.DataFrame(rng.normal(size=(n, g)) + X.iloc[:, :1].to_numpy(), index=ids,
                 columns=[f"g{j}" for j in range(g)])
CV = pd.Series(pd.factorize(groups)[0] % 5, index=ids)
LOLO = lineage.map({l: i for i, l in enumerate(sorted(lineage.unique()))}).astype(int)
EXCL = ["L3", "L10", "L11", "L40"]
CFGS = {"B1": {"model": "B1", "grid": {"kappa": [0, 2, 20]}},
        "M1": {"model": "M1", "params": {"modalities": ["X"], "residual": "B0",
                                         "lam_grid": [1.0, 10.0, 100.0], "min_rows": 20}}}


@pytest.mark.parametrize("scheme", ["cv", "lolo"])
@pytest.mark.parametrize("key", ["B1", "M1"])
def test_excluded_lines_never_train_and_cannot_change_predictions(scheme, key):
    folds = CV if scheme == "cv" else LOLO
    cfg = dict(CFGS[key], inner_folds=3, inner_min_obs=5, exclude_from_training=EXCL,
               reference_baselines={"B0": {}, "B1": {"grid": {"kappa": [1, 10]}}})
    Xk = X if key == "M1" else None
    outs = []
    for pert in (False, True):
        Yx = Y.copy()
        if pert:
            Yx.loc[EXCL] = rng.normal(size=(len(EXCL), g)) * 50
            Yx.loc[EXCL[:2], "g0"] = np.nan
        run.FIT_LOG.clear()
        outs.append(run.run_cv(Yx, folds, groups, lineage, Xk, cfg, Y.columns, set(), seed=7,
                               log=lambda *a: None))
        for e in run.FIT_LOG:
            assert not set(EXCL) & set(e["ids"])
            assert not set(folds.index[folds == e["outer_fold"]]) & set(e["ids"])
    a, b = outs
    assert set(EXCL) <= set(a["pred"].index)
    pd.testing.assert_frame_equal(a["pred"], b["pred"])
    for r in a["ref_pred"]:
        pd.testing.assert_frame_equal(a["ref_pred"][r], b["ref_pred"][r])
    assert a["chosen"] == b["chosen"]
    for f, tr in a["train_ids"].items():  # what save_fold_model records
        assert not set(EXCL) & set(tr) and not set(folds.index[folds == f]) & set(tr)


def test_m2_cache_key_is_exact():
    m = LocalStageM2(m1_params={"modalities": ["E"], "residual": "B0"})
    grp = groups.to_numpy()
    k = m._cache_key(Y, lineage, grp)
    assert k == m._cache_key(Y.copy(), lineage.copy(), grp.copy())
    Y2 = Y.copy()
    Y2.iloc[0, 0] += 1e-9
    variants = [m._cache_key(Y2, lineage, grp), m._cache_key(Y.iloc[::-1], lineage, grp[::-1]),
                m._cache_key(Y[Y.columns[::-1]], lineage, grp),
                m._cache_key(Y, lineage, np.roll(grp, 1)),
                m._cache_key(Y, lineage.sample(frac=1, random_state=0).set_axis(lineage.index), grp)]
    assert k not in variants and len(set(variants)) == len(variants)


def test_lolo_folds_marked_minus_one_are_never_held_out_but_always_train():
    folds = LOLO.copy()
    small = folds.index[folds == folds.max()]
    folds[small] = -1  # e.g. a lineage below lolo_min_lines
    cfg = dict(CFGS["B1"], inner_folds=3, inner_min_obs=5, reference_baselines={"B0": {}})
    run.FIT_LOG.clear()
    fl = sorted(int(f) for f in folds.unique() if f >= 0)
    out = run.run_cv(Y, folds, groups, lineage, None, cfg, Y.columns, set(), seed=7, folds_to_run=fl,
                     log=lambda *a: None)
    assert not set(small) & set(out["pred"].index)
    assert all(set(small) <= set(tr) for tr in out["train_ids"].values())
