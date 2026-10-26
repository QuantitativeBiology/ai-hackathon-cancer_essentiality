"""Launcher leakage tests ("probe bugs become tests", DESIGN §8).

The flawed fold-0 probe tuned lambda on the validation fold and filled NaNs with validation lines.
Here: (1) every fit (inner tuning and outer refit) is instrumented and must never receive an ID of
the held-out outer fold; (2) perturbing held-out targets (values and NaN pattern) must leave
held-out predictions and chosen hyper-parameters bit-identical; (3) the sentinel raises.
"""
import numpy as np
import pandas as pd
import pytest

import run
from essml import modelstore
from essml.baselines import GeneMean
from essml.evaluate import evaluate
from essml.genesets import gene_classes

rng = np.random.default_rng(3)


def _toy(n=150, g=30, p=40):
    ids = pd.Index([f"L{i}" for i in range(n)])
    lineage = pd.Series(rng.choice(list("abcde"), size=n), index=ids)
    groups = pd.Series([f"P{i // 2}" for i in range(n)], index=ids)  # patient pairs
    folds = pd.Series(pd.factorize(groups)[0] % 5, index=ids)
    X = pd.DataFrame(rng.normal(size=(n, p)), index=ids, columns=[f"x{j}" for j in range(p)])
    eff = {l: rng.normal(size=g) for l in "abcde"}
    Y = pd.DataFrame(np.vstack([eff[l] for l in lineage]) * 0.5 + rng.normal(size=(n, g)) * 0.5,
                     index=ids, columns=[f"g{j}" for j in range(g)])
    Y.iloc[:, :3] += X.iloc[:, :3].to_numpy()  # expression-driven genes for B2
    Y = Y.mask(rng.random(Y.shape) < 0.05)
    return Y, folds, groups, lineage, X


CFGS = {
    "B1": {"model": "B1", "grid": {"kappa": [0, 2, 20]}, "inner_folds": 3, "inner_min_obs": 5,
           "reference_baselines": {"B0": {}}},
    "B2": {"model": "B2", "params": {"n_top": 20}, "grid": {"k": [3, 10]}, "inner_folds": 3,
           "inner_min_obs": 5, "reference_baselines": {"B0": {}, "B1": {"grid": {"kappa": [1, 10]}}}},
}


def _run(cfg, Y, folds, groups, lineage, X):
    return run.run_cv(Y, folds, groups, lineage, X if cfg["model"] == "B2" else None, cfg,
                      Y.columns, set(), seed=7, log=lambda *a: None)


@pytest.mark.parametrize("key", ["B1", "B2"])
def test_no_fit_sees_held_out_ids(key):
    Y, folds, groups, lineage, X = _toy()
    run.FIT_LOG.clear()
    _run(CFGS[key], Y, folds, groups, lineage, X)
    assert len(run.FIT_LOG) > 5 * 2
    tags = {e["tag"] for e in run.FIT_LOG}
    assert "outer" in tags and any(t.startswith("inner") for t in tags)
    for e in run.FIT_LOG:
        held = set(folds.index[folds == e["outer_fold"]])
        assert not held.intersection(e["ids"]), e["tag"]


@pytest.mark.parametrize("key", ["B1", "B2"])
def test_held_out_targets_do_not_affect_predictions(key):
    Y, folds, groups, lineage, X = _toy()
    a = _run(CFGS[key], Y, folds, groups, lineage, X)
    Y2 = Y.copy()
    ho = folds.index[folds == 0]
    Y2.loc[ho] = rng.normal(size=(len(ho), Y.shape[1])) * 10
    Y2.loc[ho[:10], "g0"] = np.nan
    b = _run(CFGS[key], Y2, folds, groups, lineage, X)
    assert a["chosen"][0] == b["chosen"][0]
    pd.testing.assert_frame_equal(a["pred"].loc[ho], b["pred"].loc[ho])
    for ref in a["ref_pred"]:
        pd.testing.assert_frame_equal(a["ref_pred"][ref].loc[ho], b["ref_pred"][ref].loc[ho])


def test_sentinel_raises():
    Y, folds, _, lineage, _ = _toy()
    s = run.Sentinel(folds.index[folds == 0], 0)
    with pytest.raises(run.LeakageError):
        run.guarded_fit(GeneMean(), s, "outer", Y, lineage, None)


def test_planted_lineage_signal_end_to_end():
    Y, folds, groups, lineage, X = _toy()
    cfg = dict(CFGS["B1"], reference_baselines={"B0": {}})
    res = _run(cfg, Y, folds, groups, lineage, X)
    cls = gene_classes(Y, neg_entrez=set())  # one dev-wide class table, as in run.py
    ev = evaluate(res["pred"], Y, cls, folds, baseline_preds=res["ref_pred"], min_obs=10)
    assert ev["summary"]["all"]["P1"] > 0.4
    evb0 = evaluate(res["ref_pred"]["B0"], Y, cls, folds, min_obs=10)
    assert evb0["summary"]["all"]["P1"] == 0.0


def test_modelstore_predict_new_roundtrip(tmp_path):
    Y, folds, groups, lineage, X = _toy()
    res = _run(CFGS["B2"], Y, folds, groups, lineage, X)
    for f, m in res["models"].items():
        modelstore.save_fold_model("toy", f, m, folds.index[folds != f], X.columns, model_dir=tmp_path)
    new = pd.Index(["L0", "L1", "L2"])
    exp = sum(m.predict(new, X=X) for m in res["models"].values()) / 5
    got = modelstore.predict_new("toy", new, X=X, model_dir=tmp_path,
                                 target_loader=lambda ids: Y.loc[ids])
    pd.testing.assert_frame_equal(got, exp)
    # anchor = dev line: only the fold model that did not train on it (its OOF model)
    pa, used = modelstore.predict_anchor("toy", "L0", X=X, model_dir=tmp_path,
                                         target_loader=lambda ids: Y.loc[ids])
    assert used == [folds["L0"]]
    pd.testing.assert_frame_equal(pa, res["pred"].loc[["L0"]])


def test_registry_row_replaced(tmp_path):
    p = tmp_path / "reg.csv"
    run.write_registry({"run_id": "r1", "status": "running"}, p)
    run.write_registry({"run_id": "r2", "status": "running"}, p)
    run.write_registry({"run_id": "r1", "status": "ok", "P1_all": 0.1}, p)
    reg = pd.read_csv(p)
    assert list(reg["run_id"]) == ["r2", "r1"] and list(reg["status"]) == ["running", "ok"]
    assert list(reg.columns) == run.REG_COLS
