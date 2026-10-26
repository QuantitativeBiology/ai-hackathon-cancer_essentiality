"""evaluate(): per-fold averaging vs pooling; B0 per-fold r = 0; class handling (DESIGN §4)."""
import numpy as np
import pandas as pd

from essml.evaluate import evaluate
from essml.genesets import gene_classes

rng = np.random.default_rng(2)


def _toy(n_per=60, k=5, g=40):
    ids = [f"L{i}" for i in range(n_per * k)]
    fold = pd.Series(np.repeat(np.arange(k), n_per), index=ids)
    Y = pd.DataFrame(rng.normal(size=(len(ids), g)) * 0.5 - 0.3, index=ids,
                     columns=[f"g{j}" for j in range(g)])
    Y = Y.mask(rng.random(Y.shape) < 0.05)
    return Y, fold


def _b0_oof(Y, fold):
    P = pd.DataFrame(index=Y.index, columns=Y.columns, dtype=float)
    for f in fold.unique():
        tr = fold.index[fold != f]
        P.loc[fold.index[fold == f]] = np.tile(Y.loc[tr].mean().to_numpy(), ((fold == f).sum(), 1))
    return P


def _cls(Y, fold, neg=()):
    return {f: gene_classes(Y.loc[fold != f], neg_entrez=set(neg)) for f in fold.unique()}


def test_b0_per_fold_exact_zero_and_pooling_artefact():
    Y, fold = _toy()
    P = _b0_oof(Y, fold)
    ev = evaluate(P, Y, _cls(Y, fold), fold, min_obs=20)
    assert ev["summary"]["all"]["P1"] == 0.0
    assert all(v == 0.0 for v in ev["summary"]["all"]["P1_folds"])
    assert (ev["per_gene"]["r"] == 0).all()
    # Pooling the same OOF predictions: each line's prediction is the mean of the OTHER folds,
    # which anti-correlates with its own fold's mean -> spurious negative r for B0.
    pooled = evaluate(P, Y, gene_classes(Y, neg_entrez=set()), None, min_obs=20)
    assert pooled["summary"]["all"]["P1"] < -0.1


def test_perfect_prediction_and_baseline_deltas():
    Y, fold = _toy()
    P = Y.fillna(0.0)
    B0 = _b0_oof(Y, fold)
    ev = evaluate(P, Y, _cls(Y, fold), fold, baseline_preds={"B0": B0}, min_obs=20)
    s = ev["summary"]["all"]
    assert np.isclose(s["P1"], 1.0) and np.isclose(s["dP1_vs_B0"], 1.0) and s["frac_beat_B0"] == 1.0
    assert np.isclose(s["P2a"], 1.0) and s["dP2a_vs_B0"] > 0
    assert np.isclose(s["P2b_nap"], 1.0) and np.isclose(s["P2b_auroc"], 1.0)
    assert np.isclose(ev["headline"]["P1_all"], 1.0) and s["P1_sd"] < 1e-12
    assert len(s["P1_folds"]) == 5


def test_fixed_gene_list_and_class_intersection():
    Y, fold = _toy()
    P = Y.fillna(0.0) + rng.normal(size=Y.shape)
    ev = evaluate(P, Y, _cls(Y, fold, neg={"g0", "g1"}), fold, genes=["g0", "g5", "g6"], min_obs=20)
    pf = ev["per_fold"]
    assert pf.loc[pf["class"] == "all", "n_genes"].eq(3).all()
    assert pf.loc[pf["class"] == "nonessential", "n_genes"].eq(1).all()
    assert set(ev["per_gene"]["gene"]) == {"g0", "g5", "g6"}


def test_min_obs_excludes_gene_in_fold():
    Y, fold = _toy()
    Y.loc[fold.index[fold == 0][:50], "g3"] = np.nan  # fold 0 keeps <= 10 observations of g3
    P = Y.fillna(0.0)
    ev = evaluate(P, Y, _cls(Y, fold), fold, min_obs=20)
    pg = ev["per_gene"].set_index(["fold", "gene"])
    assert np.isnan(pg.loc[(0, "g3"), "r"]) and np.isfinite(pg.loc[(1, "g3"), "r"])
    pf = ev["per_fold"]
    assert pf.loc[(pf["class"] == "all") & (pf["fold"] == 0), "n_genes"].item() == 39
