"""M2 on the D14 base (G3): nested local-lambda CV (review finding 15) and pass-through of M1 for
non-variable genes. Toy data, frame mode, fast."""
import numpy as np
import pandas as pd

from essml.models.m2_local import LocalStageM2
from essml.priors import LocalFeatureBuilder

rng = np.random.default_rng(21)


def _toy(n=150, G=10):
    ids = pd.Index([f"L{i}" for i in range(n)])
    genes = [str(200 + j) for j in range(G)]
    E = pd.DataFrame(rng.normal(3, 1, size=(n, G)), index=ids, columns=genes)
    cn = pd.DataFrame(1.0, index=ids, columns=genes)
    dam = pd.DataFrame(0.0, index=ids, columns=genes)
    top = pd.DataFrame({"par1": [genes[(j + 1) % G] for j in range(G)], "par2": np.nan, "par3": np.nan,
                        "par1_score": 50.0, "par2_score": np.nan, "par3_score": np.nan,
                        "family_size": 1, "family_size_id20": 1, "max_identity": 50.0},
                       index=pd.Index(genes, name="entrez"))
    partners = pd.DataFrame({"a": genes, "b": genes[1:] + genes[:1], "source": "union"})
    builder = LocalFeatureBuilder(E, cn, dam, top, partners)
    Xm1 = pd.DataFrame(rng.normal(size=(n, 6)), index=ids, columns=[f"x{j}" for j in range(6)])
    Y = pd.DataFrame(-0.8 * E.to_numpy() + 0.3 * rng.normal(size=(n, G)), index=ids, columns=genes)
    Y.iloc[:, G // 2:] *= 0.05  # low-variance genes (SD < 0.25): not 'variable'
    lin = pd.Series(rng.choice(list("ab"), size=n), index=ids)
    return ids, genes, builder, Xm1, Y, lin


def _m2(**kw):
    return LocalStageM2(m1_params={"modalities": ["X"], "lam_grid": [10.0, 100.0], "min_rows": 20,
                                   "residual": "B0"},
                        priors="none", inner_folds=3, min_rows=20, cache_inner=False, **kw)


def test_nested_lambda_features_never_see_validation_rows():
    ids, genes, builder, Xm1, Y, lin = _toy()
    groups = np.arange(len(ids))
    m = _m2(nested_lambda=True)
    m.genes_ = Y.columns
    splits = m._inner_splits(ids, lin, groups)
    a0, b0 = splits[0]
    Y2 = Y.copy()
    Y2.iloc[b0] = rng.normal(size=(len(b0), Y.shape[1])) * 5  # validation targets of split 0
    X2 = Xm1.copy()
    X2.iloc[b0] = rng.normal(size=(len(b0), Xm1.shape[1])) * 5  # and their M1 features
    feats = {}
    for tag, (Yx, Xx) in {"orig": (Y, Xm1), "pert": (Y2, X2)}.items():
        oof, _, _ = m.inner_oof(Yx, lin, groups, Xx)
        for nested in (True, False):
            m.nested_lambda = nested
            feats[tag, nested] = m.lambda_cv_train_features(Yx, lin, groups, Xx, splits, oof)[0]
    # nested: the training-row features of split 0 do not depend on split 0's validation rows
    # (equal up to BLAS rounding of differently laid-out copies: ~1e-15)
    assert np.allclose(feats["orig", True][0], feats["pert", True][0], rtol=0, atol=1e-10)
    # G2 (plain OOF): they do (the optimism the review flagged)
    assert np.abs(feats["orig", False][0] - feats["pert", False][0]).max() > 0.1


def test_variable_only_local_stage_passes_m1_through_and_learns_signal():
    ids, genes, builder, Xm1, Y, lin = _toy()
    tr, te = ids[:110], ids[110:]
    m = _m2(nested_lambda=True, local_genes="variable")
    m.fit(Y.loc[tr], lineage=lin, X={"_m1": Xm1.loc[tr], "_local": builder}, groups=np.arange(110))
    P = m.predict(te, lineage=lin, X={"_m1": Xm1.loc[te], "_local": builder})
    var = Y.loc[tr].std().to_numpy() > 0.25
    assert var[: len(genes) // 2].all() and not var[len(genes) // 2:].any()
    assert (m.ok_ == var).all() and list(m.lam_by_group_) == [1]
    m1 = m.m1_.predict(te, lineage=lin, X=Xm1.loc[te])
    pd.testing.assert_frame_equal(P.loc[:, ~var], m1.loc[:, ~var])  # exact pass-through
    r = [np.corrcoef(P[g], Y.loc[te, g])[0, 1] for g in np.array(genes)[var]]
    assert np.mean(r) > 0.8
