"""M5 (m5_avg) on new lines from supplied feature blocks (G6 external profiles, scripts/predict.py): toy data.

Input refits are logged and training-only, kept inputs are reused, M2's expression fill uses training means
for every gene (never 0 / NaN), supplied genomic blocks fill unlisted entries with 0 = diploid wild type
(DESIGN §2), and restored targets must keep the training order.
"""
import numpy as np
import pandas as pd
import pytest

from essml import priors as PR
from essml.models import g3_common, m2_local
from essml.models.m5_stack import StackM5
from test_m5_stack import Ridge, _stack, _toy


def test_input_refits_are_logged_training_only_and_kept(monkeypatch):
    X, Y, lin, groups = _toy()
    tr, te = Y.index[:180], Y.index[180:]
    n_fit = {"a": 0, "b": 0}

    def counted(name):
        class C(Ridge):
            def fit(self, *a, **k):
                n_fit[name] += 1
                return super().fit(*a, **k)
        return lambda: C(X)

    monkeypatch.setattr(StackM5, "KEEP_INPUTS", True)
    st = _stack(X, {"a": counted("a"), "b": counted("b")}, mode="average", cache=False)
    st.fit(Y.loc[tr], lineage=lin, groups=groups[:180])
    p1 = st.predict(te[:10])
    p2 = st.predict(te[10:])
    assert n_fit == {"a": 1, "b": 1}                                   # kept inputs reused for later lines
    assert [e["input"] for e in st.input_fit_log_] == ["a", "b"]
    assert all(set(e["ids"]) == set(tr) for e in st.input_fit_log_)  # training lines only
    assert np.isfinite(p1.to_numpy()).all() and np.isfinite(p2.to_numpy()).all()
    with pytest.raises(ValueError, match="outside its training fold"):
        st.predict(tr[:2])


def test_restore_targets_requires_training_order():
    X, Y, lin, groups = _toy()
    st = _stack(X, {"a": lambda: Ridge(X)}, mode="average").fit(Y.iloc[:50], lineage=lin, groups=groups[:50])
    st.restore_targets(Y.iloc[:50])
    with pytest.raises(ValueError, match="order"):
        st.restore_targets(Y.iloc[:50].iloc[::-1])


def test_m2_expression_fill_uses_training_means_for_every_gene(monkeypatch):
    genes = ["1", "2", "3"]
    E = pd.DataFrame({"1": [1.0, 3.0, 9.0], "2": [5.0, 5.0, 0.0], "3": [2.0, 4.0, 100.0]},
                     index=["T1", "T2", "OTHER"])  # 'OTHER' is not a training line: never used for the mean
    top = pd.DataFrame({"par1": [np.nan] * 3, "id1": [np.nan] * 3}, index=genes)
    builder = PR.LocalFeatureBuilder(E, pd.DataFrame(), pd.DataFrame(), top,
                                     pd.DataFrame({"a": [], "b": [], "source": []}))
    monkeypatch.setattr(m2_local, "_BUILDER", builder)
    st = StackM5([{"name": "x", "model": "M1"}], mode="average")
    st.Y_ = pd.DataFrame(np.zeros((2, 3)), index=["T1", "T2"], columns=genes)
    new = pd.DataFrame({"1": [np.nan], "2": [7.0]}, index=["NEW"])   # gene 1 NaN, gene 3 absent
    b = st._local_builder({"E": new}, None)
    row = b.E.loc["NEW"]
    assert row["1"] == pytest.approx(2.0) and row["2"] == 7.0 and row["3"] == pytest.approx(3.0)
    f = b.local_features(["NEW"], genes)
    assert np.isfinite(f["own_expr"]).all() and f["own_cn"][0, 0] == 0.0 and f["own_mut"][0, 0] == 0.0


def test_supplied_genomic_block_fills_unlisted_entries_with_diploid_wt(monkeypatch):
    rng = np.random.default_rng(0)
    ids = [f"Z{i:03d}" for i in range(40)]
    E = pd.DataFrame(rng.normal(size=(40, 5)), index=ids, columns=[f"e{j}" for j in range(5)])
    C = pd.DataFrame(rng.normal(0.5, 0.3, size=(40, 4)), index=ids, columns=[f"c{j}" for j in range(4)])
    monkeypatch.setattr(g3_common, "load_blocks", lambda names, *a, **k: {n: {"E": E, "C": C}[n] for n in names})
    cf, _ = g3_common.CellFeatures(("E", "C")).fit(ids)
    newE = E.iloc[:1].rename(index=lambda _: "NEW")
    partial = pd.DataFrame({"c0": [0.7], "c1": [np.nan]}, index=["NEW"])        # c2, c3 not listed, c1 NaN
    full = pd.DataFrame({"c0": [0.7], "c1": [0.0], "c2": [0.0], "c3": [0.0]}, index=["NEW"])
    a = cf.transform(["NEW"], {"E": newE, "C": partial})
    b = cf.transform(["NEW"], {"E": newE, "C": full})
    assert np.allclose(a, b)
    absent = cf.transform(["NEW"], {"E": newE})                                 # block absent: 0 + absent flag
    cols = cf.columns_
    assert absent[0, cols.index("C:absent")] == 1.0 and a[0, cols.index("C:absent")] == 0.0
