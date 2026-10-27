"""M5 stacking (DESIGN §5) on toy data: non-negativity, nesting, permuted-input control."""
import numpy as np
import pandas as pd
import pytest

from essml.models.m5_stack import StackM5, class_weights, gram, nnls_gram

rng = np.random.default_rng(5)


def _toy(n=240, p=8, G=30):
    ids = pd.Index([f"L{i:03d}" for i in range(n)])
    X = pd.DataFrame(rng.normal(size=(n, p)), index=ids)
    B = rng.normal(size=(p, G)) * 0.5
    Y = pd.DataFrame(X.to_numpy() @ B + rng.normal(size=(n, G)) * 0.7, index=ids,
                     columns=[str(1000 + g) for g in range(G)])
    lin = pd.Series(rng.choice(list("abc"), size=n), index=ids)
    groups = np.array([f"P{i // 2}" for i in range(n)])  # patients with 2 lines each
    return X, Y, lin, groups


class Ridge:
    """Least squares on the toy features (+ optional noise on the coefficients / sign flip)."""

    def __init__(self, X, noise=0.0, sign=1.0, seed=0):
        self.X, self.noise, self.sign, self.seed = X, noise, sign, seed

    def fit(self, Y, lineage=None, X=None, groups=None):
        Z = np.column_stack([np.ones(len(Y)), self.X.loc[Y.index].to_numpy()])
        self.B_ = np.linalg.lstsq(Z, Y.to_numpy(), rcond=None)[0]
        self.B_[1:] += np.random.default_rng(self.seed).normal(size=self.B_[1:].shape) * self.noise
        self.cols_ = Y.columns
        return self

    def predict(self, index, lineage=None, X=None):
        Z = np.column_stack([np.ones(len(index)), self.X.loc[index].to_numpy()])
        P = Z @ self.B_
        P = P.mean(0) + self.sign * (P - P.mean(0))
        return pd.DataFrame(P, index=index, columns=self.cols_)


class Memoriser:
    """Returns the training target for lines it was fitted on (perfect in-sample), noise else."""

    def fit(self, Y, lineage=None, X=None, groups=None):
        self.Y_ = Y
        return self

    def predict(self, index, lineage=None, X=None):
        mu = self.Y_.mean()
        P = pd.DataFrame(np.tile(mu.to_numpy(), (len(index), 1)), index=index, columns=self.Y_.columns)
        P += np.random.default_rng(len(index)).normal(size=P.shape) * self.Y_.std().to_numpy()
        seen = index.intersection(self.Y_.index)
        P.loc[seen] = self.Y_.loc[seen]
        return P


def _stack(X, models, **kw):
    specs = [{"name": n, "model": "toy"} for n in models]
    return StackM5(specs, input_factory=lambda spec, inner: models[spec["name"]](), **kw)


def test_nnls_gram_is_nonnegative_and_exact_when_unconstrained():
    A = np.array([[2.0, 0.5], [0.5, 1.0]])
    w_true = np.array([0.3, 0.7])
    assert np.allclose(nnls_gram(A, A @ w_true), w_true)
    w = nnls_gram(A, A @ np.array([1.0, -0.5]))  # unconstrained optimum has a negative weight
    assert (w >= 0).all() and w[1] == 0


def test_weights_nonnegative_sum_to_one_and_reject_anticorrelated_input():
    X, Y, lin, groups = _toy()
    st = _stack(X, {"good": lambda: Ridge(X), "flipped": lambda: Ridge(X, sign=-1.0),
                    "noisy": lambda: Ridge(X, noise=0.3, seed=1)})
    st.fit(Y, lineage=lin, groups=groups)
    for c, w in st.weights_.items():
        assert (w >= 0).all() and np.isclose(w.sum(), 1.0), c
    k = st.names_.index("flipped")
    assert all(w[k] == 0 for w in st.weights_.values())
    k = st.names_.index("good")
    assert all(w[k] > 0.5 for w in st.weights_.values())


def test_nesting_weights_never_see_own_fold_predictions():
    """A memoriser is perfect on lines it was fitted on. With nested OOF predictions it is pure
    noise, so it must get ~0 weight; a non-nested stacker would give it weight ~1."""
    X, Y, lin, groups = _toy()
    st = _stack(X, {"ridge": lambda: Ridge(X), "memo": Memoriser})
    st.fit(Y, lineage=lin, groups=groups)
    st.check_nesting()
    for e in st.oof_log_:
        assert not (e["fit"] & e["pred"])
        # patients are never split between the fit and the predicted lines (grouped folds)
        pid = dict(zip(Y.index, groups))
        assert not ({pid[i] for i in e["fit"]} & {pid[i] for i in e["pred"]})
    k = st.names_.index("memo")
    assert all(w[k] < 0.02 for w in st.weights_.values())
    # sanity: in-sample (non-nested) predictions would have favoured the memoriser
    full = {n: f().fit(Y).predict(Y.index).to_numpy() for n, f in
            {"ridge": lambda: Ridge(X), "memo": Memoriser}.items()}
    w_leaky, _ = class_weights([full["ridge"], full["memo"]], Y.to_numpy(), st.labels_)
    assert all(w[1] > 0.9 for w in w_leaky.values())


def test_oof_rows_do_not_depend_on_their_own_targets():
    X, Y, lin, groups = _toy()
    mk = {"ridge": lambda: Ridge(X)}
    st = _stack(X, mk).fit(Y, lineage=lin, groups=groups)
    splits = st._splitter()._inner_splits(Y.index, lin, groups)
    b0 = splits[0][1]
    Y2 = Y.copy()
    Y2.iloc[b0] = rng.normal(size=(len(b0), Y.shape[1])) * 10
    st2 = _stack(X, mk).fit(Y2, lineage=lin, groups=groups)
    assert np.allclose(st.oof_["ridge"][b0], st2.oof_["ridge"][b0])
    assert not np.allclose(st.oof_["ridge"][splits[1][1]], st2.oof_["ridge"][splits[1][1]])


def test_permuted_input_control_gets_about_zero_weight():
    X, Y, lin, groups = _toy()
    mk = {"a": lambda: Ridge(X, noise=0.1, seed=2), "b": lambda: Ridge(X, noise=0.1, seed=3)}
    st = _stack(X, mk).fit(Y, lineage=lin, groups=groups)
    kb = st.names_.index("b")
    assert all(w[kb] > 0.2 for w in st.weights_.values())  # unpermuted: both inputs are used
    ctl = _stack(X, mk, permute_input="b").fit(Y, lineage=lin, groups=groups)
    assert all(w[kb] < 0.02 for w in ctl.weights_.values())
    tr, te = Y.index[:180], Y.index[180:]
    ctl = _stack(X, mk, permute_input="b").fit(Y.loc[tr], lineage=lin, groups=groups[:180])
    P = ctl.predict(te)
    r = np.mean([np.corrcoef(P[g], Y.loc[te, g])[0, 1] for g in Y.columns])
    assert r > 0.5


def test_average_mode_and_prediction_is_weighted_sum():
    X, Y, lin, groups = _toy()
    tr, te = Y.index[:180], Y.index[180:]
    mk = {"a": lambda: Ridge(X), "b": lambda: Ridge(X, noise=0.3, seed=1)}
    avg = _stack(X, mk, mode="average").fit(Y.loc[tr], lineage=lin, groups=groups[:180])
    st = _stack(X, mk).fit(Y.loc[tr], lineage=lin, groups=groups[:180])
    pa, pb = (mk[n]().fit(Y.loc[tr]).predict(te) for n in ("a", "b"))
    assert np.allclose(avg.predict(te), (pa + pb) / 2)
    P = st.predict(te).to_numpy()
    exp = np.zeros_like(P)
    for c, w in st.weights_.items():
        m = st.labels_ == c
        exp[:, m] = w[0] * pa.to_numpy()[:, m] + w[1] * pb.to_numpy()[:, m]
    assert np.allclose(P, exp)


def test_gram_centres_and_scales_per_gene():
    Yv = np.array([[1.0, 10.0], [2.0, np.nan], [3.0, 30.0]])
    A, b = gram([Yv.copy() + 5.0], Yv)  # prediction = truth + constant -> perfect after centring
    assert np.allclose(A[:, 0, 0], b[:, 0])
    assert np.allclose(A[:, 0, 0], [2.0, 1.0])  # (n_obs - 1) after scaling by 1/SD


def test_rejects_bad_config():
    with pytest.raises(ValueError):
        StackM5([{"name": "a", "model": "M1"}], permute_input="zzz")
    with pytest.raises(ValueError):
        StackM5([{"name": "a", "model": "M1"}, {"name": "a", "model": "M2"}])
