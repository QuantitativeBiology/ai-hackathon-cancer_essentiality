"""scripts/predict.py end to end on a toy model (no real data): input parsing, harmonisation, applicability
domain, confidence / tier columns, persisted-model reuse (D22)."""
import json

import numpy as np
import pandas as pd
import pytest

import predict as PD
from essml.models.kernel_ridge import KernelRidgeM1


def _toy(seed=0, n=120, p=30, g=12):
    rng = np.random.default_rng(seed)
    ids = [f"ACH-{i:06d}" for i in range(n)]
    X = pd.DataFrame(rng.normal(5, 1, (n, p)), index=ids, columns=[str(1000 + j) for j in range(p)])
    B = np.zeros((p, g))
    B[:3, :6] = rng.normal(0, 0.4, (3, 6))
    Y = pd.DataFrame((X.to_numpy() - 5) @ B + rng.normal(0, 0.2, (n, g)) - 0.3,
                     index=ids, columns=[str(9000 + k) for k in range(g)])
    lin = pd.Series(np.where(np.arange(n) % 2, "L1", "L2"), index=ids)
    return X, Y, lin


def _table(Y, path):
    genes = list(Y.columns)
    t = pd.DataFrame({"entrez": genes, "symbol": [f"G{g}" for g in genes], "mean_r": 0.4,
                      "confidence": np.linspace(-0.1, 0.6, len(genes)),
                      "tier": (["high"] * 6 + ["low"] * 5 + ["not_evaluable"])[:len(genes)],
                      "source_run_id": "toy", "threshold_high_mean_r": 0.3})
    t.to_csv(path, index=False)
    return path


def _model(X, Y, train):
    m = KernelRidgeM1(residual="B0")
    return m.fit(Y.loc[train], X=X.loc[train], groups=np.array(train))


def test_read_expression_depmap_headers_and_rows(tmp_path):
    df = pd.DataFrame([[1.0, 2.0], [3.0, 4.0]], index=["s1", "s2"], columns=["TSPAN6 (7105)", "TNMD (64102)"])
    df.to_csv(tmp_path / "e.csv")
    x = PD.read_expression(tmp_path / "e.csv")
    assert list(x.columns) == ["7105", "64102"] and list(x.index) == ["s1", "s2"]
    df.T.to_csv(tmp_path / "t.tsv", sep="\t")
    pd.testing.assert_frame_equal(PD.read_expression(tmp_path / "t.tsv", genes_in_rows=True), x)
    pd.DataFrame([[1.0]], index=["s"], columns=["TSPAN6"]).to_csv(tmp_path / "bad.csv")
    with pytest.raises(ValueError, match="Entrez"):
        PD.read_expression(tmp_path / "bad.csv")


def test_predict_table_heldout_line_matches_model_and_carries_confidence(tmp_path):
    X, Y, lin = _toy()
    held = X.index[0]
    train = [i for i in X.index if i != held]
    m = _model(X, Y, train)
    table = PD.PR.load_table(_table(Y, tmp_path / "pg.csv"))
    ad = PD.applicability_domain(X.loc[[held]], X.loc[train], lin, k=5, n_hvg=20)
    out = PD.predict_table(X.loc[[held]], [(0, m, {})], table, x_builder=lambda E: E, ad=ad)
    assert list(out.columns) == ["sample", "entrez", "symbol", "predicted_effect", "confidence", "tier",
                                 "ad_distance", "ad_percentile_vs_dev", "nearest_lineage"]
    direct = m.predict([held], X=X.loc[[held]]).iloc[0]
    assert np.allclose(out.set_index("entrez")["predicted_effect"].to_numpy(), direct.to_numpy())
    assert out.set_index("entrez").loc[Y.columns[-1], "tier"] == "not_evaluable"
    assert 0 <= out["ad_percentile_vs_dev"].iloc[0] <= 100


def test_harmonise_quantile_maps_onto_reference():
    X, _, _ = _toy()
    shifted = X.iloc[:2] * 0.5 + 3
    h = PD.harmonise(shifted, X.iloc[2:], "quantile")
    ref = np.sort(X.iloc[2:].to_numpy(), axis=1).mean(0)
    assert np.allclose(np.sort(h.iloc[0].to_numpy()), ref)
    assert (h.iloc[0].rank() == shifted.iloc[0].rank()).all()   # order kept
    assert PD.harmonise(shifted, X, "none") is shifted


def test_main_end_to_end_on_toy_model(tmp_path, monkeypatch):
    X, Y, lin = _toy()
    held = X.index[0]
    train = [i for i in X.index if i != held]
    calls = {"fit": 0}

    def fake_get_models(name, train_ids, save=False, **kw):
        assert held not in set(train_ids)                     # the held-out line never trains
        calls["fit"] += 1
        return [(0, _model(X, Y, list(train_ids)), {})], "toy_tag"

    monkeypatch.setattr(PD, "dev_reference", lambda: (X, lin))
    monkeypatch.setattr(PD, "get_models", fake_get_models)
    monkeypatch.setattr(PD, "X_BUILDER", lambda E: E)        # toy model takes a plain feature frame
    monkeypatch.setattr(PD, "AD_K", 5)
    monkeypatch.setattr(PD, "N_HVG", 20)
    monkeypatch.setitem(PD.SUPPORTED, "m1_E_resB0", {"config": None, "table": _table(Y, tmp_path / "pg.csv")})
    X.loc[[held]].rename(columns=lambda c: f"G{c} ({c})").to_csv(tmp_path / "in.csv")
    out = tmp_path / "o" / "pred.csv"
    long = PD.main(["--expression", str(tmp_path / "in.csv"), "--out", str(out), "--exclude-lines", held])
    res = pd.read_csv(out, dtype={"entrez": str})
    assert calls["fit"] == 1 and len(res) == Y.shape[1] == len(long)
    assert {"sample", "predicted_effect", "confidence", "tier", "ad_distance", "ad_percentile_vs_dev"} <= set(res.columns)
    direct = _model(X, Y, train).predict([held], X=X.loc[[held]]).iloc[0]
    assert np.allclose(res["predicted_effect"].to_numpy(), direct.to_numpy())
    meta = json.loads(out.with_suffix(".meta.json").read_text())
    assert meta["excluded_lines"] == [held] and meta["n_training_lines"] == len(train)
    assert meta["predictable_genes_sha256"] == PD.PR.sha256_file(tmp_path / "pg.csv")


def test_m5_avg_not_supported_yet(tmp_path):
    with pytest.raises(SystemExit, match="m5_avg"):
        PD.main(["--expression", "x.csv", "--out", str(tmp_path / "o.csv"), "--model", "m5_avg"])
