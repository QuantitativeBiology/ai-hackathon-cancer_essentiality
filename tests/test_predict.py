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
    long = PD.main(["--expression", str(tmp_path / "in.csv"), "--out", str(out), "--exclude-lines", held,
                    "--model", "m1_E_resB0"])
    res = pd.read_csv(out, dtype={"entrez": str})
    assert calls["fit"] == 1 and len(res) == Y.shape[1] == len(long)
    assert {"sample", "predicted_effect", "confidence", "tier", "ad_distance", "ad_percentile_vs_dev"} <= set(res.columns)
    direct = _model(X, Y, train).predict([held], X=X.loc[[held]]).iloc[0]
    assert np.allclose(res["predicted_effect"].to_numpy(), direct.to_numpy())
    assert "note" in res.columns and res["note"].fillna("").eq("").all()       # toy symbols: no flagged gene
    meta = json.loads(out.with_suffix(".meta.json").read_text())
    assert meta["excluded_lines"] == [held] and meta["n_training_lines"] == len(train)
    assert meta["predictable_genes_sha256"] == PD.PR.sha256_file(tmp_path / "pg.csv")


def test_m5_avg_is_the_default_model_with_the_primary_tiers():
    ap_default = PD.DEFAULT_MODEL
    assert ap_default == "m5_avg" and PD.SUPPORTED["m5_avg"]["table"].name == "predictable_genes.csv"
    assert PD.SUPPORTED["m1_E_resB0"]["table"].name == "predictable_genes_m1.csv"
    assert "diploid wild type" in PD.NO_GENOMICS_WARNING and "D24" in PD.NO_GENOMICS_WARNING
    assert "post hoc" in PD.NO_GENOMICS_WARNING and "0.3135" in PD.NO_GENOMICS_WARNING


def test_default_follows_the_pi_rule_on_the_expression_only_dev_cv_result():
    """PI rule fixed before the number was seen: m5_avg stays the default only if its expression-only dev-CV P1
    (variable) beats m1_E_resB0 by > delta = 0.005 and > 1 SE of the paired fold differences (DESIGN §4)."""
    p = PD.ROOT / "results" / "predict_expression_only_cv.json"
    res = json.loads(p.read_text())
    d = res["P1_by_class"]["variable"]["delta_P1_vs_m1_E_resB0"]
    met = d["mean"] > 0.005 and d["mean"] > d["se"]
    assert res["post_hoc"] is True and res["default_rule"]["met"] == met
    assert PD.DEFAULT_MODEL == ("m5_avg" if met else "m1_E_resB0") == res["default_rule"]["predict_py_default"]
    c = PD.EXPRESSION_ONLY_CV                                   # constants quoted in the warning match the result
    v = res["P1_by_class"]["variable"]
    assert round(v["P1_mean"], 4) == c["m5_avg_P1_variable"]
    assert round(v["m1_E_resB0_P1_mean"], 4) == c["m1_E_resB0_P1_variable"]
    assert round(v["m5_avg_with_genomics_P1_mean"], 4) == c["m5_avg_with_genomics_P1_variable"]
    assert round(d["mean"], 4) == c["dP1_vs_m1"] and round(d["se"], 4) == c["dP1_vs_m1_se"]
    nras = res["hotspot_auprc"]["NRAS_hotspot"]
    assert round(nras["m5_avg_expression_only"]["ap_all_mean"], 2) == c["nras_auprc_m5_avg_expression_only"]
    assert round(nras["m1_E_resB0"]["ap_all_mean"], 2) == c["nras_auprc_m1_E_resB0"]
    assert round(nras["m5_avg_with_genomics"]["ap_all_mean"], 2) == c["nras_auprc_m5_avg_with_genomics"]


def test_mutation_notes_per_model():
    long = pd.DataFrame({"symbol": ["NRAS", "KRAS", "BRAF", "PIK3CA", "SOX10"]})
    n = PD.add_notes(long, "m1_E_resB0")["note"].tolist()
    assert n[:4] == [PD.MUTATION_NOTE] * 4 and n[4] == ""
    assert PD.MUTATION_NOTE == "mutation-driven dependency; underpredicted by both expression-only configurations"
    for text in (PD.MUTATION_NOTE, PD.M1_WARNING, PD.NO_GENOMICS_NOTE, PD.NO_GENOMICS_WARNING):
        assert "use m5_avg" not in text and "--model m5_avg" not in text   # predict.py takes no genomics
    assert "both expression-only configurations" in PD.M1_WARNING
    n5 = PD.add_notes(long, "m5_avg")["note"].tolist()             # no genomic input supplied: flagged too
    assert n5[:4] == [PD.NO_GENOMICS_NOTE] * 4 and n5[4] == ""
    assert (PD.add_notes(long, "m5_avg", genomics_given=True)["note"] == "").all()


def test_fit_model_refuses_non_dev_lines():
    """predict.py trains on DepMap dev lines only: a test-set ID is refused before any target is read."""
    sp = pd.read_csv(PD.ROOT / "configs" / "splits.csv", index_col=0)
    test_id = sp.index[sp["role"] == "test"][0]
    dev_id = sp.index[sp["role"] == "dev"][0]
    with pytest.raises(PD.guard.BlindingError):
        PD.fit_model("m1_E_resB0", [dev_id, test_id])


def test_run_id_selects_the_run_models_tier_table():
    assert PD.PR.table_for_run("20261026-222117_m1_E_resB0_91e52f9d").name == "predictable_genes_m1.csv"
    assert PD.PR.table_for_run("20261027-161041_m5_avg_6dbebcf4").name == "predictable_genes.csv"
    assert PD.PR.table_for_run("20261027-023409_m1_E_resB0_noEC_1d8386a5") is None


# ---------------------------------------------------------------- user-supplied genomics (post hoc, 2026-10-28)
def test_read_table_accepts_depmap_release_layout(tmp_path):
    df = pd.DataFrame({"SequencingID": ["CDS-1", "CDS-2", "CDS-3"], "ModelConditionID": ["MC-1", "MC-2", "MC-3"],
                       "ModelID": ["ACH-1", "ACH-1", "ACH-2"], "IsDefaultEntryForMC": ["Yes", "No", "Yes"],
                       "IsDefaultEntryForModel": ["Yes", "No", "Yes"],
                       "TSPAN6 (7105)": [1.0, 9.0, 2.0], "TNMD (64102)": [0.5, 9.0, 1.5]}, index=[0, 1, 2])
    df.to_csv(tmp_path / "rel.csv")
    x = PD.read_expression(tmp_path / "rel.csv")
    assert list(x.index) == ["ACH-1", "ACH-2"] and list(x.columns) == ["7105", "64102"]
    assert x.loc["ACH-1", "7105"] == 1.0                     # the default entry, not the other profile


def _fake_columns(monkeypatch):
    cols = {"cn": ["1", "2", "3"], "C": ["1", "2", "3"], "M": ["1", "2", "3"], "H": ["2", "3"]}
    monkeypatch.setattr(PD, "_cached_columns", lambda name: cols[name])


def test_genomic_blocks_follow_the_depmap_feature_convention(monkeypatch):
    _fake_columns(monkeypatch)
    from essml import features as F
    cn = pd.DataFrame({"1": [1.0, 4.0], "2": [np.nan, 0.5]}, index=["s1", "s2"])   # gene "3" not supplied
    dam = pd.DataFrame({"1": [2.0], "3": [1.0]}, index=["s1"])                     # s2 missing -> absent
    hot = pd.DataFrame({"2": [0.0, 2.0], "3": [1.0, 0.0], "9": [1.0, 1.0]}, index=["s1", "s2"])
    b = PD.genomic_blocks(["s1", "s2"], cn, dam, hot, None, log=lambda *_: None)
    ref = F.cn_transform(cn.reindex(columns=["1", "2"]))
    assert np.allclose(b["C"][["1", "2"]].to_numpy(), ref.to_numpy())
    assert b["C"].loc["s1", "1"] == 0.0 and b["C"].loc["s1", "2"] == 0.0 and (b["C"]["3"] == 0).all()
    assert np.isnan(b["_cn"].loc["s1", "2"]) and np.isnan(b["_cn"]["3"]).all()   # M2 fills as for DepMap NaN
    assert list(b["M"].index) == ["s1"] and b["M"].loc["s1"].tolist() == [2.0, 0.0, 1.0]
    assert b["_damaging"] is b["M"]
    assert list(b["H"].columns) == ["2", "3"] and b["H"].to_numpy().tolist() == [[0.0, 1.0], [1.0, 0.0]]
    assert "meta" not in b
    ids = PD.with_input_ids(b)
    assert list(ids["C"].index) == ["input::s1", "input::s2"]


def test_read_metadata_matches_build_meta(tmp_path):
    from essml import features as F
    m = pd.DataFrame({"ModelID": ["s1", "s2"], "OncotreeLineage": ["Skin", "Unknown"],
                      "OncotreePrimaryDisease": ["Melanoma", ""], "Sex": ["Female", "Male"], "Age": ["61", ""],
                      "CCLEName": ["X", "Y"]})
    m.to_csv(tmp_path / "meta.csv", index=False)
    got = PD.read_metadata(tmp_path / "meta.csv")
    ref = F.build_meta(m.astype(str).replace({"": pd.NA}), F.META_COLS)
    assert list(got.index) == ["s1", "s2"] and got.loc["s1", "OncotreeLineage"] == "Skin"
    assert pd.isna(got.loc["s2", "OncotreeLineage"]) and pd.isna(got.loc["s2", "Age"])
    assert got.loc["s1", "Age"] == ref.loc["s1", "Age"] == 61


def test_genomic_options_refused_for_the_expression_only_model(tmp_path):
    with pytest.raises(SystemExit):
        PD.main(["--expression", "x.csv", "--out", str(tmp_path / "o.csv"), "--model", "m1_E_resB0",
                 "--cn", "cn.csv"])


def test_hotspot_input_drops_the_no_genomics_note():
    long = pd.DataFrame({"symbol": ["NRAS", "SOX10"]})
    assert PD.add_notes(long, "m5_avg", genomics_given=True)["note"].tolist() == ["", ""]
    assert PD.add_notes(long, "m5_avg", genomics_given=False)["note"].tolist() == [PD.NO_GENOMICS_NOTE, ""]
