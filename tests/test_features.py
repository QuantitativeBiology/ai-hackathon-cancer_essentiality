"""Toy-data tests with known answers for essml.features and essml.priors (DESIGN §8)."""
import gzip

import numpy as np
import pandas as pd
import pytest

from essml import features as F
from essml import priors as P


# ----------------------------------------------------------------------------- CN
def test_cn_transform_known_values():
    cn = pd.DataFrame({"1": [1.0, 0.5, 2.0], "2": [np.nan, 1.0, 0.0], "3": [np.nan] * 3,
                       "4": [100.0, 1.0, 1.0]}, index=["a", "b", "c"])
    c = F.cn_transform(cn)
    assert list(c.columns) == ["1", "2", "4"]               # all-NaN gene dropped
    assert c.loc["a", "1"] == 0.0 and c.loc["b", "2"] == 0.0  # CN = 1 -> 0 (diploid)
    assert c.loc["a", "2"] == 0.0                            # NaN -> 0 (diploid convention)
    assert c.loc["b", "1"] == pytest.approx(np.log2(0.6 / 1.1), abs=1e-6)
    assert c.loc["c", "2"] == -3.0 and c.loc["a", "4"] == 3.0  # clipped to [-3, 3]


def test_cn_instability_autosomes_only():
    cn = pd.DataFrame({"1": [1.0, 2.0], "2": [1.0, 1.2], "3": [3.0, 3.0], "4": [np.nan, 0.5]},
                      index=["a", "b"])
    chrom = pd.Series({"1": "1", "2": "7", "3": "X", "4": "22"})
    s = F.cn_instability(cn, chrom)
    assert s["a"] == 0.0                     # X gene ignored; NaN ignored
    assert s["b"] == pytest.approx(2 / 3)    # |2-1|>0.3 and |0.5-1|>0.3; 1.2 not


def test_mutation_burden_and_hotspot():
    d = pd.DataFrame({"1": [0, 2], "2": [0, 1], "3": [1, 0]}, index=["a", "b"], dtype=np.float32)
    assert F.mutation_burden(d).tolist() == pytest.approx([np.log1p(1), np.log1p(2)])
    assert F.hotspot_binary(d)["1"].tolist() == [0.0, 1.0]


# ----------------------------------------------------------------------------- fusions
def test_fusion_filter_drops_self_and_noncoding():
    fus = pd.DataFrame({
        "ModelID": ["L1", "L1", "L2", "L2", "L3", "L3"],
        "Gene1": ["TFAM (ENSG00000108064.11)", "BCR (ENSG00000186716.21)", "KRT7 (ENSG00000135480.15)",
                  "MID1 (ENSG00000101871.16)", "ABL1 (ENSG00000097007.19)", "X (ENSG00000000002.1)"],
        "Gene2": ["TFAM (ENSG00000108064.11)", "ABL1 (ENSG00000097007.19)", "KRT8 (ENSG00000170421.12)",
                  "RP11-120D5.1 (ENSG00000999999.1)", "BCR (ENSG00000186716.21)",
                  "Y (ENSG00000000003.1)"],
    })
    e2e = {"ENSG00000108064": "7019", "ENSG00000186716": "613", "ENSG00000097007": "25",
           "ENSG00000135480": "3855", "ENSG00000170421": "3856", "ENSG00000101871": "4281",
           "ENSG00000000002": "100", "ENSG00000000003": "100"}  # X and Y map to the same Entrez
    coding = {"7019", "613", "25", "3855", "3856", "4281", "100"}
    long, rep = F.filter_fusions(fus, coding, e2e)
    pairs = set(zip(long["ModelID"], long["a"], long["b"]))
    assert pairs == {("L1", "25", "613"), ("L2", "3855", "3856"), ("L3", "25", "613")}
    assert rep["self_fusion_rows"] == 2 and rep["noncoding_rows"] == 1
    pm, gm = F.fusion_matrices(long, ["L1", "L2", "L3", "L4"], min_lines_total=2)
    assert list(pm.columns) == ["25--613"] and pm["25--613"].tolist() == [1, 0, 1, 0]
    assert set(gm.columns) == {"25", "613"}


# ----------------------------------------------------------------------------- FoldTransformer
def _toy_raw(rng, n=40):
    ids = [f"ACH-{i:06d}" for i in range(n)]
    E = pd.DataFrame(rng.normal(size=(n, 6)), index=ids, columns=list("abcdef"), dtype=np.float32)
    E["f"] = 1.0                                          # constant -> dropped
    E["e"] = rng.normal(size=n) * 0.01                    # SD < 0.1 -> dropped
    M = pd.DataFrame((rng.random((n - 5, 4)) < [0.5, 0.05, 0.3, 0.0]).astype(np.float32),
                     index=ids[:n - 5], columns=["m1", "m2", "m3", "m4"])  # 5 lines lack M
    meta = pd.DataFrame({"OncotreeLineage": rng.choice(["Lung", "Skin", "Bowel"], n).astype(object),
                         "Age": np.where(rng.random(n) < 0.2, np.nan, rng.uniform(20, 80, n))},
                        index=ids)
    meta["OncotreeLineage"] = meta["OncotreeLineage"].astype("string")
    Pm = pd.DataFrame(rng.normal(size=(30, 3)), index=ids[:30], columns=["p1", "p2", "p3"])
    Pm.iloc[:15, 2] = np.nan                              # 50% missing -> dropped at 0.2
    Pm.iloc[0, 0] = np.nan
    return ids, {"E": E, "M": M, "meta": meta, "P": Pm}


def _specs():
    return {"E": F.BlockSpec(min_sd=0.1, absent_fill="mean"),
            "M": F.BlockSpec(kind="count", min_count=5, absent_flag=True),
            "meta": F.BlockSpec(kind="categorical"),
            "P": F.BlockSpec(max_missing=0.2, absent_fill="mean", absent_flag=True),
            "score": F.BlockSpec(kind="score", source="E", gene_sets={"s": ["a", "b", "c"]})}


def test_fold_transformer_filters_and_standardises():
    rng = np.random.default_rng(0)
    ids, raw = _toy_raw(rng)
    ft = F.FoldTransformer(_specs(), check_roles=False).fit(raw, ids[:30])
    assert ft.params_["E"]["columns"] == list("abcd")
    m_cols = ft.params_["M"]["columns"]
    assert "m4" not in m_cols and "m1" in m_cols          # never-mutated gene dropped
    assert ft.params_["P"]["columns"] == ["p1", "p2"]
    X = ft.transform(raw, ids[:30])
    e = X[[f"E:{c}" for c in "abcd"]].to_numpy()
    assert np.allclose(e.mean(0), 0, atol=1e-5) and np.allclose(e.std(0), 1, atol=1e-4)
    assert X.notna().all().all()
    Xv = ft.transform(raw, ids[30:] + ["ACH-999999"])     # unseen line: absent everywhere
    assert Xv.loc["ACH-999999", "M:absent"] == 1 and Xv.loc["ACH-999999", "P:absent"] == 1
    assert Xv.loc["ACH-999999", "score:s"] == 0.0
    assert Xv.loc["ACH-000036", "M:absent"] == 1          # no mutation data -> zero-filled + flag
    assert set(c for c in X.columns if c.startswith("meta:OncotreeLineage="))  # one-hot present
    assert list(X.columns) == list(Xv.columns)


def test_fold_transformer_never_uses_non_training_rows():
    rng = np.random.default_rng(1)
    ids, raw = _toy_raw(rng)
    train, other = ids[:25], ids[25:]
    ft1 = F.FoldTransformer(_specs(), check_roles=False).fit(raw, train)
    # perturb every non-training row of every block wildly (values, NaNs, new categories)
    raw2 = {k: v.copy() for k, v in raw.items()}
    for k, v in raw2.items():
        o = v.index.intersection(other)
        if k == "meta":
            v.loc[o, "OncotreeLineage"] = "Unseen"
            v.loc[o, "Age"] = 1e6
        else:
            v.loc[o] = (rng.normal(size=(len(o), v.shape[1])) * 1e3 + 7).astype(np.float32)
            v.loc[o[:3]] = np.nan
    ft2 = F.FoldTransformer(_specs(), check_roles=False).fit(raw2, train)
    p1, p2 = ft1.fitted_params(), ft2.fitted_params()
    assert p1.keys() == p2.keys()
    for k in p1:
        if p1[k].dtype == object:
            assert list(p1[k]) == list(p2[k]), k
        else:
            assert np.array_equal(p1[k], p2[k]), k
    pd.testing.assert_frame_equal(ft1.transform(raw, train), ft2.transform(raw2, train))
    # positive control: perturbing a TRAINING row does change the fit
    raw3 = {k: v.copy() for k, v in raw.items()}
    raw3["E"].loc[train[0], "a"] += 100
    ft3 = F.FoldTransformer(_specs(), check_roles=False).fit(raw3, train)
    assert not np.array_equal(ft3.params_["E"]["mean"], ft1.params_["E"]["mean"])


def test_fit_refuses_test_ids(tmp_path):
    sp = tmp_path / "splits.csv"
    pd.DataFrame({"ModelID": ["ACH-000000", "ACH-000001"], "role": ["dev", "test"]}).to_csv(sp, index=False)
    rng = np.random.default_rng(2)
    ids, raw = _toy_raw(rng)
    with pytest.raises(ValueError, match="forbidden"):
        F.FoldTransformer({"E": F.BlockSpec()}, splits_path=sp).fit(raw, ids[:2])
    F.FoldTransformer({"E": F.BlockSpec()}, splits_path=sp).fit(raw, [ids[0]] + ids[2:])


def test_composites_known_flags():
    lk = {s: str(i + 1) for i, s in enumerate(F.COMPOSITE_GENES)}  # TP53=1, MDM2=2, ... CDKN2A=8
    cn = pd.DataFrame(1.0, index=["a", "b"], columns=[str(i) for i in range(1, 9)])
    cn.loc["b", "1"] = 0.5          # TP53 loss
    cn.loc["b", "2"] = 3.0          # MDM2 amp
    cn.loc["b", "8"] = 0.01         # CDKN2A deep deletion
    dam = pd.DataFrame(0.0, index=["a", "b", "c"], columns=[str(i) for i in range(1, 9)])
    dam.loc["c", "7"] = 1           # RB1 damaging
    hot = pd.DataFrame({"1": [1.0]}, index=["a"])
    c = F.build_composites(cn, dam, hot, lk)
    assert c.loc["a", "TP53_hot"] == 1 and c.loc["a", "TP53_mut"] == 1 and c.loc["a", "p53_altered"] == 1
    assert c.loc["b", ["TP53_cnloss", "MDM2_amp", "CDKN2A_deepdel"]].tolist() == [1, 1, 1]
    assert c.loc["b", "MDM4_amp"] == 0 and c.loc["b", "RB_altered"] == 0
    assert c.loc["c", "RB1_dam"] == 1 and c.loc["c", "TP53_cnloss"] == 0  # no CN -> diploid


# ----------------------------------------------------------------------------- paralogs
def test_top_paralogs_ordering():
    pairs = pd.DataFrame({
        "a": ["10", "10", "10", "10", "10", "20"],
        "b": ["11", "12", "13", "14", "9", "10"],
        "identity": [30, 80, 50, 50, 15, 40], "rev_identity": [60, 10, 40, 40, 5, 80],
    })
    pairs["score"] = pairs[["identity", "rev_identity"]].max(axis=1)
    pairs["type"] = "other_paralog"
    t = P.top_paralogs(pairs, k=3)
    # scores: 11->60, 12->80, 13->50, 14->50 (tie: smaller Entrez first), 9->15
    assert t.loc["10", ["par1", "par2", "par3"]].tolist() == ["12", "11", "13"]
    assert t.loc["10", "par1_score"] == 80 and t.loc["10", "family_size"] == 5
    assert t.loc["10", "family_size_id20"] == 4 and t.loc["10", "max_identity"] == 80
    assert t.loc["20", "par1"] == "10" and pd.isna(t.loc["20", "par2"])


def test_load_paralog_pairs_symmetric_max(tmp_path):
    hdr = ["gene_stable_id", "protein_stable_id", "species", "identity", "homology_type",
           "homology_gene_stable_id", "homology_protein_stable_id", "homology_species",
           "homology_identity", "dn", "ds", "goc_score", "wga_coverage", "is_high_confidence",
           "homology_id"]
    rows = [["ENSG1", "p", "homo_sapiens", "70", "other_paralog", "ENSG2", "q", "homo_sapiens", "40"],
            ["ENSG1", "p", "homo_sapiens", "20", "within_species_paralog", "ENSG3", "q", "homo_sapiens", "25"],
            ["ENSG1", "p", "homo_sapiens", "99", "ortholog_one2one", "ENSMUSG", "q", "mus_musculus", "99"],
            ["ENSG1", "p", "homo_sapiens", "50", "gene_split", "ENSG4", "q", "homo_sapiens", "50"]]
    f = tmp_path / "c.tsv.gz"
    with gzip.open(f, "wt") as fh:
        fh.write("\t".join(hdr) + "\n")
        for r in rows:
            fh.write("\t".join(r + ["NULL"] * 5 + ["1"]) + "\n")
    hgnc = pd.DataFrame({"ensembl_gene_id": ["ENSG1", "ENSG2", "ENSG3", "ENSG4"],
                         "entrez_id": ["1", "2", "3", "4"]})
    p = P.load_paralog_pairs(f, hgnc=hgnc)
    got = {(a, b): s for a, b, s in zip(p["a"], p["b"], p["score"])}
    assert got == {("1", "2"): 70, ("2", "1"): 70, ("1", "3"): 25, ("3", "1"): 25}


# ----------------------------------------------------------------------------- local features
def test_local_features_known_values():
    lines = ["L1", "L2"]
    E = pd.DataFrame({"1": [5.0, 0.5], "2": [0.0, 3.0], "3": [2.0, 2.0], "4": [4.0, 6.0]}, index=lines)
    CN = pd.DataFrame({"1": [1.0, 0.1], "2": [0.2, 1.0]}, index=lines)
    M = pd.DataFrame({"2": [0.0, 2.0]}, index=lines)
    top = P.top_paralogs(pd.DataFrame({"a": ["1", "1", "2"], "b": ["2", "3", "1"], "identity": [50, 30, 50],
                                       "rev_identity": [50, 30, 50], "score": [50, 30, 50],
                                       "type": "other_paralog"}))
    partners = pd.DataFrame({"a": ["1", "1", "3"], "b": ["3", "4", "1"], "source": "union"})
    b = P.LocalFeatureBuilder(E, CN, M, top, partners)
    f = b.local_features(["L1", "L2", "Lx"], ["1", "3"])
    assert f["own_expr"][0].tolist() == [5.0, 2.0] and np.isnan(f["own_expr"][2]).all()
    assert f["not_expressed"][1].tolist() == [1.0, 0.0]          # 0.5 < 1
    assert f["expr_hinge"][0].tolist() == [4.0, 1.0]
    assert f["deep_del"][:, 0].tolist() == [0, 1, 0]             # 0.1 < 0.25; unknown line -> diploid
    assert f["own_cn"][0, 0] == 0.0
    assert f["par1_expr"][:2, 0].tolist() == [0.0, 3.0]          # gene 1: top paralog = 2
    assert f["par2_expr"][:2, 0].tolist() == [2.0, 2.0]          # then 3
    assert f["par3_expr"][:2, 0].tolist() == [0.0, 0.0]          # none -> 0
    assert f["par_min_expr"][:2, 0].tolist() == [0.0, 2.0]       # min over existing paralogs only
    assert f["par_min_expr"][:2, 1].tolist() == [0.0, 0.0]       # gene 3 has no paralog
    assert f["par_deep_del"][:2, 0].tolist() == [1.0, 0.0]       # paralog 2 CN 0.2 in L1
    assert f["par_mut"][:2, 0].tolist() == [0.0, 2.0]
    assert f["partner_mean_expr"][:2, 0].tolist() == [3.0, 4.0]  # mean of genes 3 and 4
    assert f["partner_mean_expr"][:2, 1].tolist() == [5.0, 0.5]
    assert b.static_features(["1", "3"])["n_partners"].tolist() == [2, 1]


# ----------------------------------------------------------------------------- proteomics
def test_build_proteomics_toy(tmp_path):
    f = tmp_path / "prot.csv"
    pd.DataFrame({
        "Protein_Id": ["sp|A|X", "tr|B|X", "sp|C|Y", "sp|D|Z"], "Gene_Symbol": ["TP53", "TP53", "OLDSYM", "NOPE"],
        "Description": "", "Group_ID": 0, "Uniprot": "", "Uniprot_Acc": "", "TenPx01_Peptides": 1,
        "L1_LUNG_TenPx01": [1.0, 9.0, 2.0, 0.0], "L1_LUNG_TenPx02": [3.0, np.nan, np.nan, 0.0],
        "L2_SKIN_TenPx01": [5.0, 9.0, 4.0, 0.0], "ZZZ_BONE_TenPx02": [0.0, 0.0, 0.0, 0.0],
    }).to_csv(f, index=False)
    model = pd.DataFrame({"ModelID": ["ACH-1", "ACH-2"], "CCLEName": ["L1_LUNG", "L2_SKIN"]})
    hgnc = pd.DataFrame({"symbol": ["TP53", "NEW"], "entrez_id": ["7157", "42"],
                         "prev_symbol": [np.nan, "OLDSYM"], "alias_symbol": [np.nan, np.nan]})
    rep = {}
    d = F.build_proteomics(f, model=model, hgnc=hgnc, report=rep)
    assert sorted(d.index) == ["ACH-1", "ACH-2"] and sorted(d.columns) == ["42", "7157"]
    assert d.loc["ACH-1", "7157"] == 2.0      # Swiss-Prot row, mean of two plexes
    assert d.loc["ACH-1", "42"] == 2.0         # NaN-aware mean
    assert rep["unmatched_names"] == ["ZZZ_BONE"]


def test_provenance_meta_is_opt_in():
    assert "meta_prov" not in F.default_specs()
    assert "meta_prov" in F.default_specs(include_provenance=True)
