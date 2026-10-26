"""Toy-data tests with known answers for essml.expr_harmonise (DESIGN v3 section 7.3, section 8)."""
import gzip
import io
import tarfile

import numpy as np
import pandas as pd
import pytest

from essml import expr_harmonise as H


# ----------------------------------------------------------------------------- annotation / IDs
def test_merged_length_union_of_exons():
    assert H.merged_length([], []) == 0
    assert H.merged_length([1], [100]) == 100
    # overlapping (1-100, 50-150) -> 1-150; nested 60-70 adds nothing; disjoint 201-210 adds 10
    assert H.merged_length([50, 1, 60, 201], [150, 100, 70, 210]) == 160
    # adjacent intervals are merged without double counting
    assert H.merged_length([1, 11], [10, 20]) == 20


def test_gene_table_from_toy_gtf(tmp_path):
    gtf = tmp_path / "toy.gtf.gz"
    lines = [
        "#!genome-build toy",
        '1\tx\tgene\t1\t300\t.\t+\t.\tgene_id "ENSG00000000001"; gene_name "AAA"; gene_biotype "protein_coding";',
        '1\tx\texon\t1\t100\t.\t+\t.\tgene_id "ENSG00000000001"; transcript_id "T1";',
        '1\tx\texon\t51\t150\t.\t+\t.\tgene_id "ENSG00000000001"; transcript_id "T2";',
        '1\tx\texon\t201\t300\t.\t+\t.\tgene_id "ENSG00000000001"; transcript_id "T2";',
        '2\tx\tgene\t1\t50\t.\t-\t.\tgene_id "ENSG00000000002"; gene_biotype "lncRNA";',
        '2\tx\texon\t1\t50\t.\t-\t.\tgene_id "ENSG00000000002"; transcript_id "T3";',
    ]
    with gzip.open(gtf, "wt") as fh:
        fh.write("\n".join(lines) + "\n")
    t = H.gene_table_from_gtf(gtf)
    assert t.at["ENSG00000000001", "length"] == 250       # 1-150 union + 201-300
    assert t.at["ENSG00000000002", "length"] == 50
    assert t.at["ENSG00000000001", "gene_name"] == "AAA"
    assert t.at["ENSG00000000002", "gene_biotype"] == "lncRNA"


def test_strip_version():
    out = H.strip_version(["ENSG00000000003.14", "ENSG00000000005", "ENSG00000182378.15_PAR_Y", "ERCC-00002"])
    assert list(out) == ["ENSG00000000003", "ENSG00000000005", "ENSG00000182378", "ERCC-00002"]


# ----------------------------------------------------------------------------- unit conversion
def test_counts_to_tpm_known_lengths():
    counts = pd.DataFrame({"s": [10.0, 20.0, 30.0, 5.0]}, index=["a", "b", "c", "nolen"])
    lengths = pd.Series({"a": 1000, "b": 2000, "c": 3000})
    tpm, used = H.counts_to_tpm(counts, lengths)
    # equal rates (0.01 per base) -> equal TPM; gene without a length is excluded
    assert list(tpm.index) == ["a", "b", "c"]
    assert tpm["s"].to_numpy() == pytest.approx([1e6 / 3] * 3)
    assert used["s"] == pytest.approx(60 / 65)
    # effective length = L - fragment + 1: 1000-501+1 = 500, 2000-501+1 = 1500, 3000-501+1 = 2500
    tpm2, _ = H.counts_to_tpm(counts, lengths, fragment_length=501)
    r = np.array([10 / 500, 20 / 1500, 30 / 2500])
    assert tpm2["s"].to_numpy() == pytest.approx(r / r.sum() * 1e6)
    # effective length <= 0 -> no TPM (RSEM convention)
    tpm3, _ = H.counts_to_tpm(counts, lengths, fragment_length=1500)
    assert "a" not in tpm3.index and tpm3["s"].sum() == pytest.approx(1e6)


def test_fpkm_to_tpm_and_renormalise():
    f = pd.DataFrame({"x": [1.0, 3.0], "y": [2.0, 2.0]}, index=["g1", "g2"])
    t = H.scale_to_million(f)
    assert t["x"].tolist() == pytest.approx([2.5e5, 7.5e5])
    assert t["y"].tolist() == pytest.approx([5e5, 5e5])


def test_to_tpm_drops_spikeins_for_rsem_tables():
    v = pd.DataFrame({"s": [400000.0, 400000.0, 200000.0]},
                     index=["ENSG00000000001", "ENSG00000000002", "ERCC-00002"])
    raw = H.RawTable(v, "ensembl", "tpm")
    tpm, info = H.to_tpm(raw, genes=pd.DataFrame({"length": []}), mapper=None)
    assert list(tpm.index) == ["ENSG00000000001", "ENSG00000000002"]
    assert tpm["s"].tolist() == pytest.approx([5e5, 5e5])
    assert info["frac_tpm_spikein_or_nongene"]["s"] == pytest.approx(0.2)


def test_project_to_depmap_sum_absent_and_log():
    tpm = pd.DataFrame({"s": [3.0, 4.0, 7.0, 100.0]}, index=["E1", "E2", "E3", "Enc"])
    proj = {"E1": "10", "E2": "10", "E3": "20"}          # E1+E2 -> Entrez 10; 'Enc' not a DepMap gene
    x = H.project_to_depmap(tpm, proj, ["10", "20", "30"], absent="nan")
    assert x.loc["s", "10"] == pytest.approx(np.log2(8.0))
    assert x.loc["s", "20"] == pytest.approx(3.0)
    assert np.isnan(x.loc["s", "30"])
    z = H.project_to_depmap(tpm, proj, ["10", "20", "30"], absent="zero")
    assert z.loc["s", "30"] == 0.0
    # restricted universe: PC subset rescaled to the target sum (7 + 7 -> 28: factor 2)
    r = H.project_to_depmap(tpm, proj, ["10", "20", "30"], absent="zero", pc_target=28.0)
    assert r.loc["s", "10"] == pytest.approx(np.log2(15.0)) and r.loc["s", "20"] == pytest.approx(np.log2(15.0))


def test_parse_htseq_tar_drops_summary_rows(tmp_path):
    body = "ENSG00000000003.14\t10\nENSG00000000005.5\t0\n__no_feature\t99\n__ambiguous\t5\n"
    tar_path = tmp_path / "RAW.tar"
    with tarfile.open(tar_path, "w") as tf:
        for gsm in ["GSM1", "GSM2"]:
            data = gzip.compress(body.encode())
            ti = tarfile.TarInfo(f"{gsm}_x_htseqReadCounts.txt.gz")
            ti.size = len(data)
            tf.addfile(ti, io.BytesIO(data))
    raw = H.parse_htseq_tar(tar_path, ["GSM2", "GSM1"])
    assert list(raw.values.columns) == ["GSM2", "GSM1"]
    assert list(raw.values.index) == ["ENSG00000000003", "ENSG00000000005"]
    assert raw.values.loc["ENSG00000000003", "GSM1"] == 10 and raw.unit == "counts"


# ----------------------------------------------------------------------------- harmonisation
def test_quantile_normalisation_maps_to_reference_distribution():
    rng = np.random.default_rng(0)
    ref = rng.gamma(2.0, 2.0, size=(20, 50))
    ref_q = H.reference_quantiles(ref)
    assert np.all(np.diff(ref_q) >= 0)
    x = pd.DataFrame(rng.normal(size=(3, 50)) * 7 + 100)          # different scale and location
    q = H.quantile_normalise(x, ref_q)
    for i in range(3):
        assert np.sort(q.iloc[i].to_numpy()) == pytest.approx(ref_q)                 # exact distribution
        assert np.array_equal(np.argsort(q.iloc[i].to_numpy()), np.argsort(x.iloc[i].to_numpy()))  # order kept


def test_quantile_normalisation_ties_and_missing():
    ref_q = np.arange(10, dtype=float)
    x = pd.DataFrame([[0, 0, 0, 0, 5, 6, 7, 8, 9, 10]], dtype=float)
    q = H.quantile_normalise(x, ref_q).iloc[0].to_numpy()
    assert np.unique(q[:4]).size == 1                         # tied zeros map to one value
    assert q[:4][0] == pytest.approx(1.5)                     # the mean-rank position of ranks 1-4
    assert q[4:] == pytest.approx(ref_q[4:])
    xm = pd.DataFrame([[np.nan, 1.0, 2.0, 3.0, 4.0, np.nan, 6.0, 7.0, 8.0, 9.0]])
    qm = H.quantile_normalise(xm, ref_q).iloc[0]
    assert qm.isna().sum() == 2                               # missing stays missing
    v = qm.dropna().to_numpy()
    assert np.all(np.diff(v) > 0) and v.min() >= 0 and v.max() <= 9


def test_rank_transform():
    r = H.rank_transform(pd.DataFrame([[3.0, 1.0, 2.0, np.nan]]))
    assert r.iloc[0, :3].tolist() == pytest.approx([1.0, 1 / 3, 2 / 3]) and np.isnan(r.iloc[0, 3])


def test_fit_and_apply_location():
    ext = pd.DataFrame([[1.0, 2.0], [3.0, np.nan]], columns=["g1", "g2"])
    dep = pd.DataFrame([[2.0, 2.5], [5.0, 1.0]], columns=["g1", "g2"])
    off = H.fit_location(ext, dep)
    assert off["g1"] == pytest.approx(1.5) and off["g2"] == pytest.approx(0.5)
    assert H.apply_location(ext, off).iloc[0].tolist() == pytest.approx([2.5, 2.5])


def test_heldout_location_never_uses_left_out_line():
    rng = np.random.default_rng(1)
    idx = ["L1a", "L1b", "L2", "L3", "L4"]
    groups = ["L1", "L1", "L2", "L3", "L4"]               # L1 has two replicate profiles
    ext = pd.DataFrame(rng.normal(size=(5, 30)), index=idx)
    dep = pd.DataFrame(rng.normal(size=(5, 30)), index=idx)
    audit = []
    out = H.heldout_location(ext, dep, groups, audit=audit)
    for g, train in audit:
        held = {i for i, gg in zip(idx, groups) if gg == g}
        assert not held & set(train)                      # none of the held-out line's rows in its fit
    # perturbing the held-out line's own DepMap profile does not change its corrected profile
    dep2 = dep.copy()
    dep2.loc[["L1a", "L1b"]] += 100.0
    out2 = H.heldout_location(ext, dep2, groups)
    assert np.allclose(out.loc[["L1a", "L1b"]], out2.loc[["L1a", "L1b"]])
    assert not np.allclose(out.loc["L2"], out2.loc["L2"])          # ...but it does enter other lines' fits
    # known answer: L2's offsets = mean(dep - ext) over L1a, L1b, L3, L4
    tr = ["L1a", "L1b", "L3", "L4"]
    expect = ext.loc["L2"] + (dep.loc[tr] - ext.loc[tr]).mean()
    assert np.allclose(out.loc["L2"], expect)
    # fit_rows restricts the fitting population (e.g. exact dev matches only)
    out3 = H.heldout_location(ext, dep, groups, fit_rows=[True, True, True, False, True])
    tr = ["L1a", "L1b", "L4"]
    assert np.allclose(out3.loc["L2"], ext.loc["L2"] + (dep.loc[tr] - ext.loc[tr]).mean())


def test_location_correction_undoes_known_genewise_shift():
    """Positive control: a planted gene-wise shift is removed by the held-out correction."""
    rng = np.random.default_rng(2)
    dep = pd.DataFrame(rng.normal(5, 2, size=(8, 200)))
    shift = rng.normal(0, 1, size=200)
    ext = dep - shift + rng.normal(0, 0.05, size=(8, 200))
    out = H.heldout_location(ext, dep, list(range(8)))
    assert (out - dep).abs().to_numpy().mean() < 0.1 < (ext - dep).abs().to_numpy().mean()


# ----------------------------------------------------------------------------- metrics
def test_corr_nn_recall_and_applicability():
    rng = np.random.default_rng(3)
    pool = pd.DataFrame(rng.normal(size=(6, 100)), index=[f"P{i}" for i in range(6)])
    x = pool.loc[["P2", "P4"]] + rng.normal(0, 0.1, size=(2, 100))
    x.index = ["e2", "e4"]
    x.iloc[0, :5] = np.nan                                  # missing genes are skipped
    r = H.corr_to_pool(x, pool)
    assert r.loc["e2", "P2"] > 0.95 and r.loc["e4", "P4"] > 0.95
    nn = H.nn_recall(r, pd.Series({"e2": "P2", "e4": "P0"}))
    assert nn.loc["e2", "true_rank"] == 1 and bool(nn.loc["e2", "recall1"])
    assert nn.loc["e4", "true_rank"] > 1 and nn.loc["e4", "nearest"] == "P4"
    lin = pd.Series({"P0": "A", "P1": "A", "P2": "B", "P3": "B", "P4": "C", "P5": "A"})
    ad = H.applicability(r, lin, k=1)
    assert ad.loc["e2", "nearest_lineage"] == "B" and ad.loc["e2", "nearest_line"] == "P2"
    assert ad.loc["e2", "ad_distance"] == pytest.approx(1 - r.loc["e2", "P2"])
    ad3 = H.applicability(r, lin, k=3)
    assert ad3.loc["e2", "ad_distance"] == pytest.approx((1 - r.loc["e2"].sort_values(ascending=False)[:3]).mean())


def test_pool_self_ad_excludes_self():
    R = pd.DataFrame([[1.0, 0.9, 0.1], [0.9, 1.0, 0.5], [0.1, 0.5, 1.0]], index=list("abc"), columns=list("abc"))
    s = H.pool_self_ad(R, k=1)
    assert s.tolist() == pytest.approx([0.1, 0.1, 0.5])
