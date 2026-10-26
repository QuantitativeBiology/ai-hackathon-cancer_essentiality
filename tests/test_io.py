import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import numpy as np
import pandas as pd

from essml import io


def test_parse():
    assert io.parse_gene_column("H3-3A (3020)") == ("H3-3A", "3020")
    assert io.parse_gene_column("junk") is None


def test_convert_default_entries(tmp_path):
    p = tmp_path / "expr.csv"
    pd.DataFrame({
        "": [0, 1, 2], "SequencingID": ["a", "b", "c"], "ModelConditionID": ["m1", "m1b", "m2"],
        "ModelID": ["ACH-1", "ACH-1", "ACH-2"], "IsDefaultEntryForMC": ["Yes"] * 3,
        "IsDefaultEntryForModel": ["No", "Yes", "Yes"],
        "AAA (10)": [1.0, 2.0, 3.0], "BBB (20)": [4.0, 5.0, 6.0], "AAA2 (10)": [9, 9, 9],
        "weird": [0, 0, 0],
    }).to_csv(p, index=False)
    mat, info = io.convert_matrix(p, "expression", chunksize=2)
    assert list(mat.index) == ["ACH-1", "ACH-2"]
    assert list(mat.columns) == ["10", "20"]
    assert mat.loc["ACH-1", "10"] == 2.0 and mat.loc["ACH-2", "20"] == 6.0
    assert (mat.dtypes == np.float32).all()
    assert info["unparseable"] == ["weird"] and info["n_dup_entrez"] == 1
    assert info["symbols"]["10"] == "AAA"


def test_crispr_layout(tmp_path):
    p = tmp_path / "c.csv"
    pd.DataFrame({"": ["ACH-1", "ACH-2"], "AAA (10)": [-1.0, np.nan]}).to_csv(p, index=False)
    mat, _ = io.convert_matrix(p, "crispr")
    assert list(mat.index) == ["ACH-1", "ACH-2"] and np.isnan(mat.loc["ACH-2", "10"])
