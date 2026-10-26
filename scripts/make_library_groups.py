"""Library-proxy line groups for dev lines (D16/D17 report column; G3 stratified controls).

    python scripts/make_library_groups.py      # -> configs/library_groups.csv

Definition (reused verbatim from the independent G2 review, reviews/g2_review/recompute_p1.py
lines 94-98): the row NaN pattern of a dev line's gene effects over ALL genes; the 4 most frequent
patterns are groups 0-3 (by decreasing size) and every other line is group 4 ("other"). This is
the NaN/library line group of DESIGN §2 ("kept as a library proxy for stratified reporting").

Only the NaN mask of dev lines is used (targets read through essml.guard with dev IDs; values are
never written). The file holds ModelID -> group label only.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from essml import guard  # noqa: E402

OUT = ROOT / "configs" / "library_groups.csv"


def library_groups(Y: pd.DataFrame, n_top: int = 4) -> pd.Series:
    """Group label per line: rank of its row NaN pattern among the n_top largest, else n_top."""
    key = Y.isna().apply(lambda r: hash(r.to_numpy().tobytes()), axis=1)
    top = key.value_counts().index[:n_top]
    lab = np.where(key.isin(top), key.map({k: i for i, k in enumerate(top)}), n_top).astype(int)
    return pd.Series(lab, index=Y.index, name="lib_group")


def main():
    sp = pd.read_csv(ROOT / "configs" / "splits.csv", index_col=0)
    dev = sp[sp["role"] == "dev"]
    Y = guard.load_targets(dev.index)
    lab = library_groups(Y)
    n_nan = Y.isna().sum(1).rename("n_nan_genes")
    out = pd.concat([lab, n_nan], axis=1)
    out.index.name = "ModelID"
    out.to_csv(OUT)
    print(out["lib_group"].value_counts().sort_index().to_dict())


if __name__ == "__main__":
    main()
