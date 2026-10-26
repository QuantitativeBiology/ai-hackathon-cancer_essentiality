"""Remove held-out lines from the committed CRISPR audit tables (red-team G0 finding #3).

The pre-split audit (scripts/audit_data.py, legacy) summarised CRISPR values over all 1361 lines. Committed tables
must not carry statistics of test or non-cancerous lines, so:
- crispr_line_stats.csv keeps dev lines only;
- crispr_gene_stats.csv is recomputed on dev lines only, through the guard.
Run after scripts/audit_data.py.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from essml import guard  # noqa: E402

THR = -0.5
out = ROOT / "results" / "data_audit"
splits = pd.read_csv(ROOT / "configs" / "splits.csv", index_col=0)
dev = splits.index[splits["role"] == "dev"]

line = pd.read_csv(out / "crispr_line_stats.csv", index_col=0)
line.loc[line.index.intersection(dev)].to_csv(out / "crispr_line_stats.csv")

Y = guard.load_targets(dev)
v = Y.to_numpy(dtype=np.float64)
obs = np.isfinite(v)
gene = pd.DataFrame({
    "entrez": Y.columns,
    "mean": np.nanmean(v, 0),
    "sd": np.nanstd(v, 0, ddof=1),
    "frac_below_thr": (np.where(obs, v, 0) < THR).sum(0) / obs.sum(0),
    "n_nan": (~obs).sum(0),
})
gene.to_csv(out / "crispr_gene_stats.csv", index=False)
print(f"line stats: {len(line)} -> {len(line.index.intersection(dev))} dev rows; gene stats on {len(dev)} dev lines")
