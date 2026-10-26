"""Dev-wide gene classes for evaluation reporting -> configs/gene_classes.csv (DESIGN v3 §3, §4).

Computed once from all 1022 dev lines (training data for the test set; identical for every model
in CV reporting). Columns: entrez, symbol, the six class flags, and the dev statistics behind them
(n_obs, sd, frac_below = fraction of observed lines < -0.5, n_below). Prevalences use observed
(non-NaN) lines only. Also prints the Hart NEGv1 -> Entrez mapping report.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from essml import genesets, guard, io  # noqa: E402


def main():
    s = pd.read_csv(ROOT / "configs" / "splits.csv", index_col=0)
    dev = s.index[s.role == "dev"]
    Y = guard.load_targets(dev).astype(np.float64)
    neg, rep = genesets.load_hart_neg()
    cls = genesets.gene_classes(Y, neg)
    obs = Y.notna()
    stats = pd.DataFrame({"n_obs": obs.sum(), "sd": Y.std(ddof=1).round(5),
                          "n_below": (Y < -0.5).sum(),
                          "frac_below": ((Y < -0.5).sum() / obs.sum()).round(5)})
    gm = io.load_gene_map().set_index("entrez")["symbol"]
    out = pd.concat([cls, stats], axis=1)
    out.insert(0, "symbol", gm.reindex(out.index).to_numpy())
    out.index.name = "entrez"
    out.reset_index().to_csv(ROOT / "configs" / "gene_classes.csv", index=False)
    print(cls.sum().to_string())
    print("NEGv1 mapping:", json.dumps({k: v for k, v in rep.items() if k != "unmapped"}),
          "unmapped:", rep["unmapped"])


if __name__ == "__main__":
    main()
