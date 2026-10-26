"""Fixed evaluable gene set -> configs/eval_genes.csv (DESIGN §4; PI addendum, review round 2).

A gene is evaluable iff it has >= 40 observed (non-NaN) CRISPR lines in EVERY dev CV fold and
>= 50 observed lines in the test set. Test-set counts come from guard.observed_counts (mask only,
no values; logged as access_type='mask_only'). Only the gene list is written, not the counts.
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from essml import guard, io  # noqa: E402

MIN_FOLD, MIN_TEST = 40, 50


def main():
    s = pd.read_csv(ROOT / "configs" / "splits.csv", index_col=0)
    dev = s[s.role == "dev"]
    fold_counts = pd.concat({f: guard.observed_counts(dev.index[dev.cv_fold == f])
                             for f in sorted(dev.cv_fold.unique())}, axis=1)
    test_counts = guard.observed_counts(s.index[s.role == "test"])
    min_fold = fold_counts.min(axis=1)
    keep = (min_fold >= MIN_FOLD) & (test_counts.reindex(min_fold.index) >= MIN_TEST)
    gm = io.load_gene_map().set_index("entrez")["symbol"]
    out = pd.DataFrame({"entrez": min_fold.index[keep]})
    out["symbol"] = gm.reindex(out["entrez"]).to_numpy()
    out.to_csv(ROOT / "configs" / "eval_genes.csv", index=False)
    print(f"genes {len(min_fold)}; fail dev-fold rule {(min_fold < MIN_FOLD).sum()}; "
          f"fail test rule {(test_counts < MIN_TEST).sum()}; evaluable {keep.sum()}")


if __name__ == "__main__":
    main()
