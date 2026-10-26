"""Gene classes on a toy matrix with known membership (DESIGN §4)."""
import numpy as np
import pandas as pd

from essml.genesets import gene_classes, load_hart_neg, NEG_PATH


def test_gene_classes_toy():
    n = 100
    Y = pd.DataFrame(index=[f"L{i}" for i in range(n)])
    Y["ce"] = -1.0; Y.loc["L0", "ce"] = 0.0                        # 99% below -> common essential
    Y["sel"] = 0.0; Y.iloc[:30, 2 - 1] = -1.0                       # 30% -> selective, SD ~0.46
    Y["rare"] = 0.0; Y.iloc[:3, 2] = -0.6                           # 3 lines (3%) -> rare
    Y["two"] = 0.0; Y.iloc[:2, 3] = -0.6                            # 2 lines -> none
    Y["edge5"] = 0.0; Y.iloc[:5, 4] = -0.6                          # exactly 5% -> selective
    Y["nan"] = np.nan; Y.iloc[:10, 5] = -1.0                        # all observed below -> CE
    Y["neg"] = 0.0
    Y.iloc[50:, 1] = np.nan  # 'sel': 30/50 observed below -> 60%: prevalence over observed only
    c = gene_classes(Y, neg_entrez={"neg"})
    assert c.loc["ce"].to_dict() == {"all": True, "variable": False, "selective": False,
                                     "rare_selective": False, "common_essential": True,
                                     "nonessential": False}
    assert c.loc["sel", "selective"] and c.loc["sel", "variable"]
    assert c.loc["rare", "rare_selective"] and not c.loc["rare", "selective"]
    assert not c.loc["two"][["selective", "rare_selective", "common_essential"]].any()
    assert c.loc["edge5", "selective"] and not c.loc["edge5", "rare_selective"]
    assert c.loc["nan", "common_essential"]
    assert c.loc["neg", "nonessential"] and c["nonessential"].sum() == 1
    assert c["all"].all()


def test_training_rows_only():
    Y = pd.DataFrame({"g": [0.0] * 90 + [-1.0] * 10})
    assert gene_classes(Y, neg_entrez=set()).loc["g", "selective"]
    assert not gene_classes(Y.iloc[:90], neg_entrez=set()).loc["g", "selective"]


def test_hart_neg_mapping():
    if not NEG_PATH.exists():
        import pytest
        pytest.skip("NEGv1 not downloaded")
    s, rep = load_hart_neg()
    assert rep["n_mapped"] >= 0.97 * rep["n_neg"] and all(x.isdigit() for x in s)


def test_rare_selective_needs_one_percent():
    # 3 of 400 lines (0.75%) below -0.5 is noise-level, not rare-selective; 8 of 400 (2%) is.
    Y = pd.DataFrame(0.0, index=[f"L{i}" for i in range(400)], columns=["noise", "rare"])
    Y.iloc[:3, 0] = -0.6
    Y.iloc[:8, 1] = -0.6
    c = gene_classes(Y, neg_entrez=set())
    assert not c.loc["noise", "rare_selective"]
    assert c.loc["rare", "rare_selective"]
