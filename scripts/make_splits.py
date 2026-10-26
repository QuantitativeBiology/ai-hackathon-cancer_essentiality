"""Define the locked cancer test split and the development CV folds.

Population: screened models (CRISPRGeneEffect rows) that have default-entry expression and are not
'Non-Cancerous' (those form the separate non-transformed evaluation set). Splitting is grouped by
PatientID (models from one patient never straddle train/test) and stratified by OncotreeLineage.
Output: configs/splits.csv with ModelID, PatientID, OncotreeLineage, role (test/dev/noncancerous), cv_fold.
"""
import argparse
import json
import sys
from pathlib import Path

import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from essml.io import load_matrix  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model-csv", default="data/Model.csv")
    ap.add_argument("--seeds", default="configs/seeds.json")
    ap.add_argument("--out", default="configs/splits.csv")
    ap.add_argument("--test-folds", type=int, default=5, help="1/k of the population becomes the locked test set")
    ap.add_argument("--cv-folds", type=int, default=5)
    a = ap.parse_args()
    seed = json.load(open(a.seeds))["split_seed"]

    model = pd.read_csv(a.model_csv).set_index("ModelID")
    screened = load_matrix("crispr", unguarded=True).index  # pre-split: index only
    expr = load_matrix("expression").index
    ids = screened.intersection(expr)
    m = model.loc[ids, ["PatientID", "OncotreeLineage", "OncotreePrimaryDisease"]].copy()
    # A model with no PatientID is its own group.
    m["PatientID"] = m["PatientID"].fillna(pd.Series(m.index, index=m.index))
    nonc = m["OncotreePrimaryDisease"].eq("Non-Cancerous")
    # Patients with any non-cancerous model are kept out of the cancer population entirely.
    nonc_patients = set(m.loc[nonc, "PatientID"])
    m["role"] = "dev"
    m.loc[m["PatientID"].isin(nonc_patients), "role"] = "noncancerous"
    canc = m[m.role == "dev"]
    # Rare lineages are pooled for stratification only.
    lin = canc["OncotreeLineage"].fillna("Other")
    vc = lin.value_counts()
    strat = lin.where(lin.map(vc) >= a.test_folds * 2, "_rare")
    sgk = StratifiedGroupKFold(n_splits=a.test_folds, shuffle=True, random_state=seed)
    _, test_idx = next(sgk.split(canc, strat, groups=canc["PatientID"]))
    m.loc[canc.index[test_idx], "role"] = "test"

    dev = m[m.role == "dev"]
    lin = dev["OncotreeLineage"].fillna("Other")
    vc = lin.value_counts()
    strat = lin.where(lin.map(vc) >= a.cv_folds * 2, "_rare")
    m["cv_fold"] = -1
    cv = StratifiedGroupKFold(n_splits=a.cv_folds, shuffle=True, random_state=seed + 1)
    for k, (_, vi) in enumerate(cv.split(dev, strat, groups=dev["PatientID"])):
        m.loc[dev.index[vi], "cv_fold"] = k

    # Integrity checks: no patient spans roles or folds.
    assert m.groupby("PatientID")["role"].nunique().max() == 1
    assert m[m.role == "dev"].groupby("PatientID")["cv_fold"].nunique().max() == 1
    m.drop(columns="OncotreePrimaryDisease").sort_index().to_csv(a.out)
    print(m["role"].value_counts().to_string())
    print(m[m.role == "dev"]["cv_fold"].value_counts().sort_index().to_string())


if __name__ == "__main__":
    main()
