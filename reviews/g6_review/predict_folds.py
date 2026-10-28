"""G6 independent review, step 0: per-fold predictions of the PERSISTED fold models (no refit).

The saved table results/transfer/predictions_with_confidence_B_*.parquet holds only the 5-fold MEAN, but the
primary statistic (D7) averages skill over fold models. So the persisted fold pickles are loaded and asked to
predict the same declared inputs; the fold mean is checked against the saved table (validates both).
Model code (KernelRidgeT.predict) is used only to evaluate the frozen models; no evaluation code is imported.
Output: per-fold predictions (no observed values) in this scratch directory.
"""
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(r"C:\Users\ritab\Desktop\project")
sys.path.insert(0, str(ROOT / "src"))
OUT = Path(__file__).resolve().parent

ext = pd.read_parquet(ROOT / "cache/external_expr/harmonised.parquet")
ext.columns = ext.columns.astype(str)
E = pd.read_parquet(ROOT / "cache/features/E.parquet")
E.columns = E.columns.astype(str)
E.index = E.index.astype(str)

SETA = ["ACH-002462", "ACH-002463", "ACH-002464", "ACH-002465", "ACH-002467"]
JOBS = {  # role dir -> list of (profile_id, tp53)
    "t_E_status_24b46e2a": [("state:DepMap_RPE1", "WT"), ("state:DepMap_RPE1", "null"),
                            ("state:H1_iCas9_MEF", "WT"), ("state:H1_iCas9_MEF", "null"),
                            ("state:H1_iCas9_laminin", "WT"), ("state:H1_iCas9_laminin", "null")]
                           + [(f"line:{l}", "WT") for l in SETA],
    "t_E_status_noEC_fa0181c1": [("state:H1_iCas9_MEF", "WT"), ("state:H1_iCas9_laminin", "WT")],
}


def x_for(pid, tp53):
    row = E.loc[pid[5:]] if pid.startswith("line:") else ext.loc[pid]
    key = f"profile::{pid}::{tp53}"
    Xe = pd.DataFrame([row.to_numpy(np.float64)], index=[key], columns=row.index.astype(str))
    st = pd.DataFrame({"TP53_mut": [1.0 if tp53 == "null" else 0.0], "RB_altered": [0.0]}, index=[key])
    return key, {"E": Xe, "status": st}


rows = []
for tag, jobs in JOBS.items():
    for p in sorted((ROOT / "cache/models" / tag).glob("fold*.pkl")):
        obj = pickle.load(open(p, "rb"))
        assert obj["dropped"] == [], obj["dropped"]
        m = obj["cls"].__new__(obj["cls"])
        m.__dict__.update(obj["state"])
        k = obj["fold"]
        # leakage check: training ids of the fold model contain no set-A / non-dev line
        sp = pd.read_csv(ROOT / "configs/splits.csv", index_col=0)
        roles = sp["role"].reindex(obj["train_ids"])
        assert (roles == "dev").all(), roles.value_counts(dropna=False)
        for pid, tp in jobs:
            key, X = x_for(pid, tp)
            pr = m.predict([key], lineage=None, X=X).iloc[0]
            rows.append(pd.DataFrame({"tag": tag, "fold": k, "profile_id": pid, "tp53": tp,
                                      "entrez": pr.index.astype(str), "pred": pr.to_numpy(np.float64)}))
        print(tag, k, "n_train", len(obj["train_ids"]), flush=True)
P = pd.concat(rows, ignore_index=True)
P.to_parquet(OUT / "fold_preds.parquet")

# check fold mean vs the saved G6 tables
for f in ["B_20261028-024139", "A_20261028-014928"]:
    S = pd.read_parquet(ROOT / f"results/transfer/predictions_with_confidence_{f}.parquet")
    S = S[S.fold_models == "all"]
    for (tag, pid, tp), g in P.groupby(["tag", "profile_id", "tp53"]):
        fin = {"t_E_status_24b46e2a": "t_E_status", "t_E_status_noEC_fa0181c1": "t_E_status_noEC"}[tag]
        s = S[(S.finalist == fin) & (S.profile_id == pid) & (S.tp53_input == tp)]
        if s.empty:
            continue
        mean = g.groupby("entrez")["pred"].mean()
        sv = s.set_index("entrez")["predicted_effect"]
        d = np.abs(mean.reindex(sv.index).to_numpy() - sv.to_numpy())
        print(f, fin, pid, tp, "n", len(sv), "max|fold-mean - saved|", float(np.nanmax(d)), "n_folds", g.fold.nunique())
