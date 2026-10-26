"""Provenance + leakage check: does the COMMITTED kernel_ridge.py reproduce the stored G2 OOF
predictions, when the fit can see only the training fold's rows of the raw feature block?

The G2 runs were launched from dirty trees on commits that predate kernel_ridge.py, so the stored
predictions are not tied to committed code. Here, for one outer fold:
  1. raw expression rows of every non-training line are REMOVED from the in-process feature cache
     before fit (so neither held-out nor test features can influence any fitted statistic);
  2. the full block is restored and the held-out lines are predicted;
  3. predictions are compared with the stored cache/runs/<run_id>/pred.parquet (float32).
No registry / results writes. Dev targets only (essml.guard.load_targets).

    python reviews/g2_review/reproduce_fold.py <runs dir> <fold> <name> [<name> ...]
    names: m1_E (B1 residual), m1_E_resB0, r1_E, t_E_status
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

for v in ["OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"]:
    os.environ[v] = "6"
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).parent))
from essml import guard  # noqa: E402
from essml.models import kernel_ridge as KR  # noqa: E402
from recompute_p1 import RUNS  # noqa: E402

SPECS = {"m1_E": (KR.KernelRidgeM1, {"modalities": ["E"]}),
         "m1_E_resB0": (KR.KernelRidgeM1, {"modalities": ["E"], "residual": "B0"}),
         "r1_E": (KR.KernelRidgeR1, {"modalities": ["E"]}),
         "t_E_status": (KR.KernelRidgeT, {"modalities": ["E", "status"]})}


def main(run_dir, fold, names):
    sp = pd.read_csv(ROOT / "configs" / "splits.csv", index_col=0)
    dev = sp[sp["role"] == "dev"]
    folds = dev["cv_fold"].astype(int)
    tr, ho = folds.index[folds != fold], folds.index[folds == fold]
    Y = guard.load_targets(dev.index).astype(np.float64)
    lin = dev["OncotreeLineage"].astype(str)
    grp = dev["PatientID"].astype(str)
    out = {}
    for name in names:
        cls, params = SPECS[name]
        m = cls(**params)
        KR._RAW_CACHE.clear()
        KR.load_raw_blocks(KR.raw_names(KR.block_specs(m.modalities)))
        full = {k: v.copy() for k, v in KR._RAW_CACHE.items()}  # incl. 'composite' behind 'status'
        for k in list(KR._RAW_CACHE):  # training rows only during fit
            KR._RAW_CACHE[k] = full[k].loc[full[k].index.intersection(tr)]
        t0 = time.perf_counter()
        m.fit(Y.loc[tr], lineage=lin.loc[tr], groups=grp.loc[tr].to_numpy())
        t_fit = time.perf_counter() - t0
        for k in list(KR._RAW_CACHE):
            KR._RAW_CACHE[k] = full[k]
        P = m.predict(ho, lineage=lin.loc[ho]).to_numpy(np.float64)
        S = pd.read_parquet(Path(run_dir) / RUNS[name] / "pred.parquet").loc[ho, Y.columns].to_numpy(np.float64)
        d = np.abs(P - S)
        rel = d / np.maximum(np.abs(S), 1e-3)
        out[name] = {"fold": fold, "fit_s_6threads": round(t_fit, 1), "max_abs_diff": float(d.max()),
                     "p99_abs_diff": float(np.quantile(d, 0.99)), "max_rel_diff": float(rel.max()),
                     "lam_by_group": {str(k): v for k, v in m.lam_by_group_.items()},
                     "kappa": m.kappa_, "n_features": len(m.columns_)}
        print(name, json.dumps(out[name]), flush=True)
    p = Path(__file__).parent / f"reproduce_fold{fold}.json"
    prev = json.loads(p.read_text()) if p.exists() else {}
    prev.update(out)
    p.write_text(json.dumps(prev, indent=1))


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), sys.argv[3:])
