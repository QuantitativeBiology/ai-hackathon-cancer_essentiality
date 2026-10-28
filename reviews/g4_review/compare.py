"""Compare reviews/g4_review/reproduce_output.json with results/test/*/summary.json (G4 review).

    python reviews/g4_review/compare.py    # writes reviews/g4_review/comparison.csv and prints discrepancies > 1e-4
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
TOL = 1e-4

mine = json.loads((HERE / "reproduce_output.json").read_text())
rows = []
for n, d in mine["point"].items():
    s = json.loads((ROOT / "results/test" / n / "summary.json").read_text())
    h = s["core"]["headline"]
    sv = s["core"]["summary"]
    pairs = {}
    for c in ("all", "variable", "selective", "rare_selective"):
        pairs[f"P1_{c}"] = h.get(f"P1_{c}")
        pairs[f"n_genes_{c}"] = sv[c]["n_genes"]
    pairs["P1_median_variable"] = sv["variable"]["P1_median"]
    pairs["n_r_gt_0.3_variable"] = sv["variable"]["n_r_gt_0.3"]
    for b in ("B0", "B1", "B2"):
        pairs[f"dP1_vs_{b}_variable"] = h.get(f"dP1_vs_{b}_variable")
        pairs[f"frac_beat_{b}_variable"] = sv["variable"].get(f"frac_beat_{b}")
    for k in ("P2a", "dP2a_vs_B0", "P2b_nap_prev05", "P2b_auroc_prev05"):
        pairs[k] = h.get(k)
    pairs["P2b_n_prev05"] = sv["selective"]["P2b_n_prev05"]
    for tname, key in (("primary", "per_tier"), ("m1_own", "per_tier_own_table")):
        if f"tiers_{tname}" not in d:
            continue
        pt = s[key]
        for t in ("high", "medium", "low"):
            for k in ("n_genes", "P1", "P1_variable", "P2b_nap_prev05", "P2b_auroc_prev05", "P2b_n_prev05"):
                pairs[f"{tname}:{t}:{k}"] = pt["tiers"][t][k]
        for lab in ("evaluable", "variable"):
            pairs[f"{tname}:spearman_{lab}"] = (pt["spearman_dev_mean_r_vs_heldout_r"].get(lab) or {}).get("point")
    flat = dict(d)
    for tname in ("primary", "m1_own"):
        for t, v in d.get(f"tiers_{tname}", {}).items():
            if t == "spearman_dev_mean_r_vs_test_r":
                flat[f"{tname}:spearman_evaluable"] = v["evaluable"]
                flat[f"{tname}:spearman_variable"] = v["variable"]
            else:
                for k, x in v.items():
                    flat[f"{tname}:{t}:{k}"] = x
    for k, theirs in pairs.items():
        m = flat.get(k)
        if m is None and theirs is None:
            continue
        diff = None if (m is None or theirs is None) else float(m) - float(theirs)
        rows.append({"finalist": n, "metric": k, "reproduced": m, "project": theirs, "diff": diff})

df = pd.DataFrame(rows)
df.to_csv(HERE / "comparison.csv", index=False, float_format="%.6g")
bad = df[df["diff"].abs() > TOL]
print(f"{len(df)} comparisons; {int(df['diff'].notna().sum())} numeric; max |diff| {df['diff'].abs().max():.3g}")
print("discrepancies > 1e-4:" if len(bad) else "no discrepancy > 1e-4")
if len(bad):
    print(bad.to_string())
