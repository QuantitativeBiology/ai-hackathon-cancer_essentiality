"""Post-hoc check (2026-10-28, after G7): do DepMap's screen maps change anything already reported?

    python scripts/library_map_check.py      # -> results/library_map_check.json

The user supplied CRISPRScreenMap.csv, ScreenSequenceMap.csv and the 26Q3 README after the report was final.
Everything reported used configs/library_groups.csv, a library *proxy*: the row NaN pattern of a dev line's gene
effects (top 4 patterns = groups 0-3, the rest = group 4; scripts/make_library_groups.py). This script maps each
model to the screening libraries of the screens DepMap combined into its gene effects (CRISPRScreenMap ->
ScreenSequenceMap.Library) and cross-tabulates them with the proxy groups. No gene-effect values are read; no
result is re-run. Test lines enter only as aggregate library counts (screen metadata, no targets).
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "results" / "library_map_check.json"


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def _any_true(s):
    """True if any sequence is flagged True (the flag columns hold True / False / NaN; NaN is not True)."""
    return bool((s == True).any())  # noqa: E712


def model_libraries():
    """ModelID -> '+'-joined sorted set of libraries of the screens combined for it; plus screen-level facts."""
    sm = pd.read_csv(DATA / "CRISPRScreenMap.csv", dtype=str)
    seq = pd.read_csv(DATA / "ScreenSequenceMap.csv", dtype={"ScreenID": str, "ModelID": str})
    used = seq[seq["ScreenID"].isin(set(sm["ScreenID"]))]
    scr = used.groupby("ScreenID").agg(Library=("Library", lambda s: "+".join(sorted(set(s.dropna())))),
                                       Days=("Days", "median"), ScreenType=("ScreenType", "first"),
                                       DrugTreated=("DrugTreated", _any_true), IsEngineered=("IsEngineered", _any_true))
    m = sm.join(scr, on="ScreenID")
    lib = m.groupby("ModelID")["Library"].agg(lambda s: "+".join(sorted(set("+".join(s.dropna()).split("+")) - {""})))
    per_model = m.groupby("ModelID").agg(n_screens=("ScreenID", "nunique"), days=("Days", "median"),
                                         screen_types=("ScreenType", lambda s: "+".join(sorted(set(s.dropna())))),
                                         drug_treated=("DrugTreated", "max"), engineered=("IsEngineered", "max"))
    return lib.rename("library"), per_model, {"n_screen_map_rows": int(len(sm)), "n_models": int(sm["ModelID"].nunique()),
                                              "n_screens": int(sm["ScreenID"].nunique()),
                                              "n_sequences_of_mapped_screens": int(len(used)),
                                              "screens_without_sequence_rows": int((~sm["ScreenID"].isin(seq["ScreenID"])).sum())}


def _set_a(ids=tuple(f"ACH-00{n}" for n in range(2462, 2468))) -> dict:
    sm = pd.read_csv(DATA / "CRISPRScreenMap.csv", dtype=str)
    seq = pd.read_csv(DATA / "ScreenSequenceMap.csv", dtype={"ScreenID": str, "ModelID": str})
    x = seq[seq["ModelID"].isin(ids)].drop_duplicates(["ModelID", "ScreenID"])
    return {r.ModelID: {"screen": r.ScreenID, "library": r.Library, "days": int(r.Days), "screen_type": r.ScreenType,
                        "is_engineered": _any_true(pd.Series([r.IsEngineered])), "passes_qc": bool(r.PassesQC),
                        "in_gene_effect": bool(r.ScreenID in set(sm["ScreenID"]))} for r in x.itertuples()}


def purity(ct: pd.DataFrame) -> dict:
    """Per proxy group: size, majority library, fraction in it."""
    out = {}
    for g, row in ct.iterrows():
        out[str(g)] = {"n": int(row.sum()), "majority_library": str(row.idxmax()),
                       "majority_fraction": float(row.max() / row.sum()),
                       "libraries": {k: int(v) for k, v in row[row > 0].sort_values(ascending=False).items()}}
    return out


def main():
    sp = pd.read_csv(ROOT / "configs" / "splits.csv", index_col=0)
    sp.index = sp.index.astype(str)
    lg = pd.read_csv(ROOT / "configs" / "library_groups.csv", index_col=0)["lib_group"].astype(int)
    lib, per_model, counts = model_libraries()
    dev, test = sp.index[sp["role"] == "dev"], sp.index[sp["role"] == "test"]
    d = pd.DataFrame({"lib_group": lg.reindex(dev), "library": lib.reindex(dev)})
    missing_dev = d.index[d["library"].isna()].tolist()
    d = d.dropna()
    ct = pd.crosstab(d["lib_group"], d["library"])
    ct_lib = pd.crosstab(d["library"], d["lib_group"])
    # the predictability characterisation compares lineage means between proxy groups 0 and 1 ("cross-library")
    g01 = {str(g): purity(ct).get(str(g)) for g in (0, 1)}
    out = {
        "_note": ("POST HOC (after G7). Screen maps supplied by the user on 2026-10-28; no gene-effect values read, "
                  "no result re-run. Proxy = configs/library_groups.csv (dev row NaN-pattern groups)."),
        "inputs": {p.name: sha256(p) for p in [DATA / "CRISPRScreenMap.csv", DATA / "ScreenSequenceMap.csv",
                                                ROOT / "configs" / "library_groups.csv", ROOT / "configs" / "splits.csv"]},
        "screen_map": counts,
        "population": {
            "dev_lines": int(len(dev)), "dev_lines_mapped": int(len(d)), "dev_lines_unmapped": missing_dev,
            "test_lines": int(len(test)), "test_lines_mapped": int(lib.reindex(test).notna().sum()),
            "dev_library_counts": {k: int(v) for k, v in d["library"].value_counts().items()},
            "test_library_counts": {k: int(v) for k, v in lib.reindex(test).dropna().value_counts().items()},
            "dev_multi_library_lines": int(d["library"].str.contains(r"\+").sum()),
            "test_multi_library_lines": int(lib.reindex(test).dropna().str.contains(r"\+").sum()),
        },
        "proxy_vs_library": {
            "per_proxy_group": purity(ct),
            "per_library": {str(k): {str(g): int(v) for g, v in row[row > 0].items()} for k, row in ct_lib.iterrows()},
            "adjusted_rand": float(adjusted_rand_score(d["lib_group"], d["library"])),
            "nmi": float(normalized_mutual_info_score(d["lib_group"], d["library"])),
            "lines_in_majority_library_of_their_group": float(
                np.mean([ct.loc[g].idxmax() == l for g, l in zip(d["lib_group"], d["library"])])),
            "groups_0_1_used_for_crosslib_lineage_means": g01,
        },
        # test lines: aggregate only. G4 assigned each test line to a proxy group from its own NaN pattern
        # (results/test/b0/summary.json within_group/library_group_counts); per-line groups were not saved.
        "test_proxy_vs_library_aggregate": {
            "proxy_group_counts_g4": json.loads((ROOT / "results" / "test" / "b0" / "summary.json").read_text())
            ["within_group"]["library_group_counts"],
            "library_counts": {k: int(v) for k, v in lib.reindex(test).dropna().value_counts().items()},
            "dev_majority_library_of_group": {str(g): str(ct.loc[g].idxmax()) for g in range(4)},
        },
        # set A (G6): DepMap's RPE1-derived sublines, the non-transformed DepMap screens
        "set_A_rpe1_sublines": _set_a(),
        "screen_facts_dev": {
            "screen_types": {k: int(v) for k, v in per_model.reindex(d.index)["screen_types"].value_counts().items()},
            "median_days": {k: float(v) for k, v in
                            per_model.reindex(d.index).groupby(d["library"])["days"].median().items()},
            "lines_with_drug_treated_screen": int(per_model.reindex(d.index)["drug_treated"].sum()),
            "lines_with_engineered_screen": int(per_model.reindex(d.index)["engineered"].sum()),
        },
    }
    OUT.write_text(json.dumps(out, indent=1))
    print(json.dumps({k: out[k] for k in ("screen_map", "population")}, indent=1))
    print(json.dumps(out["proxy_vs_library"], indent=1)[:4000])
    print(json.dumps(out["screen_facts_dev"], indent=1))


if __name__ == "__main__":
    main()
