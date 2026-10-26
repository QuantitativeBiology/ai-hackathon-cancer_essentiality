"""Extract non-transformed-cell and anchor screens from the BioGRID ORCS human archive.

Input : data/external_screens/biogrid_orcs/BIOGRID-ORCS-ALL-homo_sapiens-<ver>.screens.tar.gz
        (downloaded by scripts/fetch_external_screens.py); curation lists in
        configs/external_screens.json -> "orcs".
Output: data/external_screens/biogrid_orcs/extracted/  (index + selected screen tab files; git-ignored)
        results/external_screens/orcs_nontransformed_screens.csv  (metadata only, small, committed)
        results/external_screens/orcs_anchor_screens.csv          (metadata + DepMap ModelID mapping)

ORCS is used only as a test-set / harmonisation source, never as a model input.
Archive members are read as data (written by basename only; nothing is executed).
"""
import argparse
import json
import pathlib
import re
import tarfile

import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
META = ["#SCREEN_ID", "AUTHOR", "SOURCE_ID", "CELL_LINE", "CELL_TYPE", "LIBRARY", "LIBRARY_TYPE",
        "LIBRARY_METHODOLOGY", "ENZYME", "PHENOTYPE", "SCREEN_TYPE", "EXPERIMENTAL_SETUP",
        "CONDITION_NAME", "DURATION", "MOI", "ANALYSIS", "SCORE.1_TYPE", "SCORE.2_TYPE",
        "SIGNIFICANCE_CRITERIA", "NUMBER_OF_HITS", "FULL_SIZE", "THROUGHPUT"]


def norm_name(s):
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


def depmap_map(lines, model_csv):
    m = pd.read_csv(model_csv, usecols=["ModelID", "CellLineName", "StrippedCellLineName"])
    key = {}
    for r in m.itertuples(index=False):
        for k in (norm_name(r.StrippedCellLineName), norm_name(r.CellLineName)):
            key.setdefault(k, r.ModelID)
    return {l: key.get(norm_name(l), "") for l in lines}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=str(ROOT / "configs/external_screens.json"))
    ap.add_argument("--orcs-dir", default=str(ROOT / "data/external_screens/biogrid_orcs"))
    ap.add_argument("--out", default=str(ROOT / "results/external_screens"))
    a = ap.parse_args()
    cfg = json.load(open(a.config, encoding="utf-8"))["orcs"]
    odir = pathlib.Path(a.orcs_dir)
    tgz = sorted(odir.glob("BIOGRID-ORCS-ALL-homo_sapiens-*.screens.tar.gz"))[-1]
    ext = odir / "extracted"
    ext.mkdir(exist_ok=True)
    idx_path = ext / cfg["index_member"]
    if not idx_path.exists():
        with tarfile.open(tgz, "r:gz") as t:
            idx_path.write_bytes(t.extractfile(cfg["index_member"]).read())
    idx = pd.read_csv(idx_path, sep="\t", low_memory=False)

    nt = cfg["nontransformed_cell_lines"]
    sel = idx[idx["CELL_LINE"].isin(nt)].copy()
    sel["category"] = sel["CELL_LINE"].map(nt)
    sel["fitness_like"] = (sel["PHENOTYPE"].isin(cfg["fitness_phenotypes"])
                           & (sel["EXPERIMENTAL_SETUP"] == "Timecourse"))
    sel["baseline_fitness"] = sel["fitness_like"] & (sel["CONDITION_NAME"] == "-")
    extra = idx[idx["#SCREEN_ID"].isin(cfg["extra_screen_ids"]) & ~idx["#SCREEN_ID"].isin(sel["#SCREEN_ID"])].copy()
    extra["category"] = "extra (see notes)"
    extra["fitness_like"] = False
    extra["baseline_fitness"] = False
    sel = pd.concat([sel, extra])
    anc = idx[idx["#SCREEN_ID"].isin(cfg["anchor_screen_ids"])].copy()
    anc["category"] = "anchor"
    anc["DepMapModelID"] = anc["CELL_LINE"].map(depmap_map(anc["CELL_LINE"].unique(), ROOT / "data/Model.csv"))

    want = set(sel.loc[sel["fitness_like"] | sel["#SCREEN_ID"].isin(cfg["extra_screen_ids"]), "#SCREEN_ID"]) | set(anc["#SCREEN_ID"])
    names = {cfg["screen_member_pattern"].format(id=i): i for i in want}
    missing = [n for n in names if not (ext / n).exists()]
    if missing:
        with tarfile.open(tgz, "r:gz") as t:
            for mem in t:
                base = pathlib.PurePosixPath(mem.name).name
                if base in missing and mem.isfile():
                    (ext / base).write_bytes(t.extractfile(mem).read())
    sel["extracted_file"] = sel["#SCREEN_ID"].map({v: k for k, v in names.items()}).fillna("")
    anc["extracted_file"] = anc["#SCREEN_ID"].map({v: k for k, v in names.items()}).fillna("")

    out = pathlib.Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    cols = META + ["category", "fitness_like", "baseline_fitness", "extracted_file"]
    sel[cols].sort_values(["category", "#SCREEN_ID"]).to_csv(out / "orcs_nontransformed_screens.csv", index=False)
    anc[META + ["category", "DepMapModelID", "extracted_file"]].sort_values("#SCREEN_ID").to_csv(
        out / "orcs_anchor_screens.csv", index=False)
    print(f"ORCS {tgz.name}: {len(idx)} screens; non-transformed rows {len(sel)} "
          f"(fitness-like {int(sel.fitness_like.sum())}, baseline {int(sel.baseline_fitness.sum())}); "
          f"anchors {len(anc)}; extracted {len(names)} screen files")


if __name__ == "__main__":
    main()
