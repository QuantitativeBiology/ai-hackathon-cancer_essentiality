"""Check (post hoc, 2026-10-28; dev lines only): predict.py with user-supplied genomics reproduces m5_avg's
registered dev-CV predictions.

    python scripts/predict_genomics_check.py      # -> results/predict_genomics_check.json (~30 min)

1. The held-out lines of CV fold 0 are written to scratch files in DepMap release layout (all their rows,
   including non-default entries) from data/: expression, OmicsCNGeneWGS, damaging, hotspot, Model.csv.
2. Blocks: predict.genomic_blocks on those files must equal the cached blocks the CV run used (C, M, H, meta;
   M2's gene-level CN and damaging), with the same NaN pattern.
3. End to end: predict.py --model m5_avg with --cn/--damaging/--hotspot/--metadata and --exclude-lines = the
   fold-0 held-out lines (so it trains on fold 0's training lines) must reproduce the registered out-of-fold
   predictions of the m5_avg CV run for those lines within CV_REPRO_TOL (pred.parquet is float32).
No test, external or CRISPR data are read; dev targets are read only by the fit, through the guard.
"""
from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
import predict as PD  # noqa: E402

M5_RUN = "20261027-161041_m5_avg_6dbebcf4"
CV_REPRO_TOL = 1e-3
OUT = ROOT / "results" / "predict_genomics_check.json"
FILES = {"expression": "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv", "cn": "OmicsCNGeneWGS.csv",
         "damaging": "OmicsSomaticMutationsMatrixDamaging.csv", "hotspot": "OmicsSomaticMutationsMatrixHotspot.csv"}


def subset_release_file(src: Path, dst: Path, ids: set):
    """All rows of ``ids`` (any entry), DepMap layout unchanged."""
    parts = [c[c["ModelID"].isin(ids)] for c in pd.read_csv(src, index_col=0, chunksize=300, low_memory=False)]
    pd.concat(parts).to_csv(dst)


def cached_blocks(ids):
    from essml.io import load_matrix
    from essml.models.g3_common import load_blocks
    b = {k: v.reindex(pd.Index(ids).intersection(v.index)) for k, v in load_blocks(["C", "M", "H", "meta"]).items()}
    for key, name in (("_cn", "cn"), ("_damaging", "damaging")):
        m = load_matrix(name)
        b[key] = m.loc[m.index.intersection(pd.Index(ids))]
    return b


def compare_blocks(built: dict, ref: dict) -> dict:
    out = {}
    for k, r in ref.items():
        b = built[k]
        rows = r.index.astype(str)
        b = b.reindex(index=rows)
        r = r.set_axis(rows, axis=0)
        if k == "meta":
            same = all((b[c].astype("string").fillna("<NA>") == r[c].astype("string").fillna("<NA>")).all()
                       for c in r.columns if c != "Age")
            age = np.allclose(b["Age"].astype(float), r["Age"].astype(float), equal_nan=True)
            out[k] = {"rows": int(len(r)), "equal": bool(same and age)}
            continue
        cols = r.columns.astype(str)
        B = b.reindex(columns=cols).to_numpy(np.float64)
        R = r.to_numpy(np.float64)
        nan_same = bool((np.isnan(B) == np.isnan(R)).all())
        d = np.abs(np.nan_to_num(B) - np.nan_to_num(R))
        out[k] = {"rows": int(len(r)), "cols": int(len(cols)), "same_nan_pattern": nan_same,
                  "max_abs_diff": float(d.max()) if d.size else 0.0,
                  "equal": bool(nan_same and (d.max() if d.size else 0.0) == 0.0)}
    return out


def main():
    sp = pd.read_csv(ROOT / "configs" / "splits.csv", index_col=0)
    sp.index = sp.index.astype(str)
    dev = sp[sp["role"] == "dev"]
    held = sorted(dev.index[dev["cv_fold"].astype(int) == 0])
    t0 = time.perf_counter()
    tmp = ROOT / "cache" / "predict_genomics_check"  # git-ignored; DepMap rows are never committed
    tmp.mkdir(parents=True, exist_ok=True)
    for k, f in FILES.items():
        if not (tmp / f).exists():
            subset_release_file(ROOT / "data" / f, tmp / f, set(held))
    model = pd.read_csv(ROOT / "data" / "Model.csv", dtype=str)
    model[model["ModelID"].isin(held)].to_csv(tmp / "Model.csv", index=False)
    expr = PD.read_expression(tmp / FILES["expression"])
    built = PD.genomic_blocks(expr.index, PD.read_genomic(tmp / FILES["cn"]),
                              PD.read_genomic(tmp / FILES["damaging"]), PD.read_genomic(tmp / FILES["hotspot"]),
                              PD.read_metadata(tmp / "Model.csv"))
    blocks = compare_blocks(built, cached_blocks(held))
    print(json.dumps(blocks, indent=1), flush=True)
    if not all(v["equal"] for v in blocks.values()):
        raise SystemExit("block check failed; end-to-end run skipped")
    out_csv = tmp / "pred.parquet"
    cmd = [sys.executable, str(ROOT / "scripts" / "predict.py"), "--expression", str(tmp / FILES["expression"]),
           "--cn", str(tmp / FILES["cn"]), "--damaging", str(tmp / FILES["damaging"]),
           "--hotspot", str(tmp / FILES["hotspot"]), "--metadata", str(tmp / "Model.csv"),
           "--exclude-lines", ",".join(held), "--no-ad", "--out", str(out_csv)]
    t1 = time.perf_counter()
    subprocess.run(cmd, check=True, cwd=ROOT)
    t_cli = time.perf_counter() - t1
    long = pd.read_parquet(out_csv)
    wide = long.pivot(index="sample", columns="entrez", values="predicted_effect")
    wide.columns = wide.columns.astype(str)
    reg = pd.read_parquet(ROOT / "cache" / "runs" / M5_RUN / "pred.parquet")
    reg.columns, reg.index = reg.columns.astype(str), reg.index.astype(str)
    lines = pd.Index(held).intersection(wide.index)
    g = wide.columns.intersection(reg.columns)
    d = np.abs(wide.loc[lines, g].to_numpy(np.float64) - reg.loc[lines, g].to_numpy(np.float64))
    e2e = {"n_lines": int(len(lines)), "n_genes": int(len(g)), "n_genes_registered": int(reg.shape[1]),
           "max_abs_diff_vs_registered_oof": float(d.max()), "mean_abs_diff": float(d.mean()),
           "tol": CV_REPRO_TOL, "cli_s": round(t_cli, 1)}
    e2e["ok"] = bool(len(lines) == len(held) and len(g) == reg.shape[1] and e2e["max_abs_diff_vs_registered_oof"] <= CV_REPRO_TOL)
    res = {"_note": ("POST HOC (2026-10-28, after G7), dev lines only: predict.py with user-supplied genomics in "
                     "DepMap release layout vs the registered m5_avg dev-CV out-of-fold predictions (fold 0). "
                     "A reproduction check, not a new evaluation: when it passes, the script's genomics path is the "
                     "evaluated headline (dev-CV P1 0.3175; test 0.319)."),
           "cv_run": M5_RUN, "fold": 0, "held_out_lines": int(len(held)), "blocks": blocks,
           "blocks_ok": bool(all(v["equal"] for v in blocks.values())), "end_to_end": e2e,
           "command": ("python scripts/predict.py --expression <expr> --cn <OmicsCNGeneWGS rows> --damaging <rows> "
                       "--hotspot <rows> --metadata <Model.csv rows> --exclude-lines <fold-0 held-out lines> --no-ad "
                       "--out <scratch>"),
           "total_s": round(time.perf_counter() - t0, 1)}
    res["ok"] = bool(res["blocks_ok"] and e2e["ok"])
    OUT.write_text(json.dumps(res, indent=1))
    print(json.dumps(e2e, indent=1), "\nok:", res["ok"])


if __name__ == "__main__":
    main()
