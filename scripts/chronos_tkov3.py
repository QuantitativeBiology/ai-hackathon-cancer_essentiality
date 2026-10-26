"""Chronos refit for the TKOv3 library (DESIGN v3 §7.2 step 1; PI decision on the G5a Chronos item).

    python scripts/chronos_tkov3.py [--epochs 301]

Why: the only TKOv3 anchor (HAP1, Mair 2019) fails same-line NN recall under median guide LFC, which
is the DESIGN trigger for a fresh per-library Chronos fit.

Design of the fit (stated before running):
  * one JOINT fit over all TKOv3 count screens with a usable time axis: the HAP1 anchor (T0 + T3..T18)
    and the set-B TKOv3 screens (Olivieri NT arms, Brown PACE002 T18/T30, PACE006 T25, Mair H1 MEF
    T4/T8 and laminin T4/T9/T12). Chronos shares guide efficacy across screens of a library, and a
    single screen identifies it poorly, which is why the fit is joint. Brown Durocher 'T3' is left out
    (time unit not stated in the deposit) and keeps median-LFC scoring.
  * no pDNA was sequenced for the set-B screens: each screen's own T0 cell population is its
    'pDNA' reference (pDNA_batch = screen), days counted from T0.
  * guides with < 30 reads in a batch's reference are set to NaN for that batch (same filter as the
    median-LFC route); LacZ / luciferase / EGFP guides are negative controls.
  * outputs: raw Chronos gene effects are affine-scaled per screen exactly like the median-LFC route
    (DepMap-matched targets, configs/harmonisation_frozen.json / control_split.json).
Blinding: set-B gene effects are written to cache/harmonised/ only (access-log rows); results/ get the
HAP1 anchor comparison and control-gene (evaluation-half) SSMD / NNMD only.
Pre-stated adoption rule (PI): adopt Chronos for TKOv3 only if the HAP1 deviation correlation with
DepMap HAP1 improves by >= +0.02 over median LFC. HAP1 is the single TKOv3 anchor: weak evidence.
"""
from __future__ import annotations

import argparse
import gzip
import io
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from essml import guard  # noqa: E402
from essml import harmonise as H  # noqa: E402
import harmonise_screens as HS  # noqa: E402
HS.CHRONOS_ADOPTED = False  # the comparator must be median guide LFC, never a previous Chronos result

DEP = ROOT / "data" / "external_screens"
OUT = ROOT / "results" / "harmonisation"
MIN_REF = 30


def _read_gz(p):
    raw = gzip.open(p).read().decode("latin1").replace("\r\n", "\n").replace("\r", "\n")
    return pd.read_csv(io.StringIO(raw), sep="\t")


def screens():
    """(screen id, role, table, sequence column, T0 column, {column: days})."""
    lib = _read_gz(DEP / "brown2019_rpe1_gse128210" / "GSE128210_TKOv3-Human-Library.txt.gz")
    gid2seq = dict(zip(lib["GUIDE_ID"], lib["SEQUENCE"]))
    out = []
    mair = pd.read_csv(DEP / "mair2019_h1_hesc" / "GSE156284_bmair_hESC_readCounts_matrix_12Aug20.txt.gz", sep="\t")
    hap = {"Moffat_411": 3, "Moffat_412": 3, "Moffat_413": 3, "Moffat_408": 6, "Moffat_409": 6, "Moffat_410": 6,
           "Moffat_405": 10, "Moffat_406": 10, "Moffat_407": 10, "Moffat_402": 14, "Moffat_403": 14,
           "Moffat_404": 14, "Moffat_1014": 18, "Moffat_1015": 18, "Moffat_1016": 18}
    out.append(("A_mair2019_HAP1", "anchor", mair, "SEQUENCE", "Moffat_343_814", hap))
    out.append(("mair2019_H1_MEF", "setB", mair, "SEQUENCE", "ESC_1",
                {"ESC_5": 4, "ESC_6": 4, "ESC_7": 4, "ESC_8": 8, "ESC_9": 8, "ESC_10": 8}))
    out.append(("mair2019_H1_laminin", "setB", mair, "SEQUENCE", "ESC_18",
                {**{f"ESC_{i}": 4 for i in (19, 20, 21)}, **{f"ESC_{i}": 9 for i in (22, 23, 24)},
                 **{f"ESC_{i}": 12 for i in (25, 26, 27)}}))
    oli = pd.read_csv(DEP / "olivieri2020_rpe1_tp53ko" / "Dataset_S2_readcounts.txt", sep="\t")
    oli["SEQUENCE"] = oli["sgRNA"].map(gid2seq)
    oli = oli.rename(columns={"Gene": "GENE"})
    for s, d in (("S02", 15), ("S04", 18), ("S05", 18), ("S06", 18), ("S07", 18), ("S12", 18)):
        tp = f"T{d}"
        out.append((f"olivieri2020_TKOv3_{s}", "setB", oli, "SEQUENCE", f"{s}_T0",
                    {f"{s}_NT_{tp}_A": d, f"{s}_NT_{tp}_B": d}))
    p2 = _read_gz(DEP / "brown2019_rpe1_gse128210" / "GSE128210_PACE002_RPE1-hTERT_WT-TKOv3_readcounts.txt.gz")
    out.append(("brown2019_WT_PACE002", "setB", p2, "SEQUENCE", "T0",
                {**{f"RPE1-hTERT_WT_T18{r}": 18 for r in "ABC"}, **{f"RPE1-hTERT_WT_T30{r}": 30 for r in "ABC"}}))
    p6 = _read_gz(DEP / "brown2019_rpe1_gse128210" / "GSE128210_PACE006_RPE1_WT-TKOv3_readcounts.txt.gz")
    out.append(("brown2019_WT_PACE006", "setB", p6, "SEQUENCE", "T0", {"RPE1_WT_T25A": 25, "RPE1_WT_T25B": 25}))
    return out, lib


def build_inputs(scr, lib, symmap):
    seqs = lib["SEQUENCE"]
    rows, smap, nfilt = [], [], {}
    for sid, role, tab, seqcol, t0, ends in scr:
        t = tab.dropna(subset=[seqcol]).drop_duplicates(seqcol).set_index(seqcol)
        cols = [t0] + list(ends)
        x = t[cols].apply(pd.to_numeric, errors="coerce").reindex(seqs)
        low = x[t0] < MIN_REF  # mask late timepoints only: Chronos refuses NaN in pDNA rows
        x.loc[low, list(ends)] = np.nan
        nfilt[sid] = int(low.sum())
        for c in cols:
            name = f"{sid}__{c}"
            rows.append(x[c].rename(name))
            smap.append({"sequence_ID": name, "cell_line_name": "pDNA" if c == t0 else sid,
                         "pDNA_batch": sid, "days": 0 if c == t0 else ends[c]})
    rc = pd.DataFrame(rows)
    ref_rows = [r["sequence_ID"] for r in smap if r["cell_line_name"] == "pDNA"]
    rc = rc.loc[:, rc.loc[ref_rows].notna().all(axis=0)]  # guide absent from a screen's table -> dropped
    gm = pd.DataFrame({"sgrna": lib["SEQUENCE"], "symbol": lib["GENE"]})
    gm["gene"] = gm["symbol"].map(symmap)  # LacZ / luciferase / EGFP (non-cutting) map to NaN and are dropped
    gm = gm.dropna(subset=["gene"])
    gm = gm[gm["sgrna"].isin(rc.columns)][["sgrna", "gene"]].drop_duplicates("sgrna")
    import chronos
    gm = chronos.filter_guides(gm).reset_index(drop=True)
    # Chronos wants CUTTING negative controls: guides against the NEGv1 SCALING half (evaluation half held out)
    neg_scaling = set(json.loads(HS.CONTROL_SPLIT.read_text())["neg_scaling"])
    negs = gm.loc[gm["gene"].isin(neg_scaling), "sgrna"].tolist()
    return rc[gm["sgrna"].tolist()], gm, pd.DataFrame(smap), negs, nfilt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--epochs", type=int, default=301)
    a = ap.parse_args()
    t0 = time.time()
    np.random.seed(0)
    import tensorflow as tf
    tf.compat.v1.set_random_seed(0)
    import chronos
    symmap = HS.SymbolMap()
    scr, lib = screens()
    for sid, role, *_ in scr:
        if role == "setB":
            guard._log("harmonise_write_only", pd.Series(["external"], index=[sid]), True,
                       caller="chronos_tkov3.py:main", scope=f"setB:{sid}",
                       reason="Chronos joint TKOv3 fit; set-B gene effects written to cache/harmonised only")
    rc, gm, sm, negs, nfilt = build_inputs(scr, lib, symmap)
    print(f"readcounts {rc.shape}, guides mapped {len(gm)}, genes {gm['gene'].nunique()}, neg ctrl {len(negs)}",
          flush=True)
    model = chronos.Chronos(readcounts={"TKOv3": rc}, guide_gene_map={"TKOv3": gm},
                            sequence_map={"TKOv3": sm}, negative_control_sgrnas={"TKOv3": negs})
    model.train(nepochs=a.epochs)
    ge = model.gene_effect.copy()
    ge.columns = ge.columns.astype(str)
    fit_s = time.time() - t0
    print(f"Chronos fit done in {fit_s:.0f} s; gene effect {ge.shape}", flush=True)

    # ---- affine scaling (same targets as the main pipeline) and comparison
    cs = json.loads(HS.CONTROL_SPLIT.read_text())
    sp = pd.read_csv(ROOT / "configs" / "splits.csv", index_col=0)
    Ydev = guard.load_targets(sp.index[sp["role"] == "dev"].tolist()).astype(np.float64)
    mdev = Ydev.mean()
    tS = float(Ydev[Ydev.columns.intersection(cs["scaling_set"])].median(axis=1).median())
    tN = float(Ydev[Ydev.columns.intersection(cs["neg_scaling"])].median(axis=1).median())
    ctrl = set(cs["scaling_set"]) | set(cs["neg_scaling"]) | set(cs["ceg_evaluation"]) | set(cs["neg_evaluation"])
    HS.CACHE.mkdir(parents=True, exist_ok=True)
    rows = []
    hcfg = json.loads(HS.CFG_HARM.read_text(encoding="utf-8"))
    by_id = {e["id"]: e for e in hcfg["entries"]}
    med_ids = {"A_mair2019_HAP1": "A_mair2019_HAP1_T18", "mair2019_H1_MEF": "mair2019_H1_MEF_T8",
               "mair2019_H1_laminin": "mair2019_H1_laminin_T12", "brown2019_WT_PACE002": "brown2019_WT_PACE002_T18",
               "brown2019_WT_PACE006": "brown2019_WT_PACE006_T25",
               **{f"olivieri2020_TKOv3_{s}": f"olivieri2020_TKOv3_{s}" for s in ("S02", "S04", "S05", "S06", "S07", "S12")}}
    for sid, role, *_ in scr:
        g = ge.loc[sid].dropna()
        ms, mn = H.affine_params(g, cs["scaling_set"], cs["neg_scaling"])
        gc = H.affine_apply(g, ms, mn, tS, tN)
        # median-LFC comparator through the same route (raw scores of the matching arm)
        e = by_id[med_ids[sid]]
        if role == "anchor":
            med, _ = HS.raw_arm_scores(e, symmap)
        else:
            med, _ = HS.guarded_controls(e, symmap)  # controls only for set B
        mm = H.affine_apply(med, *H.affine_params(med, cs["scaling_set"], cs["neg_scaling"]), tS, tN)
        for method, x in (("chronos", gc), ("median_lfc", mm)):
            xq = x if role == "anchor" else x[x.index.isin(ctrl)]
            q = H.qc_metrics(xq, cs["ceg_evaluation"], cs["neg_evaluation"])
            r = {"screen": sid, "role": role, "method": method, "n_ref_filtered_guides": nfilt[sid],
                 **{k: q[k] for k in ("n_genes", "median_ceg_eval", "median_neg_eval", "nnmd", "ssmd", "auroc")}}
            if role == "anchor":
                y = Ydev.loc["ACH-002475"]
                r.update(dev_corr=H.deviation_corr(x, y, mdev),
                         dev_corr_ctrlfree=H.deviation_corr(x, y, mdev, genes=set(Ydev.columns) - ctrl),
                         **H.nn_rank(x, Ydev, mdev, "ACH-002475"))
            rows.append(r)
        df = pd.DataFrame({"chronos_raw": g, "chronos_affine": gc})
        df.index.name = "entrez"
        df.astype("float32").to_parquet(HS.CACHE / f"chronos_TKOv3__{sid}.parquet")
    res = pd.DataFrame(rows)
    res.to_csv(OUT / "chronos_tkov3_comparison.csv", index=False)
    an = res[res["role"] == "anchor"].set_index("method")
    gain = float(an.loc["chronos", "dev_corr"] - an.loc["median_lfc", "dev_corr"])
    adopt = gain >= 0.02
    info = {"epochs": a.epochs, "fit_seconds": round(fit_s, 1), "readcounts_shape": list(rc.shape),
            "n_guides_mapped": int(len(gm)), "n_negative_control_guides": len(negs),
            "screens": [s[0] for s in scr], "excluded": {"brown2019_WT_Durocher_T3": "time unit not stated"},
            "affine_targets": [tS, tN], "hap1_dev_corr_gain": gain, "adoption_margin": 0.02,
            "adopt_chronos_for_TKOv3": bool(adopt),
            "chronos_version": __import__("importlib.metadata").metadata.version("crispr_chronos"),
            "tensorflow_version": tf.__version__}
    (OUT / "chronos_tkov3_decision.json").write_text(json.dumps(info, indent=1))
    pd.set_option("display.width", 220)
    print(res.round(3).to_string())
    print(json.dumps(info, indent=1))


if __name__ == "__main__":
    main()
