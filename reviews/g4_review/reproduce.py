"""Independent reproduction of the G4 test-set numbers (DESIGN §9 "independent reproduction").

Written from DESIGN.md §4 definitions by the G4 reviewer. Uses only pandas / numpy / scipy; nothing is
imported from essml. Test targets are read directly from the raw DepMap file data/CRISPRGeneEffect.csv
(post-unblinding reproduction read; logged in results/test_access_log.csv before this script ran).

    python reviews/g4_review/reproduce.py            # writes reviews/g4_review/reproduce_output.json

Outputs aggregate metrics only (no target values, no per-line or per-gene values).

Definitions (DESIGN §4; evaluate.py consulted only where DESIGN leaves details open):
- Lines: the 255 role=='test' lines of configs/splits.csv. Genes: configs/eval_genes.csv.
- A gene is scored iff it has >= 50 observed (non-NaN) test lines (min_obs, DESIGN §2/§4).
- Classes: configs/gene_classes.csv (dev-wide, frozen), genes missing from it are in no class.
- P1(class): mean over scored class genes of the Pearson r across observed test lines. A prediction or
  truth that is exactly constant over the observed lines gives r = 0.
- dP1 vs B: mean over the same genes of r - r_B.
- P2a: per line, AUPRC over the scored 'all' genes (label truth < -0.5, score -prediction), step-wise
  average precision (sum over distinct score thresholds of (R_k - R_{k-1}) * P_k); lines with no positive or
  no negative are skipped; mean over lines. dP2a vs B0 = mean over lines of AP - AP_B0.
- P2b: per gene over observed lines, label truth < -0.5; genes with >= 3 positives and >= 3 negatives;
  nAP = (AP - prev) / (1 - prev); AUROC = Mann-Whitney with average ranks. Headline = selective genes with
  prevalence <= 0.5.
- Tiers: configs/predictable_genes.csv (primary). Per tier: P1 over scored genes of the tier (all classes)
  and over its variable genes, P2b over selective & tier genes (prev <= 0.5).
- Spearman(dev mean_r, test r) over scored, tiered genes (and the variable subset).
- Bootstrap (own scheme): patients resampled with replacement within lineage strata (stratum = the
  patient's modal lineage, ties -> alphabetical first), 2000 draws, seed 777; a resampled line enters with
  weight = number of draws of its patient (weighted Pearson); gene sets fixed at the point estimate;
  percentile 95% CI. All models share the draws, so differences are paired.
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr, chi2_contingency

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent / "reproduce_output.json"
MIN_OBS = 50
THR = -0.5
N_BOOT = 2000
SEED = 777
CLASSES = ["all", "variable", "selective", "rare_selective"]
TIERS = ["high", "medium", "low"]


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


# ------------------------------------------------------------------ metric primitives (own code)
def pearson_cols(P, T, min_obs=MIN_OBS):
    """Per-column Pearson r over rows where T is observed. Exact-constant column -> 0; n < min_obs -> NaN."""
    obs = ~np.isnan(T)
    n = obs.sum(0)
    out = np.full(P.shape[1], np.nan)
    for g in range(P.shape[1]):
        o = obs[:, g]
        if n[g] < max(min_obs, 2):
            continue
        p, t = P[o, g], T[o, g]
        if p.max() == p.min() or t.max() == t.min():
            out[g] = 0.0
            continue
        p = p - p.mean()
        t = t - t.mean()
        out[g] = float(np.dot(p, t) / np.sqrt(np.dot(p, p) * np.dot(t, t)))
    return out


def average_precision(y, s):
    """Step-wise AP (higher score = predicted positive); ties form one threshold."""
    y = np.asarray(y, bool)
    order = np.argsort(-s, kind="mergesort")
    s_sorted, y_sorted = s[order], y[order]
    last = np.r_[np.flatnonzero(np.diff(s_sorted) != 0), len(s_sorted) - 1]
    tp = np.cumsum(y_sorted)[last].astype(float)
    k = (last + 1).astype(float)
    prec = tp / k
    rec = tp / tp[-1]
    rec_prev = np.r_[0.0, rec[:-1]]
    return float(np.sum((rec - rec_prev) * prec))


def auroc(y, s):
    y = np.asarray(y, bool)
    r = rankdata(s)
    npos, nneg = y.sum(), (~y).sum()
    return float((r[y].sum() - npos * (npos + 1) / 2) / (npos * nneg))


def per_gene_p2b(P, T, genes_mask):
    """nAP, AUROC, prevalence per gene (NaN where not eligible)."""
    G = P.shape[1]
    nap, auc, prev = np.full(G, np.nan), np.full(G, np.nan), np.full(G, np.nan)
    obs = ~np.isnan(T)
    for g in np.flatnonzero(genes_mask):
        o = obs[:, g]
        y = T[o, g] < THR
        npos, nneg = int(y.sum()), int((~y).sum())
        prev[g] = npos / o.sum()
        if npos >= 3 and nneg >= 3:
            s = -P[o, g]
            ap = average_precision(y, s)
            nap[g] = (ap - prev[g]) / (1 - prev[g])
            auc[g] = auroc(y, s)
    return nap, auc, prev


def per_line_ap(P, T, gmask):
    out = np.full(P.shape[0], np.nan)
    for i in range(P.shape[0]):
        o = ~np.isnan(T[i]) & gmask
        y = T[i, o] < THR
        if y.any() and not y.all():
            out[i] = average_precision(y, -P[i, o])
    return out


def weighted_pearson(W, P, T):
    """W (B x n) line weights; P, T (n x G), T NaN = unobserved. Returns (B x G) r. Exact-constant -> 0."""
    O = (~np.isnan(T)).astype(np.float64)
    T0 = np.where(O > 0, T, 0.0)
    # centre each gene on its unweighted observed mean (shift-invariant; constant columns become exactly 0)
    n = O.sum(0)
    Pc = (P - (P * O).sum(0) / n) * O
    Tc = (T0 - T0.sum(0) / n) * O
    Sw = W @ O
    Sp, St = W @ Pc, W @ Tc
    Spp, Stt, Spt = W @ (Pc * Pc), W @ (Tc * Tc), W @ (Pc * Tc)
    with np.errstate(invalid="ignore", divide="ignore"):
        cov = Spt - Sp * St / Sw
        vp = Spp - Sp * Sp / Sw
        vt = Stt - St * St / Sw
        r = cov / np.sqrt(vp * vt)
    const = (np.abs(Pc).sum(0) == 0) | (np.abs(Tc).sum(0) == 0)
    r[:, const] = 0.0
    r[(vp <= 0) | (vt <= 0)] = 0.0
    return r


def patient_weights(lineage, patient, n_boot, seed):
    """(n_boot x n_lines) weights: patients resampled with replacement within modal-lineage strata."""
    df = pd.DataFrame({"lin": lineage.to_numpy(), "pat": patient.to_numpy()})
    modal = df.groupby("pat")["lin"].agg(lambda s: sorted(s.value_counts()[lambda v: v == v.max()].index)[0])
    pats = modal.index.to_numpy()
    pat_idx = {p: i for i, p in enumerate(pats)}
    line_pat = np.array([pat_idx[p] for p in df["pat"]])
    rng = np.random.default_rng(seed)
    strata = [np.flatnonzero(modal.to_numpy() == s) for s in sorted(modal.unique())]
    W = np.zeros((n_boot, len(df)))
    for b in range(n_boot):
        cnt = np.zeros(len(pats))
        for st in strata:
            np.add.at(cnt, rng.choice(st, size=len(st), replace=True), 1)
        W[b] = cnt[line_pat]
    return W


def ci(point, reps):
    reps = reps[np.isfinite(reps)]
    return {"point": float(point), "ci_lo": float(np.percentile(reps, 2.5)), "ci_hi": float(np.percentile(reps, 97.5)),
            "se": float(reps.std(ddof=1)), "n": int(len(reps))}


# ------------------------------------------------------------------ main
def main():
    t0 = time.time()
    res = {"script": "reviews/g4_review/reproduce.py", "definitions": __doc__.split("Definitions", 1)[1].strip()}
    sp = pd.read_csv(ROOT / "configs/splits.csv", index_col=0)
    sp.index = sp.index.astype(str)
    test_ids = sp.index[sp["role"] == "test"]
    eval_genes = pd.read_csv(ROOT / "configs/eval_genes.csv", dtype=str)["entrez"].tolist()
    gc = pd.read_csv(ROOT / "configs/gene_classes.csv", dtype={"entrez": str}).set_index("entrez")
    cls = {c: gc[c].reindex(eval_genes).fillna(False).astype(bool).to_numpy() for c in CLASSES}
    tiers = {"primary": pd.read_csv(ROOT / "configs/predictable_genes.csv", dtype={"entrez": str}).set_index("entrez"),
             "m1_own": pd.read_csv(ROOT / "configs/predictable_genes_m1.csv", dtype={"entrez": str}).set_index("entrez")}
    fin = json.loads((ROOT / "configs/finalists.json").read_text())
    names = [fin["headline"]["name"]] + [e["name"] for e in fin["additional"]] + [fin["T"]["name"],
             fin["H3_noEC"]["name"]] + [fin["baselines"][b]["name"] for b in ("B0", "B1", "B2")]
    base_name = {b: fin["baselines"][b]["name"] for b in ("B0", "B1", "B2")}
    headline, additional = fin["headline"]["name"], fin["additional"][0]["name"]

    # 1. prediction files vs the committed manifest
    man = json.loads((ROOT / "results/test/predictions_manifest.json").read_text())
    shas = {}
    for n, e in man.items():
        files = [(e["file"], e["sha256"])]
        if "predictions_with_confidence" in e:
            files.append((e["predictions_with_confidence"]["file"], e["predictions_with_confidence"]["sha256"]))
        files += [(v["file"], v["sha256"]) for v in (e.get("fold_preds") or {}).values()]
        for f, s in files:
            shas[f] = sha256(ROOT / f) == s
    res["sha256_match"] = shas
    res["sha256_all_match"] = all(shas.values())

    # 2. raw targets (direct CSV read, float64)
    hdr = pd.read_csv(ROOT / "data/CRISPRGeneEffect.csv", nrows=0).columns
    col_ent = {}
    for c in hdr[1:]:
        m = re.match(r"^(.*) \((\d+)\)$", c)
        if m and m.group(2) not in col_ent.values():
            col_ent[c] = m.group(2)
    raw = pd.read_csv(ROOT / "data/CRISPRGeneEffect.csv", index_col=0, usecols=[hdr[0], *col_ent])
    raw.index = raw.index.astype(str)
    raw.columns = [col_ent[c] for c in raw.columns]
    res["raw_shape"] = list(raw.shape)
    res["test_ids_missing_from_raw"] = sorted(set(test_ids) - set(raw.index))
    res["eval_genes_missing_from_raw"] = len(set(eval_genes) - set(raw.columns))
    truth = raw.loc[test_ids, eval_genes].to_numpy(np.float64)
    dev_ids = sp.index[sp["role"] == "dev"]
    dev_nan = raw.loc[dev_ids].isna()
    test_nan = raw.loc[test_ids].isna()
    del raw
    n_obs = (~np.isnan(truth)).sum(0)
    evalg = n_obs >= MIN_OBS
    res["n_eval_genes"] = len(eval_genes)
    res["n_eval_genes_ge50_obs_in_test"] = int(evalg.sum())
    res["n_lines"] = int(truth.shape[0])
    truth32 = truth.astype(np.float32).astype(np.float64)
    res["n_entries_label_flip_float32"] = int(np.sum((truth < THR) != (truth32 < THR)))

    # 3. predictions
    preds = {}
    for n in names:
        p = pd.read_parquet(ROOT / man[n]["file"])
        p.columns = p.columns.astype(str)
        p.index = p.index.astype(str)
        P = p.loc[test_ids, eval_genes].to_numpy(np.float64)
        if not np.isfinite(P[~np.isnan(truth)]).all():
            raise ValueError(f"{n}: non-finite prediction at observed entries")
        preds[n] = P

    # 4. point metrics
    rs, p2b, ap_line = {}, {}, {}
    for n, P in preds.items():
        rs[n] = pearson_cols(P, truth)
        p2b[n] = per_gene_p2b(P, truth, cls["selective"] & evalg)
        ap_line[n] = per_line_ap(P, truth, cls["all"] & evalg)

    # per-gene cross-check against the project's per_gene tables (only the max |difference| is kept)
    res["per_gene_max_abs_diff_vs_project"] = {}
    for n in names:
        pg = pd.read_csv(ROOT / "results/test" / n / "per_gene.csv.gz", dtype={"gene": str}).set_index("gene")
        mine = pd.DataFrame({"r": rs[n], "nap": p2b[n][0], "auroc": p2b[n][1]}, index=eval_genes)
        res["per_gene_max_abs_diff_vs_project"][n] = {
            k: float(np.nanmax(np.abs(mine[k] - pg[k].reindex(eval_genes)).to_numpy())) for k in mine}
        res["per_gene_max_abs_diff_vs_project"][n]["n_obs_identical"] = bool(
            (pg["n_obs"].reindex(eval_genes).to_numpy() == n_obs).all())

    out = {}
    for n, P in preds.items():
        r = rs[n]
        d = {}
        for c in CLASSES:
            m = cls[c] & evalg
            d[f"P1_{c}"] = float(np.nanmean(r[m]))
            d[f"n_genes_{c}"] = int(m.sum())
        m = cls["variable"] & evalg
        d["P1_median_variable"] = float(np.nanmedian(r[m]))
        d["n_r_gt_0.3_variable"] = int((r[m] > 0.3).sum())
        for b, bn in base_name.items():
            if bn == n:
                continue
            d[f"dP1_vs_{b}_variable"] = float(np.nanmean((r - rs[bn])[m]))
            d[f"frac_beat_{b}_variable"] = float(np.mean(r[m] > rs[bn][m]))
        d["P2a"] = float(np.nanmean(ap_line[n]))
        if n != base_name["B0"]:
            d["dP2a_vs_B0"] = float(np.nanmean(ap_line[n] - ap_line[base_name["B0"]]))
        nap, auc, prev = p2b[n]
        e05 = cls["selective"] & evalg & np.isfinite(nap) & (prev <= 0.5)
        d["P2b_nap_prev05"] = float(np.mean(nap[e05]))
        d["P2b_auroc_prev05"] = float(np.mean(auc[e05]))
        d["P2b_n_prev05"] = int(e05.sum())
        for tname, tab in tiers.items():
            if tname == "m1_own" and n != additional:
                continue
            tier = tab["tier"].reindex(eval_genes).fillna("not_evaluable").to_numpy()
            tt = {}
            for t in TIERS:
                mt = (tier == t) & evalg
                mv = mt & cls["variable"]
                e = cls["selective"] & mt & np.isfinite(nap) & (prev <= 0.5)
                tt[t] = {"n_genes": int(mt.sum()), "P1": float(np.nanmean(r[mt])),
                         "P1_variable": float(np.nanmean(r[mv])), "n_genes_variable": int(mv.sum()),
                         "P2b_nap_prev05": float(np.mean(nap[e])) if e.any() else None,
                         "P2b_auroc_prev05": float(np.mean(auc[e])) if e.any() else None, "P2b_n_prev05": int(e.sum())}
            mean_r = tab["mean_r"].reindex(eval_genes).to_numpy(np.float64)
            tiered = np.isin(tier, TIERS) & evalg
            sp_ = {}
            for lab, mm in (("evaluable", tiered), ("variable", tiered & cls["variable"])):
                ok = mm & np.isfinite(r) & np.isfinite(mean_r)
                sp_[lab] = (None if np.ptp(r[ok]) == 0 else float(spearmanr(r[ok], mean_r[ok])[0]))
                sp_[lab + "_n"] = int(ok.sum())
            tt["spearman_dev_mean_r_vs_test_r"] = sp_
            d["tiers_" + tname] = tt
        out[n] = d
    res["point"] = out

    # 5. bootstrap: P1 variable, dP1 vs B1/B2, and the paired headline - additional difference
    lineage, patient = sp.loc[test_ids, "OncotreeLineage"].astype(str), sp.loc[test_ids, "PatientID"].astype(str)
    W = patient_weights(lineage, patient, N_BOOT, SEED)
    mv = cls["variable"] & evalg
    Tm = truth[:, mv]
    R = {n: np.nanmean(weighted_pearson(W, preds[n][:, mv], Tm), axis=1) for n in names}
    boot = {}
    for n in names:
        b = {"P1_variable": ci(out[n]["P1_variable"], R[n])}
        for bb in ("B1", "B2"):
            bn = base_name[bb]
            if bn != n:
                b[f"dP1_vs_{bb}_variable"] = ci(out[n][f"dP1_vs_{bb}_variable"], R[n] - R[bn])
        boot[n] = b
    diff_point = out[headline]["P1_variable"] - out[additional]["P1_variable"]
    dd = R[headline] - R[additional]
    paired = ci(diff_point, dd)
    paired["frac_boot_le_0"] = float(np.mean(dd <= 0))
    paired["frac_genes_headline_better"] = float(np.mean(rs[headline][mv] > rs[additional][mv]))
    res["bootstrap"] = {"scheme": "patients within modal-lineage strata, weights = draw counts, seed 777, 2000 draws, "
                                  "percentile CI, gene set fixed", "n_patients": int(patient.nunique()),
                        "n_strata": int(lineage.nunique()), "per_finalist": boot,
                        f"paired_P1_variable_{headline}_minus_{additional}": paired}

    # 6. population diagnostics (lineage mix, NaN-pattern library groups)
    lin_dev = sp.loc[dev_ids, "OncotreeLineage"].value_counts()
    lin_test = sp.loc[test_ids, "OncotreeLineage"].value_counts()
    tab = pd.concat([lin_dev, lin_test], axis=1).fillna(0)
    chi2, pval, dof, _ = chi2_contingency(tab.to_numpy().T)
    res["lineage_mix"] = {"n_lineages_dev": int(len(lin_dev)), "n_lineages_test": int(len(lin_test)),
                          "lineages_test_not_in_dev": sorted(set(lin_test.index) - set(lin_dev.index)),
                          "chi2": float(chi2), "dof": int(dof), "p": float(pval),
                          "max_abs_share_diff": float((tab.iloc[:, 0] / tab.iloc[:, 0].sum()
                                                       - tab.iloc[:, 1] / tab.iloc[:, 1].sum()).abs().max())}
    keys_dev = pd.Series([m.tobytes() for m in dev_nan.to_numpy()], index=dev_ids)
    top = list(keys_dev.value_counts().index[:4])
    keys_test = [m.tobytes() for m in test_nan.to_numpy()]
    pos = {k: i for i, k in enumerate(top)}
    lib_test = pd.Series([pos.get(k, 4) for k in keys_test]).value_counts().sort_index()
    lib_dev = keys_dev.map(pos).fillna(4).astype(int).value_counts().sort_index()
    res["library_groups"] = {"dev": {str(k): int(v) for k, v in lib_dev.items()},
                             "test": {str(k): int(v) for k, v in lib_test.items()}}
    chi2, pval, dof, _ = chi2_contingency(np.vstack([lib_dev.reindex(range(5)).fillna(0), lib_test.reindex(range(5)).fillna(0)]))
    res["library_groups"].update(chi2=float(chi2), dof=int(dof), p=float(pval))
    res["runtime_s"] = round(time.time() - t0, 1)
    OUT.write_text(json.dumps(res, indent=1))
    print(f"wrote {OUT} in {res['runtime_s']} s")


if __name__ == "__main__":
    sys.exit(main())
