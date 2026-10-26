"""Headline evaluation (FROZEN: hashed in configs/frozen_hashes.json; DESIGN.md §4).

``evaluate`` computes every metric WITHIN each held-out fold and then averages over folds
(mean and fold SD, ddof=1). Out-of-fold predictions are never pooled across folds: pooling gives
B0 a spurious negative per-gene r (tested in tests/test_evaluate.py). With ``fold_of_line=None``
the whole set is one fold (the test-set case).

Metrics per gene class c. For reporting, classes are ONE dev-wide table (configs/gene_classes.csv,
identical for every model; DESIGN v3 §3); a dict of per-fold tables is also accepted. Classes are
intersected with the fixed evaluable gene list ``genes`` (configs/eval_genes.csv) when given, and
with genes having >= min_obs observed held-out lines in the fold:
- P1          mean per-gene Pearson r across held-out lines   (+ median, n genes r > 0.3,
              mean Spearman)
- P2a         per-line AUPRC across the class genes (label true < thr, score -pred), mean over
              lines; dP2a_vs_<B> = mean over lines of AP - AP_B
- P2b         per-gene normalised AP and AUROC across lines, genes with >= 3 pos and >= 3 neg;
              *_prev05 = restricted to held-out prevalence <= prev_max (the headline subset)
- vs baseline dP1_vs_<B> = mean over genes of r - r_B; frac_beat_<B> = fraction with r > r_B
The headline P2b is class 'selective', *_prev05.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import metrics as M

P1_KEYS = ["P1", "P1_median", "n_r_gt_0.3", "n_genes", "P1_spearman"]


def _mean(x):
    x = np.asarray(x, dtype=np.float64)
    x = x[np.isfinite(x)]
    return float(x.mean()) if len(x) else np.nan


def _classes_for(gene_classes, fold, genes):
    gc = gene_classes[fold] if isinstance(gene_classes, dict) else gene_classes
    return gc.reindex(genes).fillna(False).astype(bool)


def evaluate(pred_df: pd.DataFrame, true_df: pd.DataFrame, gene_classes, fold_of_line=None,
             baseline_preds=None, genes=None, min_obs=40, thr=-0.5, prev_max=0.5,
             spearman=True) -> dict:
    """Evaluate predictions against held-out truth.

    pred_df, true_df: [ModelID x gene]; pred must cover every (line, gene) of true_df/genes and be
      finite wherever truth is observed. true_df NaNs are excluded, never filled.
    gene_classes: boolean DataFrame (gene x class) or dict {fold: DataFrame} (per training fold).
    fold_of_line: Series ModelID -> fold, or None (single fold).
    baseline_preds: optional dict name -> DataFrame like pred_df (e.g. {'B0':..., 'B1':...}).
    genes: optional fixed evaluable gene list (configs/eval_genes.csv).
    min_obs: 40 for CV folds, 50 for the test set.

    Returns dict: summary {class: {metric: mean over folds, metric+'_sd', metric+'_folds'}},
      headline (flat dict), per_fold (DataFrame), per_gene (long DataFrame), per_line (DataFrame).
    """
    baseline_preds = baseline_preds or {}
    genes = pd.Index(true_df.columns if genes is None else [str(g) for g in genes])
    lines = true_df.index
    for nm, df in [("pred", pred_df)] + list(baseline_preds.items()):
        if not genes.isin(df.columns).all() or not lines.isin(df.index).all():
            raise ValueError(f"{nm} does not cover all evaluated lines/genes")
    if not genes.isin(true_df.columns).all():
        raise ValueError("genes not in true_df")
    T_all = true_df.loc[lines, genes]
    P_all = pred_df.loc[lines, genes]
    B_all = {b: df.loc[lines, genes] for b, df in baseline_preds.items()}
    fol = (pd.Series(0, index=lines) if fold_of_line is None
           else pd.Series(fold_of_line).reindex(lines))
    if fol.isna().any():
        raise ValueError("fold missing for some lines")
    classes = (list(gene_classes.values())[0] if isinstance(gene_classes, dict)
               else gene_classes).columns.tolist()

    fold_rows, gene_tabs, line_tabs = [], [], []
    for f in sorted(fol.unique()):
        ids = fol.index[fol == f]
        T = T_all.loc[ids].to_numpy(np.float64)
        P = P_all.loc[ids].to_numpy(np.float64)
        Bs = {b: d.loc[ids].to_numpy(np.float64) for b, d in B_all.items()}
        n_obs = (~np.isnan(T)).sum(0)
        r = M.per_gene_pearson(P, T, min_obs)
        rho = M.per_gene_spearman(P, T, min_obs) if spearman else np.full(len(genes), np.nan)
        p2b = M.per_gene_ap_norm(P, T, thr=thr, min_obs=min_obs)
        rB = {b: M.per_gene_pearson(B, T, min_obs) for b, B in Bs.items()}
        cls = _classes_for(gene_classes, f, genes)
        evalg = n_obs >= min_obs

        gt = pd.DataFrame({"fold": f, "gene": genes, "n_obs": n_obs, "r": r, "rho": rho,
                           "ap": p2b["ap"], "nap": p2b["nap"], "auroc": p2b["auroc"],
                           "prev": p2b["prev"], "n_pos": p2b["n_pos"]})
        for b in rB:
            gt[f"r_{b}"] = rB[b]
        for c in classes:
            gt[f"is_{c}"] = cls[c].to_numpy()
        gene_tabs.append(gt)

        lt = pd.DataFrame({"fold": f, "ModelID": ids,
                           "r_line": M.per_line_pearson(P[:, evalg], T[:, evalg])})
        for c in classes:
            m = cls[c].to_numpy() & evalg
            row = {"fold": f, "class": c, "n_lines": len(ids)}
            rc = r[m]
            row["P1"] = _mean(rc)
            row["P1_median"] = float(np.median(rc)) if len(rc) else np.nan
            row["n_r_gt_0.3"] = int((rc > 0.3).sum())
            row["n_genes"] = int(m.sum())
            row["P1_spearman"] = _mean(rho[m])
            # P2a: per-line AP across the class genes
            ap_l = (M.per_line_ap(P[:, m], T[:, m], thr) if m.any()
                    else np.full(len(ids), np.nan))
            row["P2a"] = _mean(ap_l)
            row["P2a_n_lines"] = int(np.isfinite(ap_l).sum())
            if c == "all":
                lt["ap"] = ap_l
            # P2b: per-gene normalised AP / AUROC
            e = m & np.isfinite(p2b["nap"])
            e05 = e & (p2b["prev"] <= prev_max)
            row.update({"P2b_nap": _mean(p2b["nap"][e]), "P2b_auroc": _mean(p2b["auroc"][e]),
                        "P2b_n": int(e.sum()), "P2b_nap_prev05": _mean(p2b["nap"][e05]),
                        "P2b_auroc_prev05": _mean(p2b["auroc"][e05]), "P2b_n_prev05": int(e05.sum())})
            for b, B in Bs.items():
                row[f"dP1_vs_{b}"] = _mean(rc - rB[b][m])
                row[f"frac_beat_{b}"] = float(np.mean(rc > rB[b][m])) if len(rc) else np.nan
                apb = (M.per_line_ap(B[:, m], T[:, m], thr) if m.any()
                       else np.full(len(ids), np.nan))
                row[f"dP2a_vs_{b}"] = _mean(ap_l - apb)
                if c == "all":
                    lt[f"ap_{b}"] = apb
            fold_rows.append(row)
        line_tabs.append(lt)

    per_fold = pd.DataFrame(fold_rows)
    metric_cols = [c for c in per_fold.columns if c not in ("fold", "class", "n_lines")]
    summary = {}
    for c in classes:
        sub = per_fold[per_fold["class"] == c]
        d = {}
        for k in metric_cols:
            v = sub[k].to_numpy(np.float64)
            d[k] = _mean(v)
            d[k + "_sd"] = float(np.nanstd(v, ddof=1)) if np.isfinite(v).sum() > 1 else np.nan
            d[k + "_folds"] = [None if not np.isfinite(x) else float(x) for x in v]
        summary[c] = d
    headline = {}
    for c in classes:
        for k in ["P1", "P1_median", "n_r_gt_0.3", "n_genes"] + [
                k for k in metric_cols if k.startswith(("dP1_vs_", "frac_beat_"))]:
            headline[f"{k}_{c}"] = summary[c][k]
            headline[f"{k}_{c}_sd"] = summary[c][k + "_sd"]
    for k in ["P2a"] + [k for k in metric_cols if k.startswith("dP2a_vs_")]:
        headline[k] = summary["all"][k]
        headline[k + "_sd"] = summary["all"][k + "_sd"]
    if "selective" in summary:
        for k in ["P2b_nap_prev05", "P2b_auroc_prev05", "P2b_n_prev05"]:
            headline[k] = summary["selective"][k]
            headline[k + "_sd"] = summary["selective"][k + "_sd"]
    return {"summary": summary, "headline": headline, "per_fold": per_fold,
            "per_gene": pd.concat(gene_tabs, ignore_index=True),
            "per_line": pd.concat(line_tabs, ignore_index=True)}
