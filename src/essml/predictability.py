"""Per-gene predictability tiers and the confidence score (DESIGN v3 §4 "Per-gene predictability", D22).

NOT frozen evaluation code: a diagnostic layer on top of the frozen ``essml.evaluate``. Methods and
thresholds are declared in ``configs/predictability.json`` (committed before any tier was computed).

* ``fold_r_table``: per gene and outer fold, the frozen evaluate's per-gene Pearson r (NaN where the
  gene has < min_obs observed held-out lines: evaluate's min_obs rule) plus P2b (nap, auroc, prev).
* ``fold_stats``: mean_r, min_r, sd_r (ddof 1) and n_folds_evaluable over the folds used.
* ``assign_tiers``: high (mean_r >= T_high and min_r >= 0.20), medium (not high, mean_r >= 0.15 and
  min_r > 0), low (other evaluable), not_evaluable (not evaluable in every fold used).
* ``confidence``: clip(mean_r - sd_r, -1, 1); NaN if not evaluable.
* ``null_calibrate``: the single pre-declared escalation of T_high (0.30 -> smallest of 0.35/0.40/0.50 at
  which the null control's high count is <= 5% of the observed high count).
* ``load_table`` / ``gene_annotation`` / ``annotate_long``: attach confidence and tier to predictions.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "configs" / "predictability.json"
TABLE = ROOT / "configs" / "predictable_genes.csv"

T_HIGH = 0.30          # DESIGN §4: high = mean_r >= 0.30 ...
HIGH_MIN_R = 0.20      # ... and min_r >= 0.20
MED_MEAN_R = 0.15      # medium = not high, mean_r >= 0.15 and min_r > 0
MED_MIN_R = 0.0        # strict
LADDER = (0.35, 0.40, 0.50)
MAX_FALSE_HIGH = 0.05
TIERS = ("high", "medium", "low", "not_evaluable")
EVALUABLE_TIERS = ("high", "medium", "low")


# ============================================================================ per-fold r (frozen evaluate)
def fold_r_table(pred: pd.DataFrame, truth: pd.DataFrame, folds: pd.Series, gene_classes: pd.DataFrame,
                 genes=None, min_obs=40) -> pd.DataFrame:
    """Long table (fold, gene, n_obs, r, nap, auroc, prev, n_pos) from ONE frozen ``evaluate`` call.

    pred/truth: [ModelID x gene] over the lines of ``folds``; genes default = truth's columns.
    """
    from .evaluate import evaluate
    genes = pd.Index(truth.columns if genes is None else [str(g) for g in genes])
    ids = pd.Index(folds.index)
    ev = evaluate(pred.loc[ids], truth.loc[ids], gene_classes, folds, genes=genes, min_obs=min_obs,
                  spearman=False)
    pg = ev["per_gene"]
    return pg[["fold", "gene", "n_obs", "r", "nap", "auroc", "prev", "n_pos"]].assign(gene=pg["gene"].astype(str))


def wide_r(long: pd.DataFrame, value="r") -> pd.DataFrame:
    """[gene x fold] of one per-fold value."""
    w = long.pivot(index="gene", columns="fold", values=value)
    w.columns = [int(c) for c in w.columns]
    return w


def fold_stats(r: pd.DataFrame, folds=None) -> pd.DataFrame:
    """mean_r, min_r, sd_r (ddof 1) and n_folds_evaluable over the given folds of an [gene x fold] r table.

    Statistics are computed over the evaluable folds; a gene is evaluable only if every fold is finite.
    """
    folds = list(r.columns) if folds is None else [int(f) for f in folds]
    R = r[folds].to_numpy(np.float64)
    fin = np.isfinite(R)
    n = fin.sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        Rm = np.where(fin, R, np.nan)
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            mean = np.nanmean(Rm, 1)
            mn = np.nanmin(Rm, 1)
            sd = np.nanstd(Rm, 1, ddof=1)
    out = pd.DataFrame({"mean_r": mean, "min_r": mn, "sd_r": sd, "n_folds_evaluable": n}, index=r.index)
    out["evaluable"] = n == len(folds)
    return out


def assign_tiers(stats: pd.DataFrame, t_high=T_HIGH, high_min_r=HIGH_MIN_R, med_mean_r=MED_MEAN_R,
                 med_min_r=MED_MIN_R) -> pd.Series:
    """Tier per gene (DESIGN §4). ``stats`` needs mean_r, min_r and evaluable (fold_stats output)."""
    ev = stats["evaluable"].to_numpy(bool)
    m = stats["mean_r"].to_numpy(np.float64)
    mn = stats["min_r"].to_numpy(np.float64)
    with np.errstate(invalid="ignore"):
        high = ev & (m >= t_high) & (mn >= high_min_r)
        med = ev & ~high & (m >= med_mean_r) & (mn > med_min_r)
    tier = np.where(~ev, "not_evaluable", np.where(high, "high", np.where(med, "medium", "low")))
    return pd.Series(tier, index=stats.index, name="tier")


def confidence(stats: pd.DataFrame) -> pd.Series:
    """clip(mean_r - sd_r, -1, 1) for evaluable genes, NaN otherwise."""
    c = np.clip(stats["mean_r"].to_numpy(np.float64) - stats["sd_r"].to_numpy(np.float64), -1.0, 1.0)
    c = np.where(stats["evaluable"].to_numpy(bool), c, np.nan)
    return pd.Series(c, index=stats.index, name="confidence")


def null_calibrate(obs_stats: pd.DataFrame, null_stats: pd.DataFrame, t0=T_HIGH, ladder=LADDER,
                   max_frac=MAX_FALSE_HIGH, **tier_kw) -> dict:
    """The single pre-declared escalation (DESIGN §4).

    Escalate only if n_high(null, t0) > max_frac * n_high(obs, t0); then T_high = the smallest ladder value
    with n_high(null, T) <= max_frac * n_high(obs, T). If none passes: the top of the ladder,
    status 'ladder_exhausted'. Returns {threshold, status, steps: [{t, n_high_obs, n_high_null, fdr, pass}]}.
    """
    steps = []
    for t in (t0, *ladder):
        no = int((assign_tiers(obs_stats, t_high=t, **tier_kw) == "high").sum())
        nn = int((assign_tiers(null_stats, t_high=t, **tier_kw) == "high").sum())
        steps.append({"t_high": float(t), "n_high_obs": no, "n_high_null": nn,
                      "est_fdr": (nn / no) if no > 0 else (0.0 if nn == 0 else float("inf")),
                      "pass": bool(nn <= max_frac * no)})
    if steps[0]["pass"]:
        return {"threshold": float(t0), "status": "not_triggered", "escalated": False, "steps": steps}
    for s in steps[1:]:
        if s["pass"]:
            return {"threshold": s["t_high"], "status": "escalated", "escalated": True, "steps": steps}
    return {"threshold": float(ladder[-1]), "status": "ladder_exhausted", "escalated": True, "steps": steps}


def build_table(r: pd.DataFrame, folds, t_high: float, symbols: pd.Series, source_run_id: str,
                all_genes=None) -> pd.DataFrame:
    """configs/predictable_genes.csv rows: one per gene of ``all_genes`` (default r's index)."""
    folds = [int(f) for f in folds]
    genes = pd.Index(r.index if all_genes is None else [str(g) for g in all_genes])
    rr = r.reindex(genes)
    st = fold_stats(rr, folds)
    out = pd.DataFrame(index=genes)
    out.index.name = "entrez"
    out["symbol"] = symbols.reindex(genes).to_numpy()
    for f in sorted(r.columns):
        out[f"r_fold{f}"] = rr[f].to_numpy()
    out["mean_r"] = np.where(st["evaluable"], st["mean_r"], np.nan)
    out["min_r"] = np.where(st["evaluable"], st["min_r"], np.nan)
    out["sd_r"] = np.where(st["evaluable"], st["sd_r"], np.nan)
    out["n_folds_evaluable"] = st["n_folds_evaluable"].astype(int)
    out["confidence"] = confidence(st)
    out["tier"] = assign_tiers(st, t_high=t_high)
    out["source_run_id"] = source_run_id
    out["folds_used"] = ";".join(map(str, folds))
    out["threshold_high_mean_r"] = float(t_high)
    return out.reset_index()


# ============================================================================ annotation of predictions
def sha256_file(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_table(path=TABLE) -> pd.DataFrame:
    """predictable_genes table indexed by Entrez string."""
    t = pd.read_csv(path, dtype={"entrez": str, "symbol": str, "tier": str, "source_run_id": str})
    if not {"entrez", "confidence", "tier"} <= set(t.columns):
        raise ValueError(f"{path}: not a predictable-genes table")
    bad = set(t["tier"]) - set(TIERS)
    if bad:
        raise ValueError(f"{path}: unknown tiers {bad}")
    return t.set_index("entrez")


def gene_annotation(genes, table: pd.DataFrame) -> pd.DataFrame:
    """[gene] -> symbol, confidence, tier, mean_r (genes absent from the table: not_evaluable, NaN)."""
    genes = pd.Index([str(g) for g in genes])
    t = table.reindex(genes)
    out = pd.DataFrame({"entrez": genes,
                        "symbol": t["symbol"].to_numpy() if "symbol" in t else None,
                        "confidence": t["confidence"].astype(np.float64).to_numpy(),
                        "tier": t["tier"].fillna("not_evaluable").to_numpy(),
                        "mean_r_dev_cv": t["mean_r"].astype(np.float64).to_numpy() if "mean_r" in t else np.nan})
    return out


def annotate_long(pred: pd.DataFrame, table: pd.DataFrame, value_name="predicted_effect") -> pd.DataFrame:
    """Wide [line x gene] predictions -> tidy rows (ModelID, entrez, symbol, predicted_effect, confidence, tier)."""
    P = pred.copy()
    P.columns = P.columns.astype(str)
    ann = gene_annotation(P.columns, table).set_index("entrez")
    idx_name = P.index.name or "ModelID"
    long = P.rename_axis(index=idx_name, columns="entrez").stack(future_stack=True).rename(value_name).reset_index()
    long = long.rename(columns={idx_name: "ModelID"})
    long = long.join(ann[["symbol", "confidence", "tier"]], on="entrez")
    return long[["ModelID", "entrez", "symbol", value_name, "confidence", "tier"]]


def table_provenance(path=TABLE) -> dict:
    t = load_table(path)
    return {"file": str(path), "sha256": sha256_file(path),
            "source_run_id": sorted(set(t["source_run_id"].dropna()))[:3],
            "threshold_high_mean_r": sorted(set(t["threshold_high_mean_r"].dropna().astype(float))),
            "tier_counts": t["tier"].value_counts().to_dict()}


def table_for_run(run_id: str, config=CONFIG):
    """Tier table of the model behind ``run_id`` (a CV run id '<stamp>_<name>_<hash8>' or a predict.py tag
    'predict_<name>_<hash12>'): the primary table for the G4 headline, the secondary table for a secondary source
    (configs/predictability.json); None if the run's model has no tier table of its own."""
    import re
    cfg = load_config(config)
    for src in [cfg["source"]] + list(cfg.get("secondary_sources", [])):
        if re.search(rf"(^|_){re.escape(src['run_name'])}_[0-9a-f]{{8,12}}$", str(run_id)):
            return ROOT / src["table"].split(" ")[0]
    return None


def load_config(path=CONFIG) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))
