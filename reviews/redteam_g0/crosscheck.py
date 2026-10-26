"""Independent red-team cross-check of the committed B1 run (G0/G1).

Does NOT import essml.metrics / essml.evaluate / essml.baselines / scripts/run.py. Uses only
essml.guard.load_targets for DEV ModelIDs (no non-dev access) and sklearn's StratifiedGroupKFold
(only to reproduce the inner split indices; the scoring is my own code).

    python reviews/redteam_g0/crosscheck.py

Computes, per outer dev fold, from (a) the saved OOF predictions in cache/runs/<b1>/pred.parquet
and (b) an independent re-implementation of B1 (shrunk lineage mean, kappa re-tuned by my own
inner-CV code):
  - P1 = mean per-gene Pearson r over held-out lines, variable genes (configs/gene_classes.csv)
    restricted to configs/eval_genes.csv and >= 40 observed held-out lines; constant -> 0;
    averaged over folds.
  - P2a = per-line average precision (label truth < -0.5, score -pred, ties grouped as in a
    step-wise PR curve) across the class-'all' evaluable genes, mean over lines, mean over folds.
  - P2b_nap_prev05 (selective, >= 3 pos and >= 3 neg, held-out prevalence <= 0.5), extra.
and compares with results/registry.csv.
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from essml import guard  # noqa: E402  (dev-only target loader)

RUN = "20261026-151647_b1_0e2a7a27"
THR = -0.5
MIN_OBS = 40
GRID = [0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 10000]


# ------------------------------------------------------------------ own metric code
def pearson_one(p, t):
    """Pearson r of two 1-D float64 vectors (already NaN-free). Constant -> 0."""
    if p.max() == p.min() or t.max() == t.min():
        return 0.0
    pc, tc = p - p.mean(), t - t.mean()
    d = np.sqrt((pc @ pc) * (tc @ tc))
    return 0.0 if d == 0 else float(np.clip((pc @ tc) / d, -1, 1))


def per_gene_r(P, T, min_obs):
    out = np.full(P.shape[1], np.nan)
    for g in range(P.shape[1]):
        o = ~np.isnan(T[:, g])
        if o.sum() >= max(min_obs, 2):
            out[g] = pearson_one(P[o, g], T[o, g])
    return out


def ap_stepwise(y, s):
    """Average precision, sum over distinct thresholds of (R_k - R_{k-1}) * P_k (ties grouped)."""
    y = np.asarray(y, bool)
    order = np.argsort(-s, kind="mergesort")
    s, y = s[order], y[order]
    # last index of each block of tied scores
    ends = np.r_[np.flatnonzero(np.diff(s) != 0), len(s) - 1]
    tp = np.cumsum(y)[ends].astype(float)
    fp = (ends + 1) - tp
    prec = tp / (tp + fp)
    rec = tp / y.sum()
    rec_prev = np.r_[0.0, rec[:-1]]
    return float(np.sum((rec - rec_prev) * prec))


def per_line_ap(P, T):
    out = np.full(P.shape[0], np.nan)
    for i in range(P.shape[0]):
        o = ~np.isnan(T[i])
        y = T[i, o] < THR
        if 0 < y.sum() < y.size:
            out[i] = ap_stepwise(y, -P[i, o])
    return out


def per_gene_nap(P, T):
    nap = np.full(P.shape[1], np.nan)
    prev = np.full(P.shape[1], np.nan)
    for g in range(P.shape[1]):
        o = ~np.isnan(T[:, g])
        y = T[o, g] < THR
        if o.sum():
            prev[g] = y.mean()
        if y.sum() >= 3 and (~y).sum() >= 3 and o.sum() >= MIN_OBS:
            ap = ap_stepwise(y, -P[o, g])
            nap[g] = (ap - prev[g]) / (1 - prev[g])
    return nap, prev


# ------------------------------------------------------------------ own B1
def b1_fit_predict(Ytr, lin_tr, lin_te, kappa):
    Yv = Ytr.to_numpy(np.float64)
    obs = ~np.isnan(Yv)
    Y0 = np.where(obs, Yv, 0.0)
    mu = Y0.sum(0) / obs.sum(0)
    out = np.empty((len(lin_te), Yv.shape[1]))
    groups = {}
    for lv in np.unique(lin_tr):
        m = lin_tr == lv
        groups[lv] = (Y0[m].sum(0), obs[m].sum(0))
    for r, lv in enumerate(lin_te):
        if lv in groups:
            S, N = groups[lv]
            den = N + kappa
            with np.errstate(invalid="ignore", divide="ignore"):
                v = (S + kappa * mu) / den
            out[r] = np.where(den > 0, v, mu)
        else:
            out[r] = mu
    return out


def own_variable(Y):
    """variable = SD (ddof=1, observed lines) > 0.25."""
    return Y.std(axis=0, ddof=1, skipna=True) > 0.25


def tune_kappa(Ytr, lin, grp, eval_genes, seed):
    from sklearn.model_selection import StratifiedGroupKFold
    ids = Ytr.index
    sgk = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=seed)
    scores = np.zeros((len(GRID), 5))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        splits = list(sgk.split(np.zeros(len(ids)), lin.loc[ids].to_numpy(), grp.loc[ids].to_numpy()))
    for j, (a, b) in enumerate(splits):
        itr, iva = ids[a], ids[b]
        Yi = Ytr.loc[itr]
        var = own_variable(Yi)
        g = var.index[var.to_numpy()].intersection(eval_genes)
        T = Ytr.loc[iva, g].to_numpy(np.float64)
        for i, k in enumerate(GRID):
            P = b1_fit_predict(Yi[g], lin.loc[itr].to_numpy(), lin.loc[iva].to_numpy(), k)
            r = per_gene_r(P, T, 30)
            scores[i, j] = np.nanmean(r)
    return GRID[int(np.argmax(scores.mean(1)))], scores.mean(1)


# ------------------------------------------------------------------ main
def evaluate_own(pred, Y, folds, gcl, eval_genes):
    var_g = gcl.index[gcl["variable"]].intersection(eval_genes)
    sel_g = gcl.index[gcl["selective"]].intersection(eval_genes)
    rows = []
    for f in sorted(folds.unique()):
        ids = folds.index[folds == f]
        T = Y.loc[ids, eval_genes].to_numpy(np.float64)
        P = pred.loc[ids, eval_genes].to_numpy(np.float64)
        nobs = (~np.isnan(T)).sum(0)
        ev = nobs >= MIN_OBS
        isvar = eval_genes.isin(var_g) & ev
        issel = eval_genes.isin(sel_g) & ev
        r = per_gene_r(P[:, isvar], T[:, isvar], MIN_OBS)
        ap_all = per_line_ap(P[:, ev], T[:, ev])
        nap, prev = per_gene_nap(P[:, issel], T[:, issel])
        e05 = np.isfinite(nap) & (prev <= 0.5)
        rows.append({"fold": f, "n_var_genes": int(isvar.sum()), "P1_variable": np.nanmean(r),
                     "P2a": np.nanmean(ap_all), "P2b_nap_prev05": float(np.mean(nap[e05])),
                     "n_P2b_prev05": int(e05.sum())})
    return pd.DataFrame(rows)


def main():
    sp = pd.read_csv(ROOT / "configs" / "splits.csv", index_col=0)
    dev = sp[sp["role"] == "dev"]
    assert (guard.roles_of(dev.index) == "dev").all()
    folds = dev["cv_fold"].astype(int)
    lin = dev["OncotreeLineage"].astype(str)
    grp = dev["PatientID"].astype(str)
    Y = guard.load_targets(dev.index).astype(np.float64)  # DEV ONLY
    eval_genes = pd.Index(pd.read_csv(ROOT / "configs" / "eval_genes.csv", dtype={"entrez": str})["entrez"])
    gcl = pd.read_csv(ROOT / "configs" / "gene_classes.csv", dtype={"entrez": str}).set_index("entrez")
    for c in ["variable", "selective"]:
        gcl[c] = gcl[c].astype(str).str.lower().eq("true")

    out = {}
    # sanity: committed variable class vs my own dev-wide recomputation
    myvar = own_variable(Y)
    out["variable_class_mismatches"] = int((myvar.reindex(gcl.index).fillna(False) != gcl["variable"]).sum())
    # sanity: eval_genes satisfy the >= 40 per dev fold rule (dev half of the rule only)
    cnt = pd.concat({f: Y.loc[folds.index[folds == f]].notna().sum() for f in sorted(folds.unique())}, axis=1)
    out["eval_genes_failing_dev_rule"] = int((cnt.loc[eval_genes].min(axis=1) < 40).sum())
    out["genes_passing_dev_rule_not_in_eval"] = int(((cnt.min(axis=1) >= 40) & ~cnt.index.isin(eval_genes)).sum())

    # (a) saved OOF predictions of the committed run
    pred_saved = pd.read_parquet(ROOT / "cache" / "runs" / RUN / "pred.parquet").astype(np.float64)
    pred_saved.columns = pred_saved.columns.astype(str)
    a = evaluate_own(pred_saved, Y, folds, gcl, eval_genes)

    # (b) independent B1 with my own inner-CV kappa tuning (inner seed = split_seed + fold)
    seed = json.loads((ROOT / "configs" / "seeds.json").read_text())["split_seed"]
    preds, chosen = [], {}
    for f in sorted(folds.unique()):
        tr, ho = folds.index[folds != f], folds.index[folds == f]
        k, _ = tune_kappa(Y.loc[tr], lin, grp, eval_genes, seed + int(f))
        chosen[int(f)] = k
        P = b1_fit_predict(Y.loc[tr], lin.loc[tr].to_numpy(), lin.loc[ho].to_numpy(), k)
        preds.append(pd.DataFrame(P, index=ho, columns=Y.columns))
    pred_own = pd.concat(preds)
    b = evaluate_own(pred_own, Y, folds, gcl, eval_genes)
    out["own_chosen_kappa"] = chosen
    out["max_abs_pred_diff_saved_vs_own"] = float(np.nanmax(np.abs(
        pred_saved.loc[pred_own.index, eval_genes].to_numpy() - pred_own.loc[:, eval_genes].to_numpy())))

    reg = pd.read_csv(ROOT / "results" / "registry.csv")
    rr = reg[reg["run_id"] == RUN].iloc[0]
    out["registry"] = {k: float(rr[k]) for k in ["P1_variable", "P1_variable_sd", "P2a", "P2a_sd",
                                                 "P2b_nap_prev05"]}
    for nm, t in [("saved_pred", a), ("own_B1", b)]:
        out[nm] = {"per_fold": t.round(6).to_dict("records"),
                   "P1_variable": float(t["P1_variable"].mean()),
                   "P1_variable_sd": float(t["P1_variable"].std(ddof=1)),
                   "P2a": float(t["P2a"].mean()), "P2a_sd": float(t["P2a"].std(ddof=1)),
                   "P2b_nap_prev05": float(t["P2b_nap_prev05"].mean())}
    print(json.dumps(out, indent=1, default=float))
    (Path(__file__).parent / "crosscheck_output.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
