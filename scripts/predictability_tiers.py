"""Per-gene predictability tiers (DESIGN §4 "Per-gene predictability", D22; configs/predictability.json).

    python scripts/predictability_tiers.py --threads 4                         # frozen table (headline run)
    python scripts/predictability_tiers.py --exclude-folds 0 --dry-run-tiers   # G4 dry-run tiers (folds 1-4)
    python scripts/predictability_tiers.py --run <run_id> --null-run <run_id>  # e.g. if M5 becomes the headline

Inputs: dev-CV out-of-fold predictions cache/runs/<run_id>/pred.parquet of the headline run and of the
null control (re-run through scripts/run.py if absent), dev targets through essml.guard (dev IDs only).
Per gene and fold, r comes from the frozen essml.evaluate (min_obs rule); tiers, confidence and the single
null-calibration escalation from essml.predictability.

Outputs (frozen run): configs/predictable_genes.csv, results/predictability/distribution.json,
docs/figures/predictability_*.png. Dry-run tiers: results/dryrun/predictability/.
"""
from __future__ import annotations

import os
import sys

if __name__ == "__main__":
    _n = sys.argv[sys.argv.index("--threads") + 1] if "--threads" in sys.argv else "4"
    for _v in ["OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"]:
        os.environ[_v] = _n

import argparse  # noqa: E402
import datetime as dt  # noqa: E402
import json  # noqa: E402
import subprocess  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from essml import genesets, guard  # noqa: E402
from essml import predictability as PR  # noqa: E402

QUANTILES = [0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99]
BINS = np.round(np.arange(-0.30, 1.0001, 0.05), 2)


def clean(o):
    """Strict JSON: non-finite floats -> None, numpy scalars -> Python."""
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, (float, np.floating)):
        return float(o) if np.isfinite(o) else None
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return o


def git_head():
    return subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()


def dev_frame():
    sp = pd.read_csv(ROOT / "configs" / "splits.csv", index_col=0)
    sp.index = sp.index.astype(str)
    return sp[sp["role"] == "dev"]


def oof_predictions(run_id: str, config: str, threads: int, log=print) -> tuple[pd.DataFrame, str]:
    p = ROOT / "cache" / "runs" / run_id / "pred.parquet"
    if not p.exists():  # DESIGN §4 source: re-run the config through the launcher
        log(f"{p} missing: re-running {config} through scripts/run.py")
        sys.path.insert(0, str(ROOT / "scripts"))
        import run as RUN
        run_id = RUN.main(["--config", config, "--threads", str(threads)])
        p = ROOT / "cache" / "runs" / run_id / "pred.parquet"
    P = pd.read_parquet(p).astype(np.float64)
    P.columns = P.columns.astype(str)
    P.index = P.index.astype(str)
    return P, run_id


def per_fold_r(run_id, P, Y, folds, gc, min_obs, log=print):
    long = PR.fold_r_table(P, Y, folds, gc, genes=P.columns, min_obs=min_obs)
    cached = ROOT / "cache" / "runs" / run_id / "per_gene_fold.parquet"
    check = None
    if cached.exists():  # the frozen run's own evaluation on the eval genes must be reproduced
        c = pd.read_parquet(cached, columns=["fold", "gene", "r"])
        c["gene"] = c["gene"].astype(str)
        m = c.merge(long[["fold", "gene", "r"]], on=["fold", "gene"], suffixes=("_run", "_here"))
        a, b = m["r_run"].to_numpy(), m["r_here"].to_numpy()
        same_nan = bool((np.isnan(a) == np.isnan(b)).all())
        d = float(np.nanmax(np.abs(a - b))) if len(m) else np.nan
        check = {"n_compared": int(len(m)), "max_abs_diff": d, "same_nan_pattern": same_nan,
                 "ok": bool(same_nan and len(m) == len(c) and d <= 1e-6)}  # pred.parquet is float32
        if not check["ok"]:
            raise AssertionError(f"{run_id}: recomputed per-fold r differs from the run's cached evaluation: {check}")
        log(f"  {run_id}: per-fold r reproduces the run's cached evaluation ({len(m)} gene-folds, max diff {d:.1e})")
    return long, check


def qdict(x):
    x = np.asarray(x, np.float64)
    x = x[np.isfinite(x)]
    return {str(q): float(np.quantile(x, q)) for q in QUANTILES} if len(x) else {}


def hist(x):
    x = np.asarray(x, np.float64)
    x = x[np.isfinite(x)]
    c, _ = np.histogram(np.clip(x, BINS[0], BINS[-1]), bins=BINS)
    return {"bin_edges": BINS.tolist(), "counts": c.astype(int).tolist(),
            "note": "values outside [-0.30, 1.00] are counted in the edge bins"}


def describe(x):
    x = np.asarray(x, np.float64)
    x = x[np.isfinite(x)]
    if not len(x):
        return {"n": 0}
    return {"n": int(len(x)), "mean": float(x.mean()), "median": float(np.median(x)),
            "q25": float(np.quantile(x, 0.25)), "q75": float(np.quantile(x, 0.75))}


def p2b_by_tier(long: pd.DataFrame, table: pd.DataFrame, gc: pd.DataFrame, folds_used) -> dict:
    """Per-gene P2b (mean over folds where defined) by tier, selective genes; also restricted to genes
    whose fold prevalence is <= 0.5 in every fold where P2b is defined (headline-style)."""
    sub = long[long["fold"].isin(folds_used)]
    nap = sub.pivot(index="gene", columns="fold", values="nap")
    auc = sub.pivot(index="gene", columns="fold", values="auroc")
    prev = sub.pivot(index="gene", columns="fold", values="prev")
    defined = nap.notna()
    prev_ok = (prev.where(defined) <= 0.5) | ~defined
    g = pd.DataFrame({"nap": nap.mean(axis=1), "auroc": auc.mean(axis=1), "n_folds_defined": defined.sum(axis=1),
                      "prev05_all_folds": prev_ok.all(axis=1) & defined.any(axis=1)})
    t = table.set_index("entrez")["tier"]
    sel = gc["selective"].reindex(g.index).fillna(False).astype(bool)
    out = {"definition": "selective genes (configs/gene_classes.csv); per gene the mean over folds of the frozen "
                         "evaluate's nap / AUROC where defined (>= 3 positives and >= 3 negatives in the fold); "
                         "'prev05' = fold prevalence <= 0.5 in every fold where defined (headline subset)"}
    for name, mask in (("selective", sel), ("selective_prev05", sel & g["prev05_all_folds"])):
        d = {}
        for tier in PR.TIERS:
            m = mask & (t.reindex(g.index) == tier)
            d[tier] = {"nap": describe(g.loc[m, "nap"]), "auroc": describe(g.loc[m, "auroc"])}
        out[name] = d
    return out


def make_figures(table: pd.DataFrame, null_table: pd.DataFrame, gc: pd.DataFrame, t_high: float, fig_dir: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    fig_dir.mkdir(parents=True, exist_ok=True)
    ink, muted, grid = "#1f2328", "#59636e", "#d8dee4"
    tier_col = {"high": "#2a78d6", "medium": "#eb6834", "low": "#1baf7a"}  # dataviz reference slots 1-3 (all-pairs safe)
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": muted, "axes.labelcolor": ink, "xtick.color": muted,
                         "ytick.color": muted, "axes.spines.top": False, "axes.spines.right": False})
    ev = table[table["tier"] != "not_evaluable"].set_index("entrez")
    var = gc["variable"].reindex(ev.index).fillna(False).astype(bool)
    nev = null_table[null_table["tier"] != "not_evaluable"].set_index("entrez")
    nvar = gc["variable"].reindex(nev.index).fillna(False).astype(bool)

    # 1. distribution of mean_r (all evaluable | variable genes), each with the permuted-expression control
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.2), dpi=150, sharex=True)
    bins = np.arange(-0.3, 1.0001, 0.02)
    for ax, (title, sel, nsel) in zip(axes, (("all evaluable genes", slice(None), slice(None)),
                                             ("variable genes (SD > 0.25)", var, nvar))):
        x = ev.loc[sel, "mean_r"].clip(-0.3, 1)
        ax.hist(x, bins=bins, color="#2a78d6", label=f"m1_E_resB0 (n={len(x)})")
        ax.hist(nev.loc[nsel, "mean_r"].clip(-0.3, 1), bins=bins, histtype="step", color=ink, lw=1.0,
                label="permuted-expression control")
        top = ax.get_ylim()[1]
        for xv, lab, yf in ((0.15, "0.15", 0.97), (t_high, f"{t_high:.2f}", 0.97)):
            ax.axvline(xv, color=muted, lw=0.8, ls="--")
            ax.text(xv + 0.01, top * yf, lab, color=muted, va="top", fontsize=8)
        ax.set_title(title, fontsize=9, color=ink, loc="left")
        ax.set_xlabel("mean per-gene r over 5 dev-CV folds")
        ax.grid(axis="y", color=grid, lw=0.5)
        ax.legend(frameon=False, fontsize=7.5, loc="upper right", bbox_to_anchor=(1.0, 0.82))
    axes[0].set_ylabel("genes")
    fig.tight_layout()
    fig.savefig(fig_dir / "predictability_mean_r_hist.png")
    plt.close(fig)

    # 2. mean_r vs min_r coloured by tier
    fig, ax = plt.subplots(figsize=(4.6, 4.0), dpi=150)
    for tier in ("low", "medium", "high"):
        s = ev[ev["tier"] == tier]
        ax.scatter(s["mean_r"], s["min_r"], s=2, alpha=0.35, color=tier_col[tier], lw=0, rasterized=True)
    ax.plot([t_high, 1], [0.2, 0.2], color=ink, lw=0.8)
    ax.plot([t_high, t_high], [0.2, 1], color=ink, lw=0.8)
    ax.plot([0.15, 1], [0, 0], color=muted, lw=0.8, ls="--")
    ax.plot([0.15, 0.15], [0, 1], color=muted, lw=0.8, ls="--")
    ax.set_xlim(-0.3, 1)
    ax.set_ylim(-0.6, 1)
    ax.set_xlabel("mean_r (5 folds)")
    ax.set_ylabel("min_r (worst fold)")
    handles = [Line2D([], [], marker="o", ls="", color=tier_col[t], label=f"{t} (n={int((ev['tier'] == t).sum())})")
               for t in ("high", "medium", "low")]
    ax.legend(handles=handles, frameon=False, fontsize=8, loc="upper left")
    ax.grid(color=grid, lw=0.5)
    fig.tight_layout()
    fig.savefig(fig_dir / "predictability_tiers_scatter.png")
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", default=str(ROOT / "configs" / "predictability.json"))
    ap.add_argument("--run", default=None, help="headline CV run id (default: config source.run_id)")
    ap.add_argument("--run-config", default=None, help="config to re-run if the run's predictions are missing")
    ap.add_argument("--null-run", default=None)
    ap.add_argument("--exclude-folds", default="", help="comma-separated outer folds left out (dry run: 0)")
    ap.add_argument("--dry-run-tiers", action="store_true", help="write to results/dryrun/predictability (G4 dry run)")
    ap.add_argument("--no-figures", action="store_true")
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args(argv)
    cfg = PR.load_config(a.config)
    run_id = a.run or cfg["source"]["run_id"]
    null_id = a.null_run or cfg["null_control"]["run_id"]
    min_obs = int(cfg["per_fold_r"]["min_obs"])
    excl = [int(x) for x in a.exclude_folds.split(",") if x.strip()]
    folds_used = [f for f in cfg["per_fold_r"]["folds"] if f not in excl]
    dry = a.dry_run_tiers or bool(excl)
    if dry:
        out_dir = ROOT / cfg["dry_run"]["out_dir"]
        table_path = out_dir / f"predictable_genes_excl_fold{'_'.join(map(str, excl)) or 'none'}.csv"
        dist_path = out_dir / f"distribution_excl_fold{'_'.join(map(str, excl)) or 'none'}.json"
    else:
        out_dir = ROOT / "results" / "predictability"
        table_path = ROOT / "configs" / "predictable_genes.csv"
        dist_path = out_dir / "distribution.json"
    out_dir.mkdir(parents=True, exist_ok=True)

    dev = dev_frame()
    folds = dev["cv_fold"].astype(int)
    print(f"dev lines {len(dev)}; folds used {folds_used}; headline run {run_id}; null run {null_id}")
    Y = guard.load_targets(list(dev.index)).astype(np.float64)
    Y.columns = Y.columns.astype(str)
    gc = genesets.load_gene_classes()
    gc.index = gc.index.astype(str)
    sym = pd.read_csv(ROOT / "configs" / "gene_classes.csv", dtype={"entrez": str}).set_index("entrez")["symbol"]

    P, run_id = oof_predictions(run_id, a.run_config or cfg["source"]["config"], a.threads)
    long, chk = per_fold_r(run_id, P, Y, folds, gc, min_obs)
    del P
    Pn, null_id = oof_predictions(null_id, cfg["null_control"]["config"], a.threads)
    long_n, chk_n = per_fold_r(null_id, Pn, Y, folds, gc, min_obs)
    del Pn, Y

    r, rn = PR.wide_r(long), PR.wide_r(long_n)
    st, stn = PR.fold_stats(r, folds_used), PR.fold_stats(rn, folds_used)
    cal = PR.null_calibrate(st, stn, ladder=tuple(cfg["null_calibration"]["ladder"]),
                            max_frac=float(cfg["null_calibration"]["max_false_high_fraction"]))
    t_high = cal["threshold"]
    print(f"null calibration: {cal['status']}, T_high = {t_high}; steps {cal['steps']}")
    genes_all = list(r.index)
    table = PR.build_table(r, folds_used, t_high, sym, run_id, all_genes=genes_all)
    null_table = PR.build_table(rn, folds_used, t_high, sym, null_id, all_genes=genes_all)
    keep_cols = ["entrez", "symbol"] + [f"r_fold{f}" for f in sorted(r.columns)] + [
        "mean_r", "min_r", "sd_r", "n_folds_evaluable", "confidence", "tier", "source_run_id", "folds_used",
        "threshold_high_mean_r"]
    out = table[keep_cols].copy()
    num = [c for c in out.columns if c.startswith("r_fold")] + ["mean_r", "min_r", "sd_r", "confidence"]
    out[num] = out[num].round(6)
    out.to_csv(table_path, index=False, lineterminator="\n")
    print(f"wrote {table_path.relative_to(ROOT)}: {out['tier'].value_counts().to_dict()}")

    # ---- distribution summary
    t = table.set_index("entrez")
    tn = null_table.set_index("entrez")
    ev = t["tier"] != "not_evaluable"
    var = gc["variable"].reindex(t.index).fillna(False).astype(bool)
    classes = ["variable", "selective", "rare_selective", "common_essential", "nonessential"]
    tier_counts = {"all": t["tier"].value_counts().reindex(PR.TIERS, fill_value=0).astype(int).to_dict()}
    for c in classes:
        m = gc[c].reindex(t.index).fillna(False).astype(bool)
        tier_counts[c] = t.loc[m, "tier"].value_counts().reindex(PR.TIERS, fill_value=0).astype(int).to_dict()
    null_counts = {"all": tn["tier"].value_counts().reindex(PR.TIERS, fill_value=0).astype(int).to_dict(),
                   "variable": tn.loc[var, "tier"].value_counts().reindex(PR.TIERS, fill_value=0).astype(int).to_dict()}
    nh, oh = null_counts["all"]["high"], tier_counts["all"]["high"]
    nm, om = null_counts["all"]["medium"], tier_counts["all"]["medium"]
    dist = {
        "generated": dt.datetime.now().isoformat(timespec="seconds"), "git_commit": git_head(),
        "config": "configs/predictability.json", "source_run_id": run_id, "null_run_id": null_id,
        "folds_used": folds_used, "excluded_folds": excl, "min_obs": min_obs,
        "dry_run_note": (cfg["dry_run"]["purpose"] if dry else None),
        "table": str(table_path.relative_to(ROOT).as_posix()), "table_sha256": PR.sha256_file(table_path),
        "reproduces_cached_run_evaluation": {"headline": chk, "null": chk_n},
        "null_calibration": cal, "threshold_high_mean_r": t_high,
        "n_genes": int(len(t)), "n_evaluable": int(ev.sum()), "n_variable_evaluable": int((ev & var).sum()),
        "tier_counts": tier_counts, "null_tier_counts": null_counts,
        "estimated_fdr_high": (nh / oh) if oh else None, "null_over_observed_medium": (nm / om) if om else None,
        "mean_r": {"all_evaluable": {"quantiles": qdict(t.loc[ev, "mean_r"]), "mean": float(t.loc[ev, "mean_r"].mean()),
                                     "histogram": hist(t.loc[ev, "mean_r"])},
                   "variable": {"quantiles": qdict(t.loc[ev & var, "mean_r"]),
                                "mean": float(t.loc[ev & var, "mean_r"].mean()),
                                "histogram": hist(t.loc[ev & var, "mean_r"])},
                   "null_variable": {"quantiles": qdict(tn.loc[ev & var, "mean_r"])}},
        "fraction_mean_r_ge_0.3": {"all_evaluable": float((t.loc[ev, "mean_r"] >= 0.3).mean()),
                                   "variable": float((t.loc[ev & var, "mean_r"] >= 0.3).mean())},
        "fraction_mean_r_ge_0.2": {"all_evaluable": float((t.loc[ev, "mean_r"] >= 0.2).mean()),
                                   "variable": float((t.loc[ev & var, "mean_r"] >= 0.2).mean())},
        "per_tier_mean_r": {tier: describe(t.loc[t["tier"] == tier, "mean_r"]) for tier in PR.EVALUABLE_TIERS},
        "per_tier_confidence": {tier: describe(t.loc[t["tier"] == tier, "confidence"]) for tier in PR.EVALUABLE_TIERS},
        "p2b_by_tier": p2b_by_tier(long, table, gc, folds_used),
    }
    dist_path.write_text(json.dumps(clean(dist), indent=1))
    print(f"wrote {dist_path.relative_to(ROOT)}; est. FDR(high) = {dist['estimated_fdr_high']}")
    if not dry and not a.no_figures:
        make_figures(table, null_table, gc, t_high, ROOT / "docs" / "figures")
        print("wrote docs/figures/predictability_*.png")
    return dist


if __name__ == "__main__":
    main()
