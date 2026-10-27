"""G3 nonlinear-track tables (docs/g3_nonlinear_results.md) from the registry and stored per-gene results.

    python scripts/g3_report.py > cache/logs/g3/report.md

All metrics come from the frozen ``essml.evaluate`` outputs of each run (results/runs/<id>/per_fold.csv,
cache/runs/<id>/per_gene_fold.parquet). Nothing is recomputed from predictions except the pooled
(out-of-fold, all dev lines) numbers of the published-comparable block, which use the stored OOF
predictions and dev targets read through essml.guard (dev IDs only).

- Main table: P1 variable (mean ± fold SD), paired ΔP1 vs m1_E_resB0 (mean ± SE = SD/√n_folds,
  same folds), P1 all, P2b nAP/AUROC (selective, prevalence ≤ 0.5), fit + predict seconds.
- Subset runs (gene_subset in the config): P1 and ΔP1 recomputed per fold as the mean of the stored
  per-gene r over the subset genes that are dev-wide variable and evaluable in that fold, for the run
  and for m1_E_resB0 on exactly the same genes.
- Published-comparable block (loosely comparable; not used for selection): per-fold mean Spearman
  on the DeepDEP/LEAP sets, per-fold count of genes with r >= 0.2 (all eval genes), Shi-style count
  (r > 0.2 among genes with >= 5 held-out lines with effect < -0.5 and >= 5 with effect >= -0.5),
  each averaged over folds; and the pooled-OOF versions.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

REF = {  # D17 clean re-runs on the freeze-eval tree
    "b1": "20261026-220921_b1_0e2a7a27", "b2": "20261026-221456_b2_289f5170",
    "r1_E": "20261026-223037_r1_E_1df71977", "m1_E_resB0": "20261026-222117_m1_E_resB0_91e52f9d",
}
BASE = "m1_E_resB0"
SETS = ROOT / "configs" / "published_genesets.csv"


def registry():
    r = pd.read_csv(ROOT / "results" / "registry.csv", dtype=str)
    return r


def latest(reg, name, statuses=("ok",)):
    s = reg[(reg["name"] == name) & reg["status"].isin(statuses)]
    return None if s.empty else s.sort_values("timestamp").iloc[-1]


def per_fold(run_id):
    return pd.read_csv(ROOT / "results" / "runs" / run_id / "per_fold.csv")


def per_gene(run_id):
    return pd.read_parquet(ROOT / "cache" / "runs" / run_id / "per_gene_fold.parquet")


def summary(run_id):
    return json.loads((ROOT / "results" / "runs" / run_id / "summary.json").read_text())


def cfg_of(run_id):
    return summary(run_id)["config"]


def _subset_genes(cfg):
    p = cfg.get("params", {}).get("gene_subset")
    if not p:
        return None
    sub = pd.read_csv(ROOT / p, dtype={"entrez": str})["entrez"]
    n = cfg["params"].get("subset_first_n")
    return set(sub.iloc[: int(n)] if n else sub)


def fold_p1(run_id, genes=None, cls="variable"):
    """Per-fold P1 (Series fold -> value): from per_fold.csv, or over ``genes`` from per-gene r."""
    if genes is None:
        pf = per_fold(run_id)
        return pf[pf["class"] == cls].set_index("fold")["P1"]
    pg = per_gene(run_id)
    pg = pg[pg["gene"].astype(str).isin(genes) & pg[f"is_{cls}"]]
    return pg.groupby("fold")["r"].mean()


def row(name, run_id, base_id):
    cfg = cfg_of(run_id)
    sm = summary(run_id)
    genes = _subset_genes(cfg)
    p1 = fold_p1(run_id, genes)
    b = fold_p1(base_id, genes).reindex(p1.index)
    d = p1 - b
    n = len(d)
    pf = per_fold(run_id)
    pfa = pf[pf["class"] == "all"].set_index("fold")
    pfs = pf[pf["class"] == "selective"].set_index("fold")
    out = {"name": name, "run_id": run_id, "model": cfg["model"], "leaky": int(cfg.get("leaky", 0)),
           "folds": n, "genes": "subset500" if genes and len(genes) == 500 else
           (f"subset{len(genes)}" if genes else "all variable"),
           "P1_var": p1.mean(), "P1_var_sd": p1.std(ddof=1) if n > 1 else np.nan,
           "base_P1_var": b.mean(),
           "dP1": d.mean(), "dP1_se": d.std(ddof=1) / np.sqrt(n) if n > 1 else np.nan,
           "folds_better": int((d > 0).sum()),
           "P1_all": pfa["P1"].mean(), "P2b_nap": pfs["P2b_nap_prev05"].mean(),
           "P2b_auroc": pfs["P2b_auroc_prev05"].mean(),
           "fit_s": float(sm["fit_s"]), "predict_s": float(sm["predict_s"]), "threads": sm["threads"]}
    out["selected"] = bool(n == 5 and not out["leaky"] and out["dP1"] > 0.005
                           and out["dP1"] > out["dP1_se"])
    return out


def published_block(run_ids: dict, pooled=True):
    sets = pd.read_csv(SETS, dtype=str).dropna(subset=["entrez"])
    gs = {s: set(d["entrez"]) for s, d in sets.groupby("set")}
    rows = []
    Y = None
    if pooled:
        from essml import guard
        sp = pd.read_csv(ROOT / "configs" / "splits.csv", index_col=0)
        dev = sp.index[sp["role"] == "dev"]
        Y = guard.load_targets(dev).astype(np.float64)
        ev = pd.read_csv(ROOT / "configs" / "eval_genes.csv", dtype={"entrez": str})["entrez"]
        Y = Y[list(ev)]
    for name, rid in run_ids.items():
        pg = per_gene(rid)
        pg["gene"] = pg["gene"].astype(str)
        r = {"name": name, "run_id": rid}
        for s, g in gs.items():
            sub = pg[pg["gene"].isin(g)]
            r[f"rho_{s}"] = sub.groupby("fold")["rho"].mean().mean()
            r[f"n_{s}"] = int(sub.groupby("fold")["rho"].count().mean())
        r["n_r_ge_0.2"] = pg.assign(k=pg["r"] >= 0.2).groupby("fold")["k"].sum().mean()
        shi = pg[(pg["n_pos"] >= 5) & (pg["n_obs"] - pg["n_pos"] >= 5)]
        r["shi_eligible"] = shi.groupby("fold").size().mean()
        r["shi_n_r_gt_0.2"] = shi.assign(k=shi["r"] > 0.2).groupby("fold")["k"].sum().mean()
        if pooled:
            P = pd.read_parquet(ROOT / "cache" / "runs" / rid / "pred.parquet")
            idx = P.index.intersection(Y.index)
            T = Y.loc[idx].to_numpy()
            Pv = P.loc[idx, Y.columns].to_numpy(np.float64)
            from essml.metrics import per_gene_pearson, per_gene_spearman
            rp = per_gene_pearson(Pv, T, 40)
            rs = per_gene_spearman(Pv, T, 40)
            cols = pd.Index(Y.columns)
            for s, g in gs.items():
                m = cols.isin(g)
                r[f"pooled_rho_{s}"] = float(np.nanmean(rs[m]))
            r["pooled_n_r_ge_0.2"] = int(np.nansum(rp >= 0.2))
            obs = ~np.isnan(T)
            npos = ((T < -0.5) & obs).sum(0)
            nneg = ((T >= -0.5) & obs).sum(0)
            el = (npos >= 5) & (nneg >= 5)
            r["pooled_shi_eligible"] = int(el.sum())
            r["pooled_shi_n_r_gt_0.2"] = int(np.nansum(rp[el] > 0.2))
            r["pooled_n_lines"] = int(len(idx))
        rows.append(r)
    return pd.DataFrame(rows)


def fmt(x, d=4):
    return "—" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{d}f}"


def main(names=None, pooled=True):
    reg = registry()
    names = names or []
    base_id = REF[BASE]
    rows = [row(k, v, base_id) for k, v in REF.items()]
    found = {}
    for n in names:
        lr = latest(reg, n, ("ok", "partial_folds"))
        if lr is not None:
            found[n] = lr["run_id"]
            rows.append(row(n, lr["run_id"], base_id))
    T = pd.DataFrame(rows)
    print("| run | run id | genes | folds | P1 var ± SD | ΔP1 vs m1_E_resB0 ± SE (folds better) | "
          "P1 all | P2b nAP | P2b AUROC | fit+pred s | replaces |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    for _, r in T.iterrows():
        dag = "" if r["genes"] == "all variable" else " †"
        print(f"| {r['name']}{' (leaky)' if r['leaky'] else ''} | {r['run_id'][:15]} | {r['genes']} | "
              f"{r['folds']} | {fmt(r['P1_var'])} ± {fmt(r['P1_var_sd'])} | "
              f"{r['dP1']:+.4f} ± {fmt(r['dP1_se'])} ({r['folds_better']}/{r['folds']}) | "
              f"{fmt(r['P1_all'])}{dag} | {fmt(r['P2b_nap'])}{dag} | {fmt(r['P2b_auroc'])}{dag} | "
              f"{r['fit_s'] + r['predict_s']:.0f} | {'yes' if r['selected'] else 'no'} |")
    print("\n† subset run: P1 all and P2b cover every gene, most of them passed through (M3: M1; "
          "R2/R4: training mean); ΔP1 and 'replaces' refer to the subset genes only.\n")
    pb = published_block({**REF, **{n: found[n] for n in names
                                    if n in found and per_fold(found[n])['fold'].nunique() == 5}},
                         pooled=pooled)
    cols = [c for c in pb.columns if c not in ("run_id",)]
    print("| " + " | ".join(cols) + " |")
    print("|" + "---|" * len(cols))
    for _, r in pb.iterrows():
        print("| " + " | ".join(fmt(r[c], 4) if isinstance(r[c], float) else str(r[c]) for c in cols) + " |")
    out = {"main": T.to_dict(orient="records"), "published": pb.to_dict(orient="records")}
    (ROOT / "results" / "g3_nonlinear_summary.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    main(args, pooled="--no-pooled" not in sys.argv)
