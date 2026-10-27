"""Mutation-driven dependency diagnostic (DESIGN.md §4, RESEARCH_TREE D24). Dev CV only.

    python scripts/mutation_dependencies.py pairs   # -> configs/oncogene_pairs.csv (features only)
    python scripts/mutation_dependencies.py run     # -> results/mutation_dependencies.json,
                                                    #    docs/figures/mutation_dependencies.png

``pairs`` counts biomarker-positive dev lines from the raw feature matrices only (no CRISPR) and
must be committed, with configs/mutation_dependencies.json, before ``run`` computes any metric.
``run`` reads the out-of-fold predictions (cache/runs/<run_id>/pred.parquet) of the latest ok
registry rows named in configs/mutation_dependencies.json, dev targets through essml.guard, and
computes the metrics of essml.mutdeps per pair, model and outer fold (configs/splits.csv
cv_fold), the biomarker-only reference, the "clearly helps" rule and a mechanism check.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

for _v in ["OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"]:
    os.environ.setdefault(_v, "4")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
from essml import features as F  # noqa: E402
from essml import mutdeps as MD  # noqa: E402

PAIR_CSV = ROOT / "configs" / "oncogene_pairs.csv"
DECL = ROOT / "configs" / "mutation_dependencies.json"
OUT_JSON = ROOT / "results" / "mutation_dependencies.json"
FIG = ROOT / "docs" / "figures" / "mutation_dependencies.png"


def dev_ids():
    sp = pd.read_csv(ROOT / "configs" / "splits.csv", index_col=0)
    return sp[sp["role"] == "dev"]


def raw_status_mats():
    from essml.io import load_matrix
    return load_matrix("hotspot"), load_matrix("damaging"), load_matrix("cn")


def entrez_of(symbols):
    from essml.io import load_gene_map
    gm = load_gene_map()
    m = dict(zip(gm["symbol"], gm["entrez"].astype(str)))
    missing = [s for s in symbols if s not in m]
    if missing:
        raise KeyError(f"symbols not in gene map: {missing}")
    return {s: m[s] for s in symbols}


def threshold_text(btype):
    t = {"hotspot": "hotspot matrix > 0", "damaging": "damaging dosage > 0",
         "amp": f"relative CN > {F.AMP_THR}", "deep_del": f"relative CN < {F.DEEP_DEL_THR}",
         "TP53-WT": "no TP53 damaging and no TP53 hotspot (mutant = TP53 wild type)"}
    return " or ".join(t[p] for p in btype.split("|"))


def pair_table(min_mut=5):
    dev = dev_ids()
    hot, dam, cn = raw_status_mats()
    syms = sorted({s for p in MD.PAIRS for s in (p[2], p[3])})
    ent = entrez_of(syms)
    rows = []
    for pid, group, dep, bm, btype in MD.PAIRS:
        s = MD.biomarker_status(btype, ent[bm], dev.index, hot, dam, cn)
        n_mut = int((s == 1).sum())
        rows.append({
            "pair_id": pid, "group": group, "dependency_gene": dep,
            "dependency_entrez": ent[dep], "biomarker_gene": bm, "biomarker_entrez": ent[bm],
            "biomarker_type": btype, "threshold": threshold_text(btype),
            "amp_thr": F.AMP_THR if "amp" in btype else "",
            "deep_del_thr": F.DEEP_DEL_THR if "deep_del" in btype else "",
            "n_mutant_dev": n_mut, "n_wt_dev": int((s == 0).sum()),
            "n_unknown_dev": int(s.isna().sum()), "eligible": n_mut >= min_mut})
    return pd.DataFrame(rows)


# ============================================================================ run
def load_preds(run_ids, genes, ids):
    out = {}
    for name, rid in run_ids.items():
        p = pd.read_parquet(ROOT / "cache" / "runs" / rid / "pred.parquet", columns=list(genes))
        out[name] = p.reindex(ids).astype(np.float64)
        if out[name].isna().any().any():
            raise ValueError(f"{name}: missing OOF predictions")
    return out


def latest_ok_ids(names):
    from g2_report import latest_ok
    reg = pd.read_csv(ROOT / "results" / "registry.csv", dtype=str)
    ids = {n: latest_ok(reg, n) for n in names}
    miss = [n for n, r in ids.items() if r is None or not
            (ROOT / "cache" / "runs" / r / "pred.parquet").exists()]
    if miss:
        raise FileNotFoundError(f"no ok run with predictions for {miss}; re-run via scripts/run.py")
    return ids


def compute(decl, pairs, preds, Y, status, folds, models):
    """Per pair x model x fold metrics (long DataFrame); biomarker-only preds added to ``preds``."""
    rows = []
    ids = Y.index
    bo = {}
    for _, pr in pairs.iterrows():
        g, pid = pr["dependency_entrez"], pr["pair_id"]
        s = status[pid]
        b = pd.Series(np.nan, index=ids)
        for k in sorted(folds.unique()):
            tr, te = folds.index[folds != k], folds.index[folds == k]
            b.loc[te] = MD.biomarker_only(Y.loc[tr, g], s.loc[tr], s.loc[te])
        bo[pid] = b
        for m in models:
            P = b if m == "biomarker_only" else preds[m][g]
            for k in sorted(folds.unique()):
                te = folds.index[folds == k]
                r = MD.fold_metrics(P.loc[te].to_numpy(), Y.loc[te, g].to_numpy(),
                                    s.loc[te].to_numpy(), **decl["metric_params"])
                rows.append({"pair_id": pid, "model": m, "fold": int(k), **r})
    return pd.DataFrame(rows), bo


def summarise(long, pairs, models, ref, comparators, decl):
    """Per pair x model: mean over folds, n folds defined, paired delta vs ref with CI."""
    mets = ["r_all", "r_mut", "r_wt", "ap_all", "nap_all", "ap_mut", "nap_mut", "ap_wt", "nap_wt",
            "auroc_bm"]
    out = {}
    for pid in pairs["pair_id"]:
        sub = long[long["pair_id"] == pid]
        out[pid] = {}
        base = sub[sub["model"] == ref].set_index("fold").sort_index()
        for m in models:
            x = sub[sub["model"] == m].set_index("fold").sort_index()
            d = {"n": {c: int(x[c].iloc[0]) for c in ["n_obs", "n_dep"]},
                 "n_fold_mut": x["n_mut"].tolist(), "n_fold_dep_mut": x["n_dep_mut"].tolist()}
            for c in mets:
                v = x[c].to_numpy(np.float64)
                d[c] = {"mean": float(np.nanmean(v)) if np.isfinite(v).any() else None,
                        "n_folds": int(np.isfinite(v).sum()), "folds": np.round(v, 5).tolist()}
                if m != ref:
                    dd = v - base[c].to_numpy(np.float64)
                    ci = MD.paired_ci(dd)
                    d[c]["delta_vs_ref"] = {**ci, "folds": np.round(dd, 5).tolist()}
            out[pid][m] = d
    return out


def pooled_ci(long, eligible, model, ref, metric):
    """Fold-wise mean over eligible pairs defined in all folds of metric(model) - metric(ref)."""
    a = long[long["model"] == model].pivot(index="pair_id", columns="fold", values=metric)
    b = long[long["model"] == ref].pivot(index="pair_id", columns="fold", values=metric)
    d = (a - b).reindex(eligible).dropna()
    ci = MD.paired_ci(d.mean(axis=0).to_numpy())
    ci["pairs"] = d.index.tolist()
    return ci


def mechanism(pairs, preds, Y, status, lineage, E, models):
    """Pooled-OOF descriptive check: what the expression-only prediction tracks (dev lines)."""
    out = {}
    lin = lineage.reindex(Y.index).astype(str)

    def corr(a, b):
        m = np.isfinite(a) & np.isfinite(b)
        return float(np.corrcoef(a[m], b[m])[0, 1]) if m.sum() > 2 else None

    def eta2(v):
        v = pd.Series(v, index=Y.index)
        fit = v.groupby(lin).transform("mean")
        return float(1 - ((v - fit) ** 2).sum() / ((v - v.mean()) ** 2).sum())

    def within_lin_r(p, t):
        d = pd.DataFrame({"p": p, "t": t, "l": lin.to_numpy()}).dropna()
        d["p"] -= d.groupby("l")["p"].transform("mean")
        d["t"] -= d.groupby("l")["t"].transform("mean")
        return float(np.corrcoef(d["p"], d["t"])[0, 1])

    for _, pr in pairs.iterrows():
        g, bm, pid = pr["dependency_entrez"], pr["biomarker_entrez"], pr["pair_id"]
        e_own = E[g].reindex(Y.index).to_numpy(np.float64) if g in E.columns else None
        e_bm = E[bm].reindex(Y.index).to_numpy(np.float64) if bm in E.columns else None
        t = Y[g].to_numpy(np.float64)
        s = status[pid].reindex(Y.index).to_numpy(np.float64)
        d = {"r_true_own_expr": corr(t, e_own) if e_own is not None else None,
             "r_true_status": corr(t, s), "eta2_true_lineage": eta2(np.where(np.isnan(t), np.nanmean(t), t)),
             "r_status_own_expr": corr(s, e_own) if e_own is not None else None,
             "r_status_biomarker_expr": corr(s, e_bm) if e_bm is not None else None}
        for m in models:
            p = preds[m][g].to_numpy(np.float64)
            d[m] = {"r_pred_own_expr": corr(p, e_own) if e_own is not None else None,
                    "r_pred_biomarker_expr": corr(p, e_bm) if e_bm is not None else None,
                    "r_pred_status": corr(p, s), "eta2_pred_lineage": eta2(p),
                    "r_pooled": corr(p, t), "r_within_lineage": within_lin_r(p, t)}
        out[pid] = d
    return out


def figure(pairs_show, pairs, preds, Y, status, path, models):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    nr, nc = len(models), len(pairs_show)
    fig, axes = plt.subplots(nr, nc, figsize=(2.9 * nc, 2.8 * nr), dpi=110, squeeze=False)
    pt = pairs.set_index("pair_id")
    for i, model in enumerate(models):
        for j, pid in enumerate(pairs_show):
            ax = axes[i, j]
            g = pt.loc[pid, "dependency_entrez"]
            p, t = preds[model][g], Y[g]
            s = status[pid].reindex(Y.index)
            for lab, m, col, sz in (("WT", s == 0, "#9aa5b1", 6), ("unknown", s.isna(), "#d9b26f", 6),
                                    ("mutant", s == 1, "#c0392b", 10)):
                if m.sum():
                    ax.scatter(p[m], t[m], s=sz, c=col, alpha=0.75, lw=0, label=f"{lab} ({int(m.sum())})")
            ax.axhline(MD.THR, color="0.6", lw=0.6, ls="--")
            if i == 0:
                ax.set_title(f"{pt.loc[pid, 'dependency_gene']} | {pt.loc[pid, 'biomarker_gene']} "
                             f"{pt.loc[pid, 'biomarker_type']}", fontsize=8.5)
            ax.set_xlabel(f"predicted, {model} (OOF)", fontsize=7.5)
            if j == 0:
                ax.set_ylabel("observed Chronos effect", fontsize=7.5)
            ax.tick_params(labelsize=7)
            if i == 0:
                ax.legend(fontsize=6, frameon=False, loc="upper left", markerscale=1.5)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)


def run():
    from essml import guard
    decl = json.loads(DECL.read_text())
    pairs = pd.read_csv(PAIR_CSV, dtype={"dependency_entrez": str, "biomarker_entrez": str})
    dev = dev_ids()
    folds = dev["cv_fold"].astype(int)
    hot, dam, cn = raw_status_mats()
    status = {r["pair_id"]: MD.biomarker_status(r["biomarker_type"], r["biomarker_entrez"],
                                                dev.index, hot, dam, cn)
              for _, r in pairs.iterrows()}
    # the committed pair file must match the status code exactly
    for _, r in pairs.iterrows():
        if int((status[r["pair_id"]] == 1).sum()) != int(r["n_mutant_dev"]):
            raise AssertionError(f"{r['pair_id']}: n mutant differs from {PAIR_CSV.name}")
    genes = pairs["dependency_entrez"].unique().tolist()
    run_ids = latest_ok_ids(decl["models"])
    preds = load_preds(run_ids, genes, dev.index)
    Y = guard.load_targets(dev.index)[genes].astype(np.float64)  # dev IDs only
    models = decl["models"] + ["biomarker_only"]
    long, bo = compute(decl, pairs, preds, Y, status, folds, models)
    ref = decl["reference"]
    eligible = pairs.loc[pairs["eligible"], "pair_id"].tolist()
    per_pair = summarise(long, pairs, models, ref, decl["comparators"], decl)
    rule = {}
    for c in decl["comparators"]:
        a = long[long["model"] == c].pivot(index="pair_id", columns="fold", values="ap_all")
        b = long[long["model"] == ref].pivot(index="pair_id", columns="fold", values="ap_all")
        dap = {p: (a.loc[p] - b.loc[p]).to_numpy(np.float64) for p in a.index}
        rule[c] = MD.clearly_helps(dap, eligible, **decl["rule"])
    rule_met = any(v["met"] for v in rule.values())
    pooled = {m: {met: pooled_ci(long, eligible, m, ref, met)
                  for met in ["ap_all", "nap_all", "r_all", "ap_mut", "r_mut", "auroc_bm"]}
              for m in models if m != ref}
    from essml.io import load_matrix
    E = load_matrix("expression")
    lineage = dev["OncotreeLineage"]
    mech = mechanism(pairs, {m: preds[m] for m in decl["models"]}, Y, status, lineage, E,
                     decl["models"])
    figure(decl["figure_pairs"], pairs, preds, Y, status, FIG, decl["figure_models"])
    res = {"design": "DESIGN.md §4 'Mutation-driven dependencies'; RESEARCH_TREE D24",
           "declared": str(DECL.relative_to(ROOT)), "pair_file": str(PAIR_CSV.relative_to(ROOT)),
           "run_ids": run_ids, "reference": ref, "eligible_pairs": eligible,
           "rule": rule, "rule_met": rule_met, "pooled_vs_ref": pooled,
           "per_pair": per_pair, "mechanism": mech,
           "pairs": pairs.to_dict(orient="records")}
    OUT_JSON.write_text(json.dumps(res, indent=1, default=_json))
    long.to_csv(ROOT / "results" / "mutation_dependencies_folds.csv", index=False,
                float_format="%.6g")
    print(json.dumps({c: {k: v[k] for k in ["met", "cond1_ci_excludes_0", "cond2_ge3_pairs",
                                             "pairs_gain"]} | {"pooled": v["pooled"]}
                      for c, v in rule.items()}, indent=1, default=_json))
    return res


def _json(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("cmd", choices=["pairs", "run"])
    a = ap.parse_args(argv)
    if a.cmd == "pairs":
        t = pair_table()
        t.to_csv(PAIR_CSV, index=False)
        print(t[["pair_id", "biomarker_type", "n_mutant_dev", "n_wt_dev", "n_unknown_dev",
                 "eligible"]].to_string(index=False))
    else:
        run()


if __name__ == "__main__":
    main()
