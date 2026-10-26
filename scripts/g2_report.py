"""G2 report: results table, paired contrasts, Holm-adjusted ablation, selection rule (DESIGN §4, §6).

    python scripts/g2_report.py            # -> docs/g2_results.md, results/g2_summary.json

Uses the latest registry row with status ok for each run name. All numbers come from the frozen
evaluation (results/runs/<id>/per_fold.csv, written by scripts/run.py via essml.evaluate) except
the common-line-subset re-evaluation, which calls the same frozen ``essml.evaluate.evaluate`` on
the cached out-of-fold predictions restricted to the subset lines.

Paired fold-wise contrasts: d_k = P1_k(A) - P1_k(B) over the 5 outer folds (same lines and genes);
SE = sd(d)/sqrt(5); two-sided paired t-test (df 4). Holm adjustment over the declared ablation
contrasts. Selection rule (DESIGN §4): in complexity order a candidate replaces the current one
only if mean d > 0.005 and mean d > 1 SE.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DELTA = 0.005

RUNS = ["b1", "b2", "r1_E", "m1_E", "m1b_E_C", "m1b_E_M", "m1b_E_H", "m1b_E_F", "m1b_E_Meta",
        "m1b_E_comp", "m1_E_Psub", "m1b_E_P_Psub", "m1b_full", "m1b_full_minus_E",
        "m1b_full_minus_C", "m1b_full_minus_M", "m1b_full_minus_H", "m1b_full_minus_F",
        "m1b_full_minus_Meta", "m1b_full_minus_comp", "m1b_Efree", "m1b_full_prov", "t_E_status",
        "m1_E_permuted", "m1_E_resB0"]
DIAGNOSTIC = [("m1_E_resB0", "m1_E", "B0 vs B1 residual target (E only)"),
              ("r1_E", "m1_E", "R1 (raw y, single lambda) vs M1-E"),
              ("t_E_status", "m1_E", "T (E+status, B0 residual) vs M1-E"),
              ("t_E_status", "m1_E_resB0", "status flags given B0 residual (T vs M1-E-resB0)"),
              ("m1_E_resB0", "r1_E", "grouped vs single lambda (B0 residual)"),
              ("m1_E_permuted", "b1", "negative control: permuted features vs B1"),
              ("m1_E", "b2", "M1-E vs B2 (kNN)")]
FULL = ["E", "C", "M", "H", "F", "Meta", "composites"]
TAG = {"composites": "comp"}


def contrasts():
    """Declared ablation contrasts (A, B, label): effect = P1(A) - P1(B)."""
    c = [(f"m1b_E_{TAG.get(m, m)}", "m1_E", f"add {m} to E") for m in ["C", "M", "H", "F", "Meta",
                                                                       "composites"]]
    c.append(("m1b_E_P_Psub", "m1_E_Psub", "add P to E (P subset)"))
    c += [("m1b_full", f"m1b_full_minus_{TAG.get(m, m)}", f"{m} in full (leave-one-out)") for m in FULL]
    c.append(("m1b_full", "m1b_Efree", "E vs E-free (full - E-free)"))
    c.append(("m1b_full_prov", "m1b_full", "add provenance Meta fields to full"))
    return c


def latest_ok(reg, name):
    r = reg[(reg["name"] == name) & (reg["status"] == "ok")]
    return None if r.empty else r.iloc[-1]["run_id"]


def load_run(rid):
    d = ROOT / "results" / "runs" / rid
    s = json.loads((d / "summary.json").read_text())
    pf = pd.read_csv(d / "per_fold.csv")
    return s, pf


def fold_vals(pf, cls, col):
    sub = pf[pf["class"] == cls].sort_values("fold")
    return sub[col].to_numpy(np.float64)


def paired(a, b):
    d = np.asarray(a) - np.asarray(b)
    n = len(d)
    m, sd = float(d.mean()), float(d.std(ddof=1))
    se = sd / np.sqrt(n)
    if se > 0:
        t = m / se
        p = float(2 * stats.t.sf(abs(t), n - 1))
    else:
        t, p = (np.inf if m != 0 else 0.0), (0.0 if m != 0 else 1.0)
    h = stats.t.ppf(0.975, n - 1) * se
    return {"mean": m, "sd": sd, "se": se, "t": float(t), "p": p, "ci95": [m - h, m + h],
            "n_folds_pos": int((d > 0).sum()), "folds": d.round(5).tolist()}


def holm(ps):
    ps = np.asarray(ps, float)
    o = np.argsort(ps)
    adj = np.empty_like(ps)
    run = 0.0
    for i, j in enumerate(o):
        run = max(run, (len(ps) - i) * ps[j])
        adj[j] = min(1.0, run)
    return adj


def subset_eval(run_ids, line_block="C", min_obs=40):
    """Per-fold P1 of cached OOF predictions on dev lines having raw block ``line_block``.

    Same definition as essml.evaluate (frozen): within each outer fold, mean of the frozen
    ``metrics.per_gene_pearson`` over class genes in configs/eval_genes.csv with >= min_obs
    observed held-out lines; only P1 is needed here, so the full evaluate() is not re-run.
    """
    from essml import genesets, guard, metrics
    sp = pd.read_csv(ROOT / "configs" / "splits.csv", index_col=0)
    dev = sp[sp["role"] == "dev"]
    have = pd.read_parquet(ROOT / "cache" / "features" / f"{line_block}.parquet", columns=[]).index
    sub = dev.index[dev.index.isin(have)]
    genes = pd.Index(pd.read_csv(ROOT / "configs" / "eval_genes.csv", dtype={"entrez": str})["entrez"])
    gc = genesets.load_gene_classes().reindex(genes).fillna(False).astype(bool)
    Y = guard.load_targets(sub).astype(np.float64)[genes]
    fold = dev.loc[sub, "cv_fold"].astype(int)
    out = {}
    for name, rid in run_ids.items():
        pred = pd.read_parquet(ROOT / "cache" / "runs" / rid / "pred.parquet")
        if not sub.isin(pred.index).all():
            continue
        res = {"n_lines": len(sub)}
        for cls in ["variable", "all"]:
            vals = []
            for f in sorted(fold.unique()):
                ids = fold.index[fold == f]
                T = Y.loc[ids].to_numpy(np.float64)
                P = pred.loc[ids, genes].to_numpy(np.float64)
                r = metrics.per_gene_pearson(P, T, min_obs)
                m = gc[cls].to_numpy() & ((~np.isnan(T)).sum(0) >= min_obs)
                rr = r[m]
                vals.append(float(rr[np.isfinite(rr)].mean()))
            res[f"P1_{cls}_folds"] = vals
        out[name] = res
    return out


def fmt(x, d=3):
    return "—" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{d}f}"


def main(argv=None):
    reg = pd.read_csv(ROOT / "results" / "registry.csv", dtype=str)
    ids = {n: latest_ok(reg, n) for n in RUNS}
    ids.update({n: latest_ok(reg, n) for n in reg["name"].unique() if str(n).startswith("m2_")})
    missing = [n for n, r in ids.items() if r is None]
    ids = {n: r for n, r in ids.items() if r is not None}
    runs = {n: load_run(r) for n, r in ids.items()}
    b2pf = runs["b2"][1] if "b2" in runs else None
    rows = []
    for n, (s, pf) in runs.items():
        sm = s["summary"]
        p1v = fold_vals(pf, "variable", "P1")
        row = {"name": n, "run_id": ids[n], "model": s["config"]["model"],
               "modalities": "+".join(s["config"].get("params", {}).get("modalities", [])
                                      or s["config"].get("params", {}).get("m1_params", {}).get("modalities", [])
                                      or s["config"].get("features", [])),
               "line_subset": s["config"].get("line_subset", ""), "n_lines": s["n_lines"],
               "P1_var": sm["variable"]["P1"], "P1_var_sd": sm["variable"]["P1_sd"],
               "P1_all": sm["all"]["P1"], "P1_sel": sm["selective"]["P1"],
               "P1_rare": sm["rare_selective"]["P1"],
               "P2a": sm["all"]["P2a"], "dP2a_B0": sm["all"].get("dP2a_vs_B0"),
               "P2b_nap": sm["selective"]["P2b_nap_prev05"],
               "P2b_auroc": sm["selective"]["P2b_auroc_prev05"],
               "fit_s": s["fit_s"], "predict_s": s["predict_s"], "peak_rss_mb": s["peak_rss_mb"],
               "threads": s["threads"]}
        if "dP1_vs_B1" in pf.columns:
            d = paired(p1v, p1v - fold_vals(pf, "variable", "dP1_vs_B1"))
            row.update({"dP1_B1": d["mean"], "dP1_B1_se": d["se"], "dP1_B1_sd": d["sd"]})
        if b2pf is not None and not row["line_subset"] and n != "b2":
            d = paired(p1v, fold_vals(b2pf, "variable", "P1"))
            row.update({"dP1_B2": d["mean"], "dP1_B2_se": d["se"]})
        rep = s.get("model_report", {})
        if rep:
            r0 = next(iter(rep.values()))
            row["lam"] = json.dumps({k: sorted({float(v["lam_by_group"][k]) for v in rep.values()})
                                     for k in r0.get("lam_by_group", {})})
            row["lam_edge"] = any(any(v.get("lam_at_grid_edge", {}).values()) for v in rep.values())
            row["n_features"] = r0.get("n_features") or r0.get("m1", {}).get("n_features")
        rows.append(row)
    tab = pd.DataFrame(rows)

    # ---- ablation contrasts (Holm)
    con = []
    for a, b, lab in contrasts():
        if a in runs and b in runs:
            d = paired(fold_vals(runs[a][1], "variable", "P1"), fold_vals(runs[b][1], "variable", "P1"))
            con.append({"contrast": lab, "A": a, "B": b, **d})
    if con:
        adj = holm([c["p"] for c in con])
        for c, pa in zip(con, adj):
            c["p_holm"] = float(pa)
            c["helps"] = bool(c["mean"] > 0 and pa < 0.05)

    diag = []
    for a, b, lab in DIAGNOSTIC:
        if a in runs and b in runs:
            d = paired(fold_vals(runs[a][1], "variable", "P1"), fold_vals(runs[b][1], "variable", "P1"))
            diag.append({"contrast": lab, "A": a, "B": b, **d})

    # ---- common-subset rule for C (967/1022 lines): re-evaluate on C lines
    c_runs = [n for n in ["m1_E", "m1b_E_C", "m1b_full", "m1b_full_minus_C", "m1b_Efree", "b1", "b2"]
              if n in ids]
    csub = subset_eval({n: ids[n] for n in c_runs}, "C") if c_runs else {}
    c_con = []
    for a, b, lab in [("m1b_E_C", "m1_E", "add C to E"), ("m1b_full", "m1b_full_minus_C", "C in full"),
                      ("m1b_E_C", "b1", "E+C vs B1"), ("m1_E", "b1", "E vs B1")]:
        if a in csub and b in csub:
            d = paired(csub[a]["P1_variable_folds"], csub[b]["P1_variable_folds"])
            c_con.append({"contrast": lab + " (C-covered lines only)", "A": a, "B": b, **d})

    # ---- selection rule (complexity order; M1 rungs on all 1022 lines)
    def nmod(n):
        return len(runs[n][0]["config"]["params"].get("modalities", []))
    cand = [n for n in runs if n.startswith(("m1_", "m1b_")) and not runs[n][0]["config"].get("line_subset")
            and "permute_seed" not in runs[n][0]["config"]["params"]
            and "residual" not in runs[n][0]["config"]["params"]]  # diagnostics are not candidates
    cand = sorted(cand, key=lambda n: (nmod(n), runs[n][0]["summary"]["variable"]["P1"] * 0 +
                                       (runs[n][0]["model_report"]["0"]["n_features"]
                                        if runs[n][0].get("model_report") else 0)))
    sel_log = []
    current = "b1"
    for n in ["m1_E"] + [c for c in cand if c != "m1_E"]:
        if n not in runs:
            continue
        d = paired(fold_vals(runs[n][1], "variable", "P1"), fold_vals(runs[current][1], "variable", "P1"))
        win = d["mean"] > DELTA and d["mean"] > d["se"]
        sel_log.append({"candidate": n, "vs": current, "dP1": d["mean"], "se": d["se"], "replaces": win})
        if win:
            current = n
    m2 = sorted([n for n in runs if n.startswith("m2_")])
    m2_log = []
    m1_sel = current
    for n in m2:
        d = paired(fold_vals(runs[n][1], "variable", "P1"), fold_vals(runs[m1_sel][1], "variable", "P1"))
        win = d["mean"] > DELTA and d["mean"] > d["se"]
        m2_log.append({"candidate": n, "vs": m1_sel, "dP1": d["mean"], "se": d["se"], "beats_M1": win})
    best_m2 = None
    winners = [x for x in m2_log if x["beats_M1"]]
    if winners:  # simplest prior set that beats M1 (none < paralogs < partners), then same rule
        order = {"none": 0, "paralogs": 1, "partners": 2}
        winners.sort(key=lambda x: order[x["candidate"].rsplit("_", 1)[1]])
        best_m2 = winners[0]["candidate"]
        for x in winners[1:]:
            d = paired(fold_vals(runs[x["candidate"]][1], "variable", "P1"),
                       fold_vals(runs[best_m2][1], "variable", "P1"))
            if d["mean"] > DELTA and d["mean"] > d["se"]:
                best_m2 = x["candidate"]
    m2_prior = []
    for a, b in [("partners", "paralogs"), ("paralogs", "none")]:
        pa = [n for n in m2 if n.endswith("_" + a)]
        pb = [n for n in m2 if n.endswith("_" + b)]
        if pa and pb:
            d = paired(fold_vals(runs[pa[0]][1], "variable", "P1"), fold_vals(runs[pb[0]][1], "variable", "P1"))
            m2_prior.append({"contrast": f"M2 {a} vs {b}", "A": pa[0], "B": pb[0], **d})
    if m2_prior:
        for c, pa in zip(m2_prior, holm([c["p"] for c in m2_prior])):
            c["p_holm"] = float(pa)

    out = {"runs": tab.to_dict("records"), "missing_runs": missing, "contrasts": con,
           "diagnostic_contrasts": diag,
           "c_subset": {"runs": csub, "contrasts": c_con}, "selection_log": sel_log,
           "selected_m1": current, "m2_log": m2_log, "selected_m2": best_m2,
           "m2_prior_contrasts": m2_prior,
           "selected_model": best_m2 or current, "delta": DELTA}
    (ROOT / "results" / "g2_summary.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    write_md(out, tab)
    print(json.dumps({"selected_m1": current, "selected_m2": best_m2, "missing": missing}, indent=1))


def write_md(out, tab):
    L = ["# G2 results — M1/M1b kernel ridge, R1, T, modality ablation, M2 (dev CV)", "",
         "Generated by `scripts/g2_report.py` from `results/registry.csv` (latest ok run per name).",
         "All metrics are computed within each of the 5 outer dev folds and averaged (± fold SD);",
         "Δ columns are paired fold-wise differences (mean, SE = SD/√5). P1 = mean per-gene Pearson r",
         "on variable genes (dev-wide `configs/gene_classes.csv`). ΔP1 vs B2 uses the B2 run's",
         "per-fold P1 on the same folds/genes (not available for subset runs). Runtime = sum over",
         "the 5 folds of fit (incl. inner tuning) and predict wall-clock; peak RSS of the run process.",
         "", "## Results table", "",
         "| run | model | features | lines | P1 var | ΔP1 vs B1 | ΔP1 vs B2 | P1 all | P1 sel | "
         "P2a | ΔP2a vs B0 | P2b nAP | P2b AUROC | fit s | pred s | RSS MB | λ (other/var) |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in out["runs"]:
        def pm(k, se=None):
            v = r.get(k)
            if v is None or (isinstance(v, float) and not np.isfinite(v)):
                return "—"
            return f"{v:+.4f} ± {r[se]:.4f}" if se and r.get(se) is not None else f"{v:.4f}"
        L.append(f"| {r['name']} | {r['model']} | {r['modalities'] or '—'} | {r['n_lines']} | "
                 f"{r['P1_var']:.4f} ± {r['P1_var_sd']:.4f} | {pm('dP1_B1', 'dP1_B1_se')} | "
                 f"{pm('dP1_B2', 'dP1_B2_se')} | {r['P1_all']:.4f} | {r['P1_sel']:.4f} | "
                 f"{r['P2a']:.4f} | {fmt(r.get('dP2a_B0'), 4)} | {r['P2b_nap']:.4f} | "
                 f"{r['P2b_auroc']:.4f} | {r['fit_s']:.0f} | {r['predict_s']:.1f} | "
                 f"{r['peak_rss_mb']:.0f} | {r.get('lam', '—') if isinstance(r.get('lam'), str) else '—'} |")
    L += ["", "## Ablation contrasts (P1 variable; Holm over the declared contrasts)", "",
          "| contrast | A | B | ΔP1 | SE | 95% CI (unadj.) | folds A>B | p | p Holm | helps |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for c in out["contrasts"]:
        L.append(f"| {c['contrast']} | {c['A']} | {c['B']} | {c['mean']:+.4f} | {c['se']:.4f} | "
                 f"[{c['ci95'][0]:+.4f}, {c['ci95'][1]:+.4f}] | {c['n_folds_pos']}/5 | {c['p']:.3g} | "
                 f"{c['p_holm']:.3g} | {'yes' if c['helps'] else 'no'} |")
    L += ["", "## Diagnostic contrasts (not in the Holm family; not selection candidates)", "",
          "| contrast | A | B | ΔP1 | SE | folds A>B | p |", "|---|---|---|---|---|---|---|"]
    for c in out["diagnostic_contrasts"]:
        L.append(f"| {c['contrast']} | {c['A']} | {c['B']} | {c['mean']:+.4f} | {c['se']:.4f} | "
                 f"{c['n_folds_pos']}/5 | {c['p']:.3g} |")
    L += ["", "## Common-subset rule for C (lines with CN data only)", "",
          "Default runs impute missing C as 0 (diploid) plus a `C:absent` flag on all 1022 lines (table above).",
          "Here the same out-of-fold predictions are re-evaluated on the C-covered dev lines only.", "",
          "| contrast | ΔP1 | SE | folds A>B | p |", "|---|---|---|---|---|"]
    for c in out["c_subset"]["contrasts"]:
        L.append(f"| {c['contrast']} | {c['mean']:+.4f} | {c['se']:.4f} | {c['n_folds_pos']}/5 | {c['p']:.3g} |")
    for n, v in out["c_subset"]["runs"].items():
        L.append(f"- {n}: P1 var on {v['n_lines']} C lines = {np.mean(v['P1_variable_folds']):.4f}")
    L += ["", "## Selection rule (DESIGN §4: ΔP1 > 0.005 and > 1 SE, complexity order)", "",
          "| candidate | vs current | ΔP1 | SE | replaces |", "|---|---|---|---|---|"]
    for s in out["selection_log"]:
        L.append(f"| {s['candidate']} | {s['vs']} | {s['dP1']:+.4f} | {s['se']:.4f} | {'yes' if s['replaces'] else 'no'} |")
    L += ["", f"**Selected M1/M1b configuration: `{out['selected_m1']}`.**", ""]
    if out["m2_log"]:
        L += ["## M2 local stage", "", "| candidate | vs | ΔP1 | SE | beats M1 by rule |", "|---|---|---|---|---|"]
        for s in out["m2_log"]:
            L.append(f"| {s['candidate']} | {s['vs']} | {s['dP1']:+.4f} | {s['se']:.4f} | {'yes' if s['beats_M1'] else 'no'} |")
        for c in out["m2_prior_contrasts"]:
            L.append(f"- {c['contrast']}: ΔP1 {c['mean']:+.4f} ± {c['se']:.4f} (Holm p {c['p_holm']:.3g})")
        L += ["", f"**Selected model after G2: `{out['selected_model']}`.**", ""]
    if out["missing_runs"]:
        L += [f"Missing (no ok run): {', '.join(out['missing_runs'])}", ""]
    md = ROOT / "docs" / "g2_results.md"
    notes = ""
    if md.exists():  # keep the hand-written interpretation between the markers
        t = md.read_text(encoding="utf-8")
        a, b = t.find("<!-- notes -->"), t.find("<!-- /notes -->")
        if a >= 0 and b > a:
            notes = t[a:b + len("<!-- /notes -->")]
    L[1:1] = ["", notes or "<!-- notes -->\n<!-- /notes -->"]
    md.write_text("\n".join(L) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
