"""G3 close-out report: M5 stacking (DESIGN §5) and the D20c lambda-criterion diagnostic.

    python scripts/g3_m5_report.py   # -> docs/g3_m5_results.md, results/g3_m5_summary.json

Latest registry row with status ok per run name. P1 / P2b / P1 all / runtime come from the frozen
evaluation outputs (results/runs/<id>/summary.json, per_fold.csv). Within-library, within-lineage
and within lineage x library P1 and the Nadeau-Bengio SE use the functions of
scripts/g3_kernel_report.py on the cached out-of-fold predictions (cache/runs/<id>/pred.parquet),
dev targets through ``essml.guard``. P1 by gene-SD decile: per outer fold, deciles of the gene SD
over that fold's TRAINING lines (all eval genes), mean per-gene r (frozen per_gene_fold.parquet)
within each decile, then averaged over folds.

Pre-registered rules (configs/runs/m5_*.json, m1_E_resB0_*.json):
* M5 rung representative = m5_stack if paired dP1(m5_stack - m5_avg) > 0 and > 1 fold SE, else
  m5_avg; it replaces m1_E_resB0 only if paired dP1 vs m1_E_resB0 > 0.005 and > 1 fold SE.
* lambda diagnostic: a criterion run replaces the headline under the same rule (D20c, like D14).
* control m5_ctl_perm: weight on the permuted input <= 0.01 in every class and fold, and P1
  variable >= m5_stack - 0.002.
"""
from __future__ import annotations

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
from g2_report import fold_vals, latest_ok, load_run, paired  # noqa: E402
from g3_kernel_report import (dev_frame, load_targets_eval, paired_nb, test_train_ratio,  # noqa: E402
                              within_p1)

DELTA = 0.005
HEADLINE = "m1_E_resB0"
INPUTS = [HEADLINE, "m2n_E_paralogs", "g3_m3n"]
M5 = ["m5_avg", "m5_stack", "m5_stack_f1", "m5_ctl_perm"]
LAM = ["m1_E_resB0_loomse", "m1_E_resB0_icvr"]
NAMES = ["b1"] + INPUTS + M5 + LAM
MD = ROOT / "docs" / "g3_m5_results.md"
OUT_JSON = ROOT / "results" / "g3_m5_summary.json"
PERM_INPUT = "m2n_E_paralogs"


def f4(x, sign=False):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "—"
    return f"{x:+.4f}" if sign else f"{x:.4f}"


def sd_deciles(Y, folds, run_ids, genes):
    """Mean per-gene r by training-fold gene-SD decile (all eval genes), per fold then averaged."""
    out = {}
    pgf = {n: pd.read_parquet(ROOT / "cache" / "runs" / rid / "per_gene_fold.parquet")
           for n, rid in run_ids.items()}
    edges_rep = []
    for f in sorted(folds.unique()):
        tr = folds.index[folds != f]
        sd = Y.loc[tr, genes].std(ddof=1)
        q = np.nanquantile(sd.to_numpy(), np.linspace(0, 1, 11))
        dec = pd.Series(np.clip(np.searchsorted(q, sd.to_numpy(), side="right") - 1, 0, 9),
                        index=genes)
        edges_rep.append(q.round(4).tolist())
        for n, t in pgf.items():
            r = t[t["fold"] == f].set_index("gene")["r"].reindex(genes)
            m = r.groupby(dec).mean()
            for d, v in m.items():
                out.setdefault(n, {}).setdefault(int(d), []).append(float(v))
    return out, edges_rep


def main(argv=None):
    reg = pd.read_csv(ROOT / "results" / "registry.csv", dtype=str)
    ids = {n: latest_ok(reg, n) for n in NAMES}
    missing = [n for n, r in ids.items() if r is None]
    ids = {n: r for n, r in ids.items() if r is not None}
    runs = {n: load_run(r) for n, r in ids.items()}
    regrow = {n: reg[reg["run_id"] == r].iloc[-1] for n, r in ids.items()}
    dev, lib = dev_frame()
    folds = dev["cv_fold"].astype(int)
    lineage = dev["OncotreeLineage"].astype(str)
    Y, genes = load_targets_eval(dev)
    from essml import genesets
    gc = genesets.load_gene_classes().reindex(genes).fillna(False).astype(bool)
    var_mask = gc["variable"].to_numpy()
    groupings = {"within_library": lib, "within_lineage": lineage,
                 "within_lin_lib": lineage + "|" + lib.astype(int).astype(str)}
    ratio = test_train_ratio(folds)
    head_pf = runs[HEADLINE][1]
    rows, within, preds = [], {}, {}
    for n in [x for x in NAMES if x in runs]:
        s, pf = runs[n]
        sm = s["summary"]
        p1v = fold_vals(pf, "variable", "P1")
        pred = pd.read_parquet(ROOT / "cache" / "runs" / ids[n] / "pred.parquet")
        preds[n] = pred
        w = within_p1(pred, Y, folds, groupings, var_mask)
        if not np.allclose(w["raw"], p1v, atol=1e-6):
            raise AssertionError(f"{n}: recomputed P1 {w['raw']} != per_fold {p1v.tolist()}")
        within[n] = w
        row = {"name": n, "run_id": ids[n], "P1_var": sm["variable"]["P1"],
               "P1_var_sd": sm["variable"]["P1_sd"], "P1_all": sm["all"]["P1"],
               "P2b_nap": sm["selective"]["P2b_nap_prev05"],
               "fit_s": s["fit_s"], "predict_s": s["predict_s"],
               "total_s": float(regrow[n]["total_s"]), "threads": s["threads"],
               "peak_rss_mb": s.get("peak_rss_mb"), "tree_clean": s.get("tree_clean")}
        for k in groupings:
            row[k] = float(np.mean(w[k]))
        if n != HEADLINE:
            d = paired_nb(p1v, fold_vals(head_pf, "variable", "P1"), ratio)
            row.update({"dP1_head": d["mean"], "dP1_head_se": d["se"], "dP1_head_se_nb": d["se_nb"],
                        "dP1_head_pos": d["n_folds_pos"], "dP1_head_p": d["p"],
                        "dP1_head_folds": d["folds"]})
            for k in groupings:
                dw = paired_nb(w[k], within[HEADLINE][k], ratio)
                row[f"d_{k}"] = dw["mean"]
                row[f"d_{k}_se"] = dw["se"]
        rows.append(row)
    tab = {r["name"]: r for r in rows}

    # ---- M5 vs plain average, and the rung representative
    contrasts = []

    def con(a, b):
        if a in runs and b in runs:
            d = paired_nb(fold_vals(runs[a][1], "variable", "P1"), fold_vals(runs[b][1], "variable", "P1"),
                          ratio)
            contrasts.append({"contrast": f"{a} vs {b}", "A": a, "B": b, **d})
            return d
        return None

    d_sa = con("m5_stack", "m5_avg")
    for x in ["m5_stack", "m5_avg"]:
        for b in INPUTS[1:]:
            con(x, b)
    con("m5_ctl_perm", "m5_stack")
    stack_name = "m5_stack_f1" if "m5_stack_f1" in runs else "m5_stack"
    if stack_name == "m5_stack_f1":
        d_sa = con("m5_stack_f1", "m5_avg")
    rep = None
    if d_sa is not None:
        rep = stack_name if (d_sa["mean"] > 0 and d_sa["mean"] > d_sa["se"]) else "m5_avg"
    sel = []
    for n in [x for x in [rep] + LAM if x in tab]:
        r = tab[n]
        win = r["dP1_head"] > DELTA and r["dP1_head"] > r["dP1_head_se"]
        sel.append({"candidate": n, "dP1": r["dP1_head"], "se": r["dP1_head_se"],
                    "se_nb": r["dP1_head_se_nb"], "replaces": bool(win)})
    winners = [s for s in sel if s["replaces"]]
    selected = winners[0]["candidate"] if len(winners) == 1 else (
        HEADLINE if not winners else "MULTIPLE: " + ", ".join(s["candidate"] for s in winners))

    # ---- weights per class and fold
    weights = {}
    for n in [x for x in ["m5_stack", "m5_stack_f1", "m5_ctl_perm"] if x in runs]:
        mr = runs[n][0]["model_report"]
        weights[n] = {f: {"inputs": [i["name"] for i in runs[n][0]["config"]["params"]["inputs"]],
                          "classes": v["weights"], "inner_oof_P1_variable": v.get("inner_oof_P1_variable"),
                          "inner_cached": v.get("inner_cached"), "inner_s": v.get("inner_s"),
                          "inner_reports": v.get("inner_reports")}
                      for f, v in mr.items()}
    ctl = {}
    if "m5_ctl_perm" in weights:
        k = weights["m5_ctl_perm"]["0"]["inputs"].index(PERM_INPUT)
        wmax = max(c["weights"][k] for v in weights["m5_ctl_perm"].values() for c in v["classes"].values())
        ok_w = wmax <= 0.01
        ok_p = ("m5_stack" in tab) and tab["m5_ctl_perm"]["P1_var"] >= tab["m5_stack"]["P1_var"] - 0.002
        ctl = {"max_weight_permuted": wmax, "weight_ok": bool(ok_w), "p1_ok": bool(ok_p),
               "pass": bool(ok_w and ok_p)}

    # ---- average reproduces the registered input predictions (outer refits are exact)
    repro = {}
    if "m5_avg" in preds and all(i in preds for i in INPUTS):
        ix = preds["m5_avg"].index
        avg = sum(preds[i].loc[ix, preds["m5_avg"].columns].to_numpy(np.float64) for i in INPUTS) / 3
        repro["max_abs_diff_avg_vs_mean_of_registered_inputs"] = float(
            np.nanmax(np.abs(preds["m5_avg"].to_numpy(np.float64) - avg)))

    # ---- lambda diagnostic
    lam = {}
    for n in [HEADLINE] + [x for x in LAM if x in runs]:
        mr = runs[n][0]["model_report"]
        lam[n] = {f: {"lam_by_group": v.get("lam_by_group"), "lam_at_grid_edge": v.get("lam_at_grid_edge"),
                      "score_by_group": v.get("score_by_group")} for f, v in mr.items()}
    dec_runs = {n: ids[n] for n in [HEADLINE] + LAM if n in ids}
    dec, edges = sd_deciles(Y, folds, dec_runs, genes)
    dec_tab = {}
    for n, dd in dec.items():
        for d, v in dd.items():
            e = dec_tab.setdefault(d, {})
            e[n] = float(np.mean(v))
            if n != HEADLINE:
                p = paired(v, dec[HEADLINE][d])
                e[f"{n}_d"] = p["mean"]
                e[f"{n}_d_se"] = p["se"]
                e[f"{n}_d_pos"] = p["n_folds_pos"]
    # P1 on the non-variable ("other") eval genes (lambda group 0 drives them)
    other = {}
    for n in dec_runs:
        t = pd.read_parquet(ROOT / "cache" / "runs" / ids[n] / "per_gene_fold.parquet")
        t = t[~t["is_variable"]]
        other[n] = t.groupby("fold")["r"].mean().sort_index().to_numpy().tolist()
    other_d = {n: paired(other[n], other[HEADLINE]) for n in other if n != HEADLINE}

    out = {"runs": rows, "missing_runs": missing, "contrasts": contrasts, "m5_representative": rep,
           "selection_log": sel, "selected_model": selected, "weights": weights, "control": ctl,
           "reproduction": repro, "lambda": lam, "sd_deciles": dec_tab, "sd_decile_edges": edges,
           "P1_nonvariable_folds": other, "P1_nonvariable_vs_head": other_d, "delta": DELTA,
           "nb_ratio": ratio}
    OUT_JSON.write_text(json.dumps(out, indent=1, default=float) + "\n")
    write_md(out, tab)
    print(json.dumps({"selected": selected, "m5_rep": rep, "missing": missing, "control": ctl,
                      "repro": repro}, indent=1, default=float))


def write_md(out, tab):
    L = ["", "Generated by `scripts/g3_m5_report.py` from `results/registry.csv` (latest ok run per name).",
         "Dev lines only, 5 outer dev folds. P1 = mean per-gene Pearson r on dev-wide variable genes",
         "within each held-out fold, then averaged (± fold SD). ΔP1 = paired fold-wise difference vs",
         "m1_E_resB0; SE = SD/√5; NB SE = Nadeau–Bengio corrected SE (disclosure only). Within-group",
         "P1 columns: scripts/g3_kernel_report.py definitions. P2b = normalised AP, selective genes,",
         "held-out prevalence ≤ 0.5. fit/pred s = sum over folds of the run's fit and predict time;",
         "wall = launcher total (6 threads; machine shared with another agent, indicative only).", "",
         "## Results table (5-fold dev CV)", "",
         "| run | P1 var ± SD | ΔP1 vs m1_E_resB0 (SE / NB SE; folds +) | P1 all | P2b nAP | within-library P1 (Δ) | "
         "within-lineage P1 (Δ) | within lin×lib P1 (Δ) | fit+pred s | wall s |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    for r in out["runs"]:
        dh = (f"{r['dP1_head']:+.4f} ({r['dP1_head_se']:.4f} / {r['dP1_head_se_nb']:.4f}; "
              f"{r['dP1_head_pos']}/5)") if "dP1_head" in r else "—"

        def w(k):
            return f"{r[k]:.4f}" + (f" ({r['d_' + k]:+.4f})" if "d_" + k in r else "")
        L.append(f"| {r['name']} | {r['P1_var']:.4f} ± {r['P1_var_sd']:.4f} | {dh} | {r['P1_all']:.4f} | "
                 f"{r['P2b_nap']:.4f} | {w('within_library')} | {w('within_lineage')} | {w('within_lin_lib')} | "
                 f"{r['fit_s'] + r['predict_s']:.0f} | {r['total_s']:.0f} |")
    L += ["", "Run IDs: " + "; ".join(f"{r['name']} `{r['run_id']}`" for r in out["runs"]), ""]
    L += ["## Paired contrasts (variable-gene P1)", "",
          "| contrast | ΔP1 | SE | NB SE | folds + | p (paired t) |", "|---|---|---|---|---|---|"]
    for c in out["contrasts"]:
        L.append(f"| {c['contrast']} | {c['mean']:+.4f} | {c['se']:.4f} | {c['se_nb']:.4f} | "
                 f"{c['n_folds_pos']}/5 | {c['p']:.3g} |")
    L += ["", f"M5 rung representative (pre-registered: stack if ΔP1 vs average > 0 and > 1 SE): "
          f"**{out['m5_representative']}**", ""]
    L += ["## Selection rule vs m1_E_resB0 (DESIGN §4: ΔP1 > 0.005 and > 1 fold SE)", "",
          "| candidate | ΔP1 | SE | NB SE | replaces |", "|---|---|---|---|---|"]
    for s in out["selection_log"]:
        L.append(f"| {s['candidate']} | {s['dP1']:+.4f} | {s['se']:.4f} | {s['se_nb']:.4f} | "
                 f"{'yes' if s['replaces'] else 'no'} |")
    L += ["", f"**Selected model after G3: `{out['selected_model']}`.**", ""]
    if out["weights"]:
        L += ["## Stacking weights per outer fold and training-fold gene class", "",
              "Weights are normalised NNLS weights (sum 1) in the input order given; n = training-fold genes in the class.", ""]
        for n, wf in out["weights"].items():
            inp = next(iter(wf.values()))["inputs"]
            L += [f"**{n}** (inputs: {', '.join(inp)})", "",
                  "| fold | " + " | ".join(next(iter(wf.values()))["classes"]) + " | inner-OOF P1 variable per input |",
                  "|---|" + "---|" * (len(next(iter(wf.values()))["classes"]) + 1)]
            for f, v in wf.items():
                cells = [f"{', '.join(f'{x:.3f}' for x in c['weights'])} (n {c['n_genes']}"
                         + ("" if c["how"] == "nnls" else f", {c['how']}") + ")" for c in v["classes"].values()]
                oo = v.get("inner_oof_P1_variable") or {}
                L.append(f"| {f} | " + " | ".join(cells) + " | "
                         + ", ".join(f"{k} {x:.4f}" for k, x in oo.items()) + " |")
            L.append("")
    if out["control"]:
        c = out["control"]
        L += ["## Permuted-input control (m5_ctl_perm; m2n_E_paralogs predictions row-permuted)", "",
              f"- max weight on the permuted input over classes and folds: {c['max_weight_permuted']:.4f} "
              f"(expected ≤ 0.01): {'ok' if c['weight_ok'] else 'FAIL'}",
              f"- P1 variable ≥ m5_stack − 0.002: {'ok' if c['p1_ok'] else 'FAIL'}",
              f"- **{'pass' if c['pass'] else 'FAIL'}**", ""]
    if out["reproduction"]:
        L += ["## Reproduction check", "",
              f"- m5_avg predictions vs the mean of the registered input runs' predictions (float32 caches): "
              f"max |Δ| = {out['reproduction']['max_abs_diff_avg_vs_mean_of_registered_inputs']:.2e}", ""]
    L += ["## λ-criterion diagnostic (D20c)", "",
          "λ chosen per training-fold gene group (0 = other, 1 = variable) and fold:", "",
          "| run | fold 0 | fold 1 | fold 2 | fold 3 | fold 4 |", "|---|---|---|---|---|---|"]
    for n, lf in out["lambda"].items():
        cells = []
        for f in sorted(lf):
            lg = lf[f]["lam_by_group"] or {}
            cells.append(" / ".join(f"{float(lg[k]):.3g}" for k in sorted(lg)))
        L.append(f"| {n} | " + " | ".join(cells) + " |")
    L += ["", "P1 by training-fold gene-SD decile (all eval genes; decile 0 = lowest SD; mean per-gene r per fold, then "
          "averaged; Δ vs m1_E_resB0 with fold SE and folds +):", ""]
    lam_runs = [n for n in out["lambda"] if n != HEADLINE]
    L += ["| SD decile | SD range (fold 0) | m1_E_resB0 | " + " | ".join(f"{n} (Δ, SE, +)" for n in lam_runs) + " |",
          "|---|---|---|" + "---|" * len(lam_runs)]
    e0 = out["sd_decile_edges"][0]
    for d in sorted(out["sd_deciles"]):
        e = out["sd_deciles"][d]
        cells = [f"{e[n]:.4f} ({e[n + '_d']:+.4f}, {e[n + '_d_se']:.4f}, {e[n + '_d_pos']}/5)" for n in lam_runs if n in e]
        L.append(f"| {d} | {e0[d]:.3f}–{e0[d + 1]:.3f} | {e[HEADLINE]:.4f} | " + " | ".join(cells) + " |")
    L += ["", "P1 on non-variable eval genes (λ group 0): " + "; ".join(
        f"{n} {np.mean(v):.4f}" for n, v in out["P1_nonvariable_folds"].items()) + ". Δ vs m1_E_resB0: "
        + "; ".join(f"{n} {d['mean']:+.4f} (SE {d['se']:.4f}, {d['n_folds_pos']}/5)"
                    for n, d in out["P1_nonvariable_vs_head"].items()), ""]
    if out["missing_runs"]:
        L += [f"Missing (no ok run): {', '.join(out['missing_runs'])}", ""]
    notes = ""
    if MD.exists():
        t = MD.read_text(encoding="utf-8")
        a, b = t.find("<!-- notes -->"), t.find("<!-- /notes -->")
        if a >= 0 and b > a:
            notes = t[a:b + len("<!-- /notes -->")]
    head = ["# G3 close-out — M5 stacking and the λ-criterion diagnostic (D20c)", "",
            notes or "<!-- notes -->\n<!-- /notes -->"]
    MD.write_text("\n".join(head + L) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
