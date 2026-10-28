"""Figures for BRIEFING.md (a plain-language research briefing), drawn from committed results files only.

    python scripts/briefing_figures.py      # -> docs/figures/briefing_*.png

No new analysis: every plotted value is read from results/ or the registry, as REPORT.md quotes it.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
OUT = ROOT / "docs" / "figures"

# reference palette (dataviz skill): categorical slots in fixed order, baselines in a neutral gray
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE, AQUA, GRAY, LIGHTGRAY = "#2a78d6", "#eb6834", "#1baf7a", "#8a8984", "#c9c8c3"
DELTA = 0.005

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "font.family": "DejaVu Sans", "font.size": 10, "text.color": INK, "axes.labelcolor": INK,
    "axes.edgecolor": INK2, "xtick.color": INK2, "ytick.color": INK2, "axes.titlesize": 11,
    "axes.titleweight": "bold", "axes.titlelocation": "left", "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "axes.axisbelow": True, "legend.frameon": False,
})


def load(p):
    return json.loads((RES / p).read_text())


def registry():
    return pd.read_csv(RES / "registry.csv").set_index("run_id")


def save(fig, name):
    fig.savefig(OUT / name, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


def dots(ax, labels, vals, lo=None, hi=None, colors=None, fmt="{:.3f}", label_dx=0.004):
    """Horizontal dot plot, first item on top; optional interval whiskers and a value label at each dot."""
    y = list(range(len(labels)))[::-1]
    for i, (yy, v) in enumerate(zip(y, vals)):
        c = colors[i] if colors else BLUE
        if lo is not None:
            ax.plot([lo[i], hi[i]], [yy, yy], color=c, lw=2, solid_capstyle="round")
        ax.plot(v, yy, "o", ms=8, color=c, mec=SURFACE, mew=2, zorder=3)
        right = hi[i] if hi is not None else v
        ax.text(right + label_dx, yy, fmt.format(v).replace("-", "−"), va="center", fontsize=9, color=INK2)
    ax.set_yticks(y, labels)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="y", visible=False)


# --------------------------------------------------------------------------- figure 1: design and timeline
def fig_design():
    fig = plt.figure(figsize=(11, 6.4))
    ax = fig.add_axes([0, 0.36, 1, 0.64])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 50)
    ax.axis("off")

    def box(x, y, w, h, title, body, edge=INK2, fill="#f3f2ee"):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.2",
                                    fc=fill, ec=edge, lw=1.2))
        ax.text(x + w / 2, y + h - 2.2, title, ha="center", va="top", fontsize=10, weight="bold")
        ax.text(x + w / 2, y + h - 6.6, body, ha="center", va="top", fontsize=8.6, color=INK2, linespacing=1.35)

    def arrow(x0, y0, x1, y1):
        ax.annotate("", (x1, y1), (x0, y0), arrowprops=dict(arrowstyle="-|>", color=INK2, lw=1.2))

    ax.text(0.5, 49, "A  What was predicted, from what, and how it was tested", fontsize=11, weight="bold", va="top")
    box(0.5, 22, 18, 22, "DepMap 26Q3",
        "1,277 screened cancer lines\nwith expression\n\ninputs: expression, copy\nnumber, mutations, lineage\n\n"
        "target: gene effect of\n18,434 genes per line")
    box(24, 31, 20, 13, "Locked test set",
        "255 lines (one fifth)\ngrouped by patient\nsealed by a code guard", edge=ORANGE, fill="#fdf0ea")
    box(24, 9, 20, 18, "Development set",
        "1,022 lines\n5-fold cross-validation\ngrouped by patient\nall model choices made here", edge=BLUE,
        fill="#eaf2fc")
    box(50, 9, 21, 18, "Model ladder",
        "baselines → ridge →\nper-gene corrections →\nboosted trees → averages\n\nkeep a model only if it beats\n"
        "the simpler one by ≥ 0.005", edge=BLUE, fill="#eaf2fc")
    box(50, 31, 21, 13, "One-time test",
        "7 frozen finalists\neach scored once\n(255 unseen lines)", edge=ORANGE, fill="#fdf0ea")
    box(77, 18, 22.5, 26, "Transfer test",
        "published screens of\nnon-cancer cells:\nRPE1 (4 studies), hPSC (1)\n\nscreens re-scored and\nrescaled first; decision\n"
        "rule frozen before\nthe model saw them", edge=AQUA, fill="#e8f6f0")
    arrow(19, 37, 23.6, 37)
    arrow(19, 24, 23.6, 24)
    arrow(44.5, 18, 49.6, 18)
    arrow(60.5, 27.4, 60.5, 30.6)
    arrow(44.5, 37, 49.6, 37)
    arrow(71.5, 37, 76.6, 33)
    ax.text(47, 3, "Every number in the briefing is labelled as development cross-validation (dev CV), test set, or "
                   "transfer.", ha="center", fontsize=8.6, color=INK2, style="italic")

    # B: timeline of gates (hours since the first commit, 2026-10-26 12:37, from git history)
    tl = fig.add_axes([0.04, 0.06, 0.93, 0.24])
    events = [(1.7, "literature,\ndesign v0"), (2.25, "design frozen\n(3 review rounds)"),
              (7.7, "transfer rule\nfrozen"), (9.2, "first models;\nmetrics frozen"),
              (35.3, "ladder closed, finalists\nfrozen; test opened once"),
              (39.8, "transfer\nevaluated"), (41.3, "report\ncomplete")]
    tl.set_xlim(-1, 43)
    tl.set_ylim(-1.6, 1.6)
    tl.axhline(0, color=INK2, lw=1.2, zorder=1)
    for i, (h, lab) in enumerate(events):
        up = i % 2 == 0
        c = ORANGE if "test" in lab else (AQUA if "transfer" in lab else BLUE)
        tl.plot(h, 0, "o", ms=8, color=c, mec=SURFACE, mew=2, zorder=3)
        tl.plot([h, h], [0, 0.45 if up else -0.45], color=GRID, lw=1, zorder=0)
        tl.text(h, 0.55 if up else -0.55, lab, ha="center", va="bottom" if up else "top", fontsize=8, color=INK)
    tl.set_yticks([])
    tl.grid(False)
    for s in ("left",):
        tl.spines[s].set_visible(False)
    tl.set_xlabel("hours since project start (wall clock, one workstation)")
    tl.set_title("B  Timeline of the agent-run project (from the git history)", pad=4)
    save(fig, "briefing_fig1_design.png")


# --------------------------------------------------------------------------- figure 2: model ladder (dev CV)
def fig_ladder():
    reg = registry()
    g3n, g3k, g3m = load("g3_nonlinear_summary.json"), load("g3_kernel_summary.json"), load("g3_m5_summary.json")
    runs = [("B0  gene mean", "20261026-220531_b0_831b9415", GRAY),
            ("B1  lineage mean (shrunk)", "20261026-220921_b1_0e2a7a27", GRAY),
            ("B2  expression nearest neighbours", "20261026-221456_b2_289f5170", GRAY),
            ("R1  published-style ridge", "20261026-223037_r1_E_1df71977", LIGHTGRAY),
            ("M1  kernel ridge on expression", "20261026-222117_m1_E_resB0_91e52f9d", ORANGE),
            ("T  M1 + TP53/RB1 status", "20261026-224125_t_E_status_24b46e2a", ORANGE),
            ("M1b  + copy number, mutations (tuned)", "20261027-002223_m1bt_full_48508032", ORANGE),
            ("M2  + per-gene paralog stage", "20261027-044450_m2n_E_paralogs_62ce1575", ORANGE),
            ("M3  + boosted trees on residuals", "20261027-014824_g3_m3n_601e3f5e", ORANGE),
            ("M5  learned stacking of M1, M2, M3", "20261027-122954_m5_stack_bf84c59d", ORANGE),
            ("M5  plain average of M1, M2, M3", "20261027-161041_m5_avg_6dbebcf4", BLUE)]
    p1 = [reg.loc[r, "P1_variable"] for _, r, _ in runs]
    sd = [0 if pd.isna(reg.loc[r, "P1_variable_sd"]) else reg.loc[r, "P1_variable_sd"] for _, r, _ in runs]

    fig, (a, b) = plt.subplots(1, 2, figsize=(13.5, 5), gridspec_kw={"width_ratios": [1.1, 1], "wspace": 0.75})
    dots(a, [n for n, _, _ in runs], p1, [v - s for v, s in zip(p1, sd)], [v + s for v, s in zip(p1, sd)],
         [c for _, _, c in runs])
    a.set_xlim(-0.01, 0.37)
    a.set_xlabel("P1: mean per-gene correlation across held-out lines\n(variable genes; dev CV; bars = ±1 fold SD)")
    a.set_title("A  Accuracy of each rung")

    # B: paired gain over M1 with SE, against the selection margin δ
    sel = {d["candidate"]: d for d in g3k["selection_log"]}
    main = {d["name"]: d for d in g3n["main"]}
    m5 = {d["name"]: d for d in g3m["runs"]}
    gains = [("T  + TP53/RB1 status", load("g2_summary.json")["diagnostic_contrasts"][3]["mean"],
              load("g2_summary.json")["diagnostic_contrasts"][3]["se"], ORANGE),
             ("M1b  + all genomics, tuned", sel["m1bt_full"]["dP1"], sel["m1bt_full"]["se"], ORANGE),
             ("M2  paralogs", sel["m2n_E_paralogs"]["dP1"], sel["m2n_E_paralogs"]["se"], ORANGE),
             ("M2  interaction partners", sel["m2n_E_partners"]["dP1"], sel["m2n_E_partners"]["se"], ORANGE),
             ("M3  boosted trees", main["g3_m3n"]["dP1"], main["g3_m3n"]["dP1_se"], ORANGE),
             ("M5  learned stacking", m5["m5_stack"]["dP1_head"], m5["m5_stack"]["dP1_head_se"], ORANGE),
             ("M5  plain average (selected)", m5["m5_avg"]["dP1_head"], m5["m5_avg"]["dP1_head_se"], BLUE)]
    b.axvspan(DELTA, 0.0125, color="#eaf2fc", zorder=0)
    b.axvline(DELTA, color=BLUE, lw=1)
    b.axvline(0, color=INK2, lw=1)
    b.text(DELTA + 0.0002, len(gains) - 0.45, "selection margin δ = 0.005", color=INK2, fontsize=8.5, va="bottom")
    dots(b, [g[0] for g in gains], [g[1] for g in gains], [g[1] - g[2] for g in gains], [g[1] + g[2] for g in gains],
         [g[3] for g in gains], fmt="{:+.4f}", label_dx=0.0004)
    b.set_xlim(-0.001, 0.0125)
    b.set_ylim(-0.6, len(gains) - 0.1)
    b.set_xlabel("gain in P1 over M1, paired by fold (dev CV; bars = ±1 SE)")
    b.set_title("B  Which additions earned their place")
    save(fig, "briefing_fig2_ladder.png")


# --------------------------------------------------------------------------- figure 3: modality ablation
def fig_ablation():
    g2, g3k = load("g2_summary.json"), load("g3_kernel_summary.json")
    c = {d["contrast"]: d for d in g2["contrasts"]}
    names = {"C": "copy number", "M": "damaging mutations", "H": "hotspot mutations", "F": "fusions",
             "Meta": "lineage, sex, age"}
    add = [(f"+ {names[k]}", c[f"add {k} to E"]["mean"], c[f"add {k} to E"]["se"]) for k in names]
    loo = [(f"− {names.get(k, 'expression')}", c[f"{k} in full (leave-one-out)"]["mean"],
            c[f"{k} in full (leave-one-out)"]["se"]) for k in ("E", "C", "M", "H", "F", "Meta")]
    tuned_lab = {"m1bt_E_C": "E + C", "m1bt_E_C_M": "E + C + M", "m1bt_E_C_M_H": "E + C + M + H",
                 "m1bt_E_C_M_H_F": "E + C + M + H + F", "m1bt_E_C_M_H_F_Meta": "… + metadata",
                 "m1bt_full": "full set (+ composites)"}
    tuned = [(tuned_lab[d["A"]], d["mean"], d["se"]) for d in g3k["m1b_contrasts"]]

    fig, axs = plt.subplots(1, 3, figsize=(13, 3.9), gridspec_kw={"wspace": 0.75})
    for ax, rows, title, col, xl in [
            (axs[0], add, "A  Add one data type to expression\n(default weights)", ORANGE, (-0.03, 0.006)),
            (axs[1], loo, "B  Remove one data type from the full\nset (default weights)", GRAY, (-0.016, 0.03)),
            (axs[2], tuned, "C  Add data types with tuned\nweights", BLUE, (-0.001, 0.0075))]:
        ax.axvline(0, color=INK2, lw=1)
        dots(ax, [r[0] for r in rows], [r[1] for r in rows], [r[1] - r[2] for r in rows],
             [r[1] + r[2] for r in rows], [col] * len(rows), fmt="{:+.4f}", label_dx=(xl[1] - xl[0]) * 0.03)
        ax.set_xlim(*xl)
        ax.set_title(title)
        ax.set_xlabel("change in P1 (dev CV; ±1 SE)")
    axs[1].text(0.02, -0.21, "positive = the full model needs it", transform=axs[1].transAxes, fontsize=8.5,
                color=INK2)
    axs[2].axvline(DELTA, color=BLUE, lw=1)
    axs[2].text(DELTA, -0.21, "δ = 0.005", transform=axs[2].get_xaxis_transform(), fontsize=8.5, color=INK2,
                ha="center")
    save(fig, "briefing_fig3_ablation.png")


# --------------------------------------------------------------------------- figure 4: test set
def fig_test():
    t = load("test/summary.json")["table"]
    order = [(0, "M5 average (headline)", BLUE), (6, "M1 expression ridge", ORANGE),
             (1, "T transfer model", ORANGE), (5, "B2 expression neighbours", GRAY),
             (4, "B1 lineage mean", GRAY), (3, "B0 gene mean", LIGHTGRAY)]
    fig, axs = plt.subplots(1, 3, figsize=(13.5, 4.1), gridspec_kw={"width_ratios": [1.1, 1, 1], "wspace": 0.55})
    a, b, c = axs
    dots(a, [n for _, n, _ in order], [t[i]["P1_variable"] for i, _, _ in order],
         [t[i]["P1_variable_ci"][0] for i, _, _ in order], [t[i]["P1_variable_ci"][1] for i, _, _ in order],
         [col for _, _, col in order])
    a.set_xlim(-0.01, 0.39)
    a.set_xlabel("P1, variable genes (95% CI)")
    a.set_title("A  Correlation across 255 unseen lines")
    dots(b, [n for _, n, _ in order], [t[i]["P2b_nap_prev05"] for i, _, _ in order],
         [t[i]["P2b_nap_prev05_ci"][0] for i, _, _ in order], [t[i]["P2b_nap_prev05_ci"][1] for i, _, _ in order],
         [col for _, _, col in order])
    b.set_yticklabels([])
    b.set_xlim(-0.01, 0.24)
    b.set_xlabel("normalised AUPRC, selective dependencies (95% CI)\n0 = gene mean, 1 = perfect")
    b.set_title("B  Finding which lines depend on a gene")

    tiers = ["high", "medium", "low"]
    series = [(0, "M5 average", BLUE), (6, "M1 ridge", ORANGE), (5, "B2 neighbours", AQUA), (4, "B1 lineage", GRAY)]
    w = 0.19
    for k, (i, lab, col) in enumerate(series):
        xs = [j + (k - 1.5) * (w + 0.02) for j in range(3)]
        c.bar(xs, [t[i][f"P1_{tt}"] for tt in tiers], width=w, color=col, label=lab)
    c.set_xticks(range(3), ["high\n(1,632)", "medium\n(6,198)", "low\n(10,085)"])
    c.set_xlabel("frozen tier (number of genes)")
    c.grid(axis="x", visible=False)
    c.set_ylabel("test-set P1 (all genes in tier)")
    c.set_title("C  Tiers fixed before the test was opened")
    c.legend(fontsize=8.5, loc="upper right")
    save(fig, "briefing_fig4_test.png")


# --------------------------------------------------------------------------- figure 5: mutation-driven dependencies
def fig_mutations():
    pp = load("mutation_dependencies.json")["per_pair"]
    pairs = [("NRAS_hotspot", "NRAS  |  NRAS hotspot"), ("BRAF_hotspot", "BRAF  |  BRAF hotspot"),
             ("KRAS_hotspot", "KRAS  |  KRAS hotspot"), ("PIK3CA_hotspot", "PIK3CA  |  PIK3CA hotspot"),
             ("ERBB2_amp", "ERBB2  |  ERBB2 amplified"), ("MDM2_TP53WT", "MDM2  |  TP53 wild type"),
             ("STAG2_STAG1", "STAG1  |  STAG2 lost"), ("SMARCA4_SMARCA2", "SMARCA2  |  SMARCA4 lost"),
             ("ARID1A_ARID1B", "ARID1B  |  ARID1A lost")]
    series = [("b1", "B1 lineage mean", GRAY), ("m1_E_resB0", "M1 expression only", ORANGE),
              ("m5_avg", "M5 headline (sees mutation calls)", BLUE)]
    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    h = 0.24
    n = len(pairs)
    for k, (key, lab, col) in enumerate(series):
        ys = [(n - 1 - j) + (1 - k) * (h + 0.03) for j in range(n)]
        ax.barh(ys, [pp[p][key]["ap_all"]["mean"] for p, _ in pairs], height=h, color=col, label=lab)
    ax.set_yticks(range(n)[::-1], [lab for _, lab in pairs])
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="y", visible=False)
    ax.set_xlim(0, 1)
    ax.axhline(4.5, color=GRID, lw=1)
    ax.axhline(2.5, color=GRID, lw=1)
    for y, txt in [(8.55, "RAS/RAF/PI3K hotspots"), (4.4, "amplification, p53"), (2.4, "paralog loss")]:
        ax.text(0.99, y, txt, ha="right", va="top", fontsize=8.5, color=INK2, style="italic")
    ax.set_xlabel("AUPRC: ranking dependent lines above non-dependent ones (dev CV; not prevalence-normalised)")
    ax.set_title("Dependency  |  biomarker")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.3), ncol=3, fontsize=8.5)
    save(fig, "briefing_fig5_mutations.png")


# --------------------------------------------------------------------------- figure 6: transfer
def fig_transfer():
    bg = {b: load(f"transfer/setB/{b}.json")["statistics"] for b in ("RPE1", "hPSC")}
    sec = load("transfer/setB/summary.json")["secondary"]
    names = {"brown2019:WT": "RPE1 Brown 2019, TP53 WT", "brown2019:null": "RPE1 Brown 2019, TP53 null",
             "haapaniemi2018:null": "RPE1 Haapaniemi 2018", "hart2015:WT": "RPE1 Hart 2015",
             "olivieri2020:null": "RPE1 Olivieri 2020", "mair2019:WT": "hPSC Mair 2019"}
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.3), gridspec_kw={"width_ratios": [1.15, 1, 1], "wspace": 0.6})
    a = axs[0]
    rows = []
    for b in ("RPE1", "hPSC"):
        s = bg[b]["Q3a_primary"]["T"]
        rows += [(names[k], s["comparator_abs_per_entry"][k], s["model_abs_per_entry"][k])
                 for k in s["model_abs_per_entry"]]
    y = list(range(len(rows)))[::-1]
    for yy, (lab, cmp_, mod) in zip(y, rows):
        a.plot([cmp_, mod], [yy, yy], color=LIGHTGRAY, lw=2)
        a.plot(cmp_, yy, "o", ms=8, color=GRAY, mec=SURFACE, mew=2, zorder=3)
        a.plot(mod, yy, "o", ms=8, color=BLUE, mec=SURFACE, mew=2, zorder=3)
    a.set_yticks(y, [r[0] for r in rows])
    a.tick_params(axis="y", length=0)
    a.grid(axis="y", visible=False)
    a.set_xlim(0.35, 0.85)
    a.plot([], [], "o", color=GRAY, label="DepMap cancer gene mean")
    a.plot([], [], "o", color=BLUE, label="model T")
    a.set_xlabel("AP for the screen's top-10% depleted genes")
    a.set_title("A  Each published screen")

    for ax, q, title, xl in [(axs[1], "Q3a_primary", "B  Core ranking (Q3a): gain over\nthe gene mean", (-0.11, 0.12)),
                             (axs[2], "Q3b_primary", "C  Context-specific dependencies (Q3b):\ngain over lineage mean",
                              (-0.08, 0.18))]:
        labs, v, lo, hi, cols = [], [], [], [], []
        for b in ("RPE1", "hPSC"):
            for who, src, col in [("T", bg[b][q]["T"], BLUE), ("headline M5", sec[b][f"{q}_headline"], ORANGE),
                                  ("B2 neighbours", sec[b][f"{q}_B2"], GRAY)]:
                labs.append(f"{b}: {who}")
                v.append(src["point"]); lo.append(src["ci_lo"]); hi.append(src["ci_hi"]); cols.append(col)
        ax.axvline(0, color=INK2, lw=1)
        dots(ax, labs, v, lo, hi, cols, fmt="{:+.3f}", label_dx=0.005)
        for k, b in enumerate(("RPE1", "hPSC")):
            top, bot = 5 - 3 * k + 0.35, 3 - 3 * k - 0.35
            ax.plot([bg[b][q]["ceiling"]["MDE"]] * 2, [bot, top], color=ORANGE, lw=2.5)
            ax.plot(bg[b][q]["ceiling"]["ceiling_skill"], top - 0.35, "D", ms=7, color=AQUA, mec=SURFACE, mew=1.5,
                    zorder=4)
        ax.set_xlim(*xl)
        ax.set_title(title)
        ax.set_xlabel("skill (difference; 95% CI)")
    hpsc_ceil = bg["hPSC"]["Q3a_primary"]["ceiling"]["ceiling_skill"]
    axs[1].annotate("ceiling below 0:\nratio condition vacuous", (hpsc_ceil, 2.0), xytext=(hpsc_ceil + 0.004, 2.3), va="bottom",
                    fontsize=8, color=INK2, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
    axs[1].plot([], [], color=ORANGE, lw=2.5, label="pre-registered MDE (threshold for \"fails\")")
    axs[1].plot([], [], "D", color=AQUA, label="agreement between independent screens of the same cells (ceiling)")
    handles = a.get_legend_handles_labels()[0] + axs[1].get_legend_handles_labels()[0]
    labels = a.get_legend_handles_labels()[1] + axs[1].get_legend_handles_labels()[1]
    fig.legend(handles, labels, fontsize=8.5, loc="lower center", bbox_to_anchor=(0.5, -0.1), ncol=4)
    save(fig, "briefing_fig6_transfer.png")


# --------------------------------------------------------------------------- figure 7: runtime vs accuracy
def fig_runtime():
    reg = registry()
    pts = [("B1", "20261026-220921_b1_0e2a7a27", GRAY), ("B2", "20261026-221456_b2_289f5170", GRAY),
           ("R1 ridge", "20261026-223037_r1_E_1df71977", LIGHTGRAY), ("M1", "20261026-222117_m1_E_resB0_91e52f9d", ORANGE),
           ("M1b tuned", "20261027-002223_m1bt_full_48508032", ORANGE),
           ("M2 paralogs", "20261027-044450_m2n_E_paralogs_62ce1575", ORANGE),
           ("M3 trees", "20261027-014824_g3_m3n_601e3f5e", ORANGE),
           ("M5 stacking", "20261027-122954_m5_stack_bf84c59d", ORANGE)]
    tot = {n: reg.loc[r, "fit_s"] + reg.loc[r, "predict_s"] for n, r, _ in pts}
    offset = {"M1": (-6, 6, "right"), "R1 ridge": (6, -12, "left"), "M1b tuned": (-4, -14, "center"),
              "M3 trees": (4, -14, "left"), "M2 paralogs": (-4, 7, "center")}
    fig, ax = plt.subplots(figsize=(7.6, 4.4))
    for n, r, col in pts:
        ax.plot(tot[n] / 60, reg.loc[r, "P1_variable"], "o", ms=8, color=col, mec=SURFACE, mew=2, zorder=3)
        dx, dy, ha = offset.get(n, (6, -3, "left"))
        ax.annotate(n, (tot[n] / 60, reg.loc[r, "P1_variable"]), xytext=(dx, dy), textcoords="offset points",
                    fontsize=8.5, color=INK2, ha=ha)
    avg = tot["M1"] + tot["M2 paralogs"] + tot["M3 trees"]
    p_avg = reg.loc["20261027-161041_m5_avg_6dbebcf4", "P1_variable"]
    ax.plot(avg / 60, p_avg, "o", ms=9, color=BLUE, mec=SURFACE, mew=2, zorder=4)
    ax.annotate("M5 average (headline)", (avg / 60, p_avg), xytext=(-8, 8), textcoords="offset points",
                fontsize=8.5, color=INK, ha="right")
    ax.set_xscale("log")
    ax.set_xlabel("compute for 5-fold cross-validation, minutes (log scale; shared machine, indicative)")
    ax.set_ylabel("dev-CV P1 (variable genes)")
    ax.set_ylim(0.17, 0.33)
    save(fig, "briefing_fig7_runtime.png")


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    fig_design()
    fig_ladder()
    fig_ablation()
    fig_test()
    fig_mutations()
    fig_transfer()
    fig_runtime()
