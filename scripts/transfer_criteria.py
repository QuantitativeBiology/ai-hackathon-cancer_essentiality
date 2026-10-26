"""Write configs/transfer_criteria.json (DESIGN v3 §7.4-7.6 as amended by PI decision D13).

    python scripts/transfer_criteria.py

No model is applied. Inputs: harmonised OBSERVED screens (cache/harmonised, written blind by
scripts/harmonise_screens.py; read here only for the replicate-ceiling computation, logged as
'replicate_ceiling_aggregate'), DepMap DEV targets (B0 and the frozen nearest-lineage B1 comparators;
guard.load_targets on dev IDs), control lists, gene classes, paralog pairs and the expression-profile
table (nearest DepMap lineage per profile). Outputs: the criteria file (every population, threshold,
gene list and lineage mapping, plus the computed ceilings, SEs and MDEs), the frozen comparator vectors
(cache/transfer/comparators.parquet, hashed in the file) and aggregate ceiling tables in
results/harmonisation/ (per entry/study and bootstrap replicates; no gene-level set-B values).
All verdict computations live in src/essml/transfer.py.
"""
from __future__ import annotations

import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from essml import guard  # noqa: E402
from essml import harmonise as H  # noqa: E402
from essml import transfer as TR  # noqa: E402
from essml.baselines import GeneMean, LineageMean  # noqa: E402

CACHE = ROOT / "cache" / "harmonised"
COMP = ROOT / "cache" / "transfer" / "comparators.parquet"
OUT = ROOT / "results" / "harmonisation"
B1_RUN = "results/runs/20261026-151647_b1_0e2a7a27/summary.json"
B1_KAPPA = 10  # modal per-fold kappa of the frozen B1 run (folds 0,1,3,4: 10; fold 2: 20)

# ---------------------------------------------------------------- D13.1 populations (arms per study x TP53 status)
BACKGROUNDS = {
    "RPE1": {
        "expression_input": {"profile_id": "state:DepMap_RPE1", "file": "cache/external_expr/harmonised.parquet",
                             "note": "DepMap 26Q3 RPE1-ss sublines (native DepMap expression, set-A input); identical "
                                     "for every RPE1 entry, so T's prediction differs between entries only via the "
                                     "TP53 flag and the gene universe"},
        "studies": {
            "hart2015": {"role": "verdict", "citation": "Hart 2015 RPE1 TKOv1",
                         "entries": {"WT": {"arms": ["hart2015_RPE1_T18"], "library": "TKOv1"}}},
            "brown2019": {"role": "verdict", "citation": "Brown 2019 (GSE128210) RPE1",
                          "entries": {"WT": {"arms": ["brown2019_WT_Durocher_T3", "brown2019_WT_PACE002_T30",
                                                      "brown2019_WT_PACE006_T25"], "library": "TKOv3",
                                             "note": "one arm per independent screen; PACE002 enters at its final "
                                                     "timepoint T30 (T18 is an earlier timepoint of the same arm)"},
                                      "null": {"arms": ["brown2019_p53null_TKOv2_T18"], "library": "TKOv2"}}},
            "olivieri2020": {"role": "verdict", "citation": "Olivieri 2020 RPE1 TP53-/- NT arms, TKOv3",
                             "entries": {"null": {"arms": ["olivieri2020_TKOv3_S02", "olivieri2020_TKOv3_S04",
                                                           "olivieri2020_TKOv3_S05", "olivieri2020_TKOv3_S06",
                                                           "olivieri2020_TKOv3_S07", "olivieri2020_TKOv3_S12"],
                                                  "library": "TKOv3"}}},
            "haapaniemi2018": {"role": "verdict", "citation": "Haapaniemi 2018 RPE1 Brunello",
                               "entries": {"null": {"arms": ["haapaniemi2018_null_D28"], "library": "Brunello"}}},
        },
        "excluded_arms": {
            "haapaniemi2018_WT_D28": "QC fail (SSMD -0.42)",
            "olivieri2020_TKOv2_S08": "TKOv2 arm: D13 names the Olivieri TKOv3 NT consensus only",
            "olivieri2020_TKOv2_S09": "TKOv2 arm: D13 names the Olivieri TKOv3 NT consensus only",
            "olivieri2020_TKOv2_S10": "TKOv2 arm: D13 names the Olivieri TKOv3 NT consensus only",
            "brown2019_WT_PACE002_T18": "earlier timepoint of the PACE002 arm (T30 used)",
            "bowden2020_WT_D19": "focused DDR library: H2 only",
            "bowden2020_KO_D19": "focused DDR library: H2 only",
            "ORCS_799, ORCS_1791, ORCS_1792, ORCS_1932, ORCS_1933": "rank-only published gene scores (no effect "
                                                                   "scale, no same-pipeline anchor); not in D13 list",
            "ORCS_1934, ORCS_1935": "QC fail", "ORCS_28": "duplicate of hart2015_RPE1_T18",
            "ORCS_1487, ORCS_1488": "117-gene DUB library"},
    },
    "hPSC": {
        "expression_input": {"per_arm": {"mair2019_H1_MEF_T8": "state:H1_iCas9_MEF",
                                         "mair2019_H1_laminin_T12": "state:H1_iCas9_laminin"},
                             "file": "cache/external_expr/harmonised.parquet",
                             "note": "study prediction = mean of the two arm predictions (each with its matched "
                                     "expression state), mirroring the label consensus"},
        "studies": {
            "mair2019": {"role": "verdict", "citation": "Mair 2019 H1 TKOv3 (unsorted MEF T8, laminin T12)",
                         "entries": {"WT": {"arms": ["mair2019_H1_MEF_T8", "mair2019_H1_laminin_T12"],
                                            "library": "TKOv3"}}},
            "ihry2019": {"role": "predictor_only", "citation": "Ihry 2019 hPSC, Novartis library (no anchor)",
                         "entries": {"WT": {"arms": ["ihry2019_hPSC_D18"], "library": "Novartis genome-wide"}}},
            "yilmaz2018": {"role": "predictor_only", "citation": "Yilmaz 2018 h-pES10, Wang library (D12)",
                           "entries": {"WT": {"arms": ["yilmaz2018_hpES10_D23"], "library": "Wang 2015 genome-wide"}}},
        },
        "excluded_arms": {"X_mair2019_H1_MEF_T12_Oct4high": "FACS-sorted arm (DESIGN §7.1)",
                          "ORCS_94, ORCS_2089, ORCS_597": "see screen_plan.csv (top-1000 deposit / differential / duplicate)"},
    },
}
NEAREST_SOURCE = {  # sample-level profiles whose majority nearest lineage is frozen as B1-nearest
    "RPE1": {"states": ["DepMap_RPE1"]},
    "hPSC": {"states": ["H1_iCas9_MEF", "H1_iCas9_laminin"]},  # verdict-eligible study (Mair 2019) only
}
ANCHORS = {  # D13.7 library-matched anchors (library, not lab)
    "TKOv3": [{"id": "A_mair2019_HAP1_T18", "depmap_id": "ACH-002475", "data_level": "counts"}],
    "TKOv1": [{"id": "A_hart2015_DLD1", "depmap_id": "ACH-001061", "data_level": "guide_lfc"},
              {"id": "A_hart2015_HCT116", "depmap_id": "ACH-000971", "data_level": "guide_lfc"}],
    "Brunello": [{"id": "A_sanson2018_A375", "depmap_id": "ACH-000219", "data_level": "counts"}],
    "TKOv2": [],
}
SET_A1 = ["ACH-002462", "ACH-002463", "ACH-002464", "ACH-002465", "ACH-002467"]
H3_GENES = ["POU5F1", "NANOG", "SOX2", "LIN28A", "DNMT3B", "PRDM14", "ZIC2", "L1TD1", "TDGF1", "DPPA4"]


def sha256(path: Path, text=None) -> str:
    b = path.read_bytes()
    if text if text is not None else path.suffix in (".json", ".csv", ".txt", ".py", ".md"):
        b = b.replace(b"\r\n", b"\n")
    return hashlib.sha256(b).hexdigest()


def r4(x):
    return None if x is None or not np.isfinite(x) else round(float(x), 4)


def main():
    seeds = json.loads((ROOT / "configs" / "seeds.json").read_text())
    seed = int(seeds["split_seed"])
    cs = json.loads((ROOT / "configs" / "control_split.json").read_text())
    gc = pd.read_csv(ROOT / "configs" / "gene_classes.csv", dtype={"entrez": str}).set_index("entrez")
    sym2ent = {s: e for e, s in gc["symbol"].items()}
    eval_genes = pd.read_csv(ROOT / "configs" / "eval_genes.csv", dtype=str)["entrez"].tolist()
    plan_csv = pd.read_csv(ROOT / "configs" / "screen_plan.csv").set_index("screen")
    hf = json.loads((ROOT / "configs" / "harmonisation_frozen.json").read_text())
    prof = pd.read_csv(ROOT / "results" / "expression_harmonisation" / "profiles.csv")
    splits = pd.read_csv(ROOT / "configs" / "splits.csv", index_col=0)
    pairs = pd.read_parquet(ROOT / "cache" / "priors" / "paralog_pairs.parquet")
    excl = TR.scaling_halves(cs)
    frac = gc["frac_below"].astype(float)

    # ---- check the D13 population against screen_plan.csv
    for bg, spec in BACKGROUNDS.items():
        for st, s in spec["studies"].items():
            for tp, e in s["entries"].items():
                for a in e["arms"]:
                    r = plan_csv.loc[a]
                    assert r["background"] == bg and r["modality_stratum"] == "KO", a
                    assert bool(r["qc_pass"]) and str(r["treatment"]).startswith("effect-scale"), a
                    assert r["library"] == e["library"], a

    # ---- frozen comparators from DepMap dev (B0; B1 shrunk lineage means, kappa = B1_KAPPA)
    dev = splits.index[splits["role"] == "dev"].tolist()
    Y = guard.load_targets(dev).astype(np.float64)
    lin = splits.loc[dev, "OncotreeLineage"].astype(str)
    b0 = pd.Series(GeneMean().fit(Y).mu_, index=Y.columns.astype(str))
    lm = LineageMean(kappa=B1_KAPPA).fit(Y, lin).lineage_means()
    lm.columns = lm.columns.astype(str)
    del Y

    # ---- D13.4 nearest DepMap lineage per background (majority over sample-level profiles), frozen here
    nearest = {}
    for bg, src in NEAREST_SOURCE.items():
        rows = prof[(prof["level"] == "sample") & prof["state"].isin(src["states"])]
        cnt = Counter(rows["nearest_lineage"])
        top, n = cnt.most_common(1)[0]
        nearest[bg] = {"lineage": top, "votes": f"{n}/{len(rows)}", "profiles": rows["profile_id"].tolist(),
                       "counts": dict(cnt), "n_dev_lines_in_lineage": int((lin == top).sum())}
    comp = pd.DataFrame({"B0": b0})
    for lg in sorted({v["lineage"] for v in nearest.values()} | {"Eye"}):
        comp[f"B1_{lg}"] = lm.loc[lg]
    COMP.parent.mkdir(parents=True, exist_ok=True)
    comp.astype("float64").to_parquet(COMP)

    # ---- observed consensus profiles (replicate-ceiling computation only; aggregates leave this function)
    results, entry_rows, boot_tabs = {}, [], []
    for bg, spec in BACKGROUNDS.items():
        arms = [a for s in spec["studies"].values() for e in s["entries"].values() for a in e["arms"]]
        guard._log("replicate_ceiling_aggregate", pd.Series(["external"] * len(arms), index=arms), True,
                   caller="transfer_criteria.py:main", scope=f"setB:ceiling_d13_{bg}",
                   reason="D13 LOSO replicate-ceiling skill, hierarchical bootstrap (aggregate only; no model)")
        eligible, everyone = TR.build_entries(
            spec, lambda a: pd.read_parquet(CACHE / f"{a}.parquet")["final"].astype(float).rename(a))
        profiles = TR.study_profiles(everyone)
        pools = TR.loso_pools(spec)
        studies = sorted({e.study for e in eligible})
        U_all = sorted(set().union(*[set(TR.gene_universe(e.observed, eval_genes, excl)) for e in eligible]))
        blocks = TR.paralog_blocks(U_all, pairs)
        plan = TR.BootstrapPlan(studies, U_all, blocks, n_boot=TR.N_BOOT, seed=seed)
        b1n = comp[f"B1_{nearest[bg]['lineage']}"]

        # harmonisation validation (control genes only): consensus QC per entry
        cqc = {}
        for e in eligible:
            q = H.qc_metrics(e.observed, cs["ceg_evaluation"], cs["neg_evaluation"])
            cqc[e.name] = {"ssmd": r4(q["ssmd"]), "median_neg_eval": r4(q["median_neg_eval"]),
                           "median_ceg_eval": r4(q["median_ceg_eval"]), "qc_pass": H.qc_pass(q)}
        arm_qc = {a: bool(plan_csv.loc[a, "qc_pass"]) for s in spec["studies"].values() if s["role"] == "verdict"
                  for e in s["entries"].values() for a in e["arms"]}

        stats = {
            "Q3a_primary": dict(q="Q3a", label="rank10", comp="B0", frac_max=TR.FRAC_MAX),
            "Q3a_effect_label": dict(q="Q3a", label="effect", comp="B0", frac_max=TR.FRAC_MAX),
            "Q3b_primary": dict(q="Q3b", label=None, comp="B1n", frac_max=TR.FRAC_MAX),
            "Q3b_vs_B0": dict(q="Q3b", label=None, comp="B0", frac_max=TR.FRAC_MAX),
            "Q3b_rare_vs_B1n": dict(q="Q3b", label=None, comp="B1n", frac_max=TR.FRAC_MAX_RARE),
        }
        res_bg = {}
        for key, d in stats.items():
            evals = {e.name: TR.EntryEval(e, d["q"], eval_genes, excl, frac_below=frac, neg_eval=cs["neg_evaluation"],
                                          label=d["label"] or "rank10", plan_genes=plan.genes, frac_max=d["frac_max"])
                     for e in eligible}
            cvec = b0 if d["comp"] == "B0" else b1n
            res = TR.ceiling_bootstrap(evals, pools, profiles, {n: cvec for n in evals}, b0, plan, keep_boot=True)
            boot_tabs.append(pd.DataFrame({"background": bg, "statistic": key, "replicate": np.arange(plan.n_boot),
                                           "ceiling_skill": res["boot"]}))
            for n, ev in evals.items():
                entry_rows.append({"background": bg, "statistic": key, **ev.summary(),
                                   "ceiling_skill": r4(res["per_entry"].get(n, np.nan)),
                                   "predictor_abs": r4(res["predictor_abs_per_entry"].get(n, np.nan)),
                                   "comparator_abs": r4(res["comparator_abs_per_entry"].get(n, np.nan)),
                                   "predictor_pool": ";".join(pools[ev.entry.study])})
            res_bg[key] = {"ceiling_skill": r4(res["point"]), "bootstrap_se": r4(res["se"]),
                           "ci95": [r4(res["ci_lo"]), r4(res["ci_hi"])], "MDE": r4(res["mde"]),
                           "n_boot_valid": res["n_boot_valid"],
                           "per_study": {s: r4(v) for s, v in res["per_study"].items()},
                           "per_entry": {s: r4(v) for s, v in res["per_entry"].items()},
                           "entries_invalid": res["entries_invalid"],
                           "verdict_bearing": key.endswith("_primary")}
            if key.endswith("_primary"):  # sensitivity: LOSO predictor pools held fixed (not verdict-bearing)
                fx = TR.ceiling_bootstrap(evals, pools, profiles, {n: cvec for n in evals}, b0, plan, resample_pool=False)
                res_bg[key]["sensitivity_fixed_pool"] = {"bootstrap_se": r4(fx["se"]), "MDE": r4(fx["mde"]),
                                                         "ci95": [r4(fx["ci_lo"]), r4(fx["ci_hi"])]}
            print(bg, key, res_bg[key]["ceiling_skill"], res_bg[key]["bootstrap_se"], res_bg[key]["MDE"], flush=True)
        results[bg] = {"ceilings": res_bg, "consensus_qc": cqc, "arm_qc_pass": arm_qc, "plan": plan,
                       "n_blocks": plan.n_blocks, "blocks": blocks,
                       "max_block": int(pd.Series(plan.block_of).value_counts().max())}

    pd.DataFrame(entry_rows).to_csv(OUT / "transfer_ceiling_entries.csv", index=False)
    pd.concat(boot_tabs, ignore_index=True).to_csv(OUT / "transfer_ceiling_bootstrap.csv", index=False)
    write_criteria(results, nearest, comp, seed, hf, sym2ent, gc, splits, eval_genes)


def write_criteria(results, nearest, comp, seed, hf, sym2ent, gc, splits, eval_genes):
    P = hf["parameters"]
    inputs = ["configs/harmonisation_frozen.json", "configs/harmonisation_screens.json", "configs/screen_plan.csv",
              "configs/control_split.json", "configs/gene_classes.csv", "configs/eval_genes.csv",
              "configs/expression_harmonisation.json", "configs/seeds.json", "configs/splits.csv",
              "results/expression_harmonisation/profiles.csv", "cache/priors/paralog_pairs.parquet",
              "cache/external_expr/harmonised.parquet", "cache/transfer/comparators.parquet", B1_RUN,
              "src/essml/transfer.py", "scripts/transfer_criteria.py"]
    ent = lambda syms: {s: sym2ent.get(s) for s in syms}  # noqa: E731
    p53 = hf["p53_ddr_genes"]
    bgs = {}
    for bg, spec in BACKGROUNDS.items():
        R = results[bg]
        anchors = sorted({a["id"] for s in spec["studies"].values() if s["role"] == "verdict"
                          for e in s["entries"].values() for a in ANCHORS.get(e["library"], [])})
        bgs[bg] = {
            "verdict_eligible": True,
            "studies": spec["studies"], "excluded_arms": spec["excluded_arms"],
            "expression_input_T": spec["expression_input"],
            "nearest_depmap_lineage": nearest[bg],
            "B1_nearest_column": f"B1_{nearest[bg]['lineage']}",
            "library_matched_anchors": anchors,
            "expression_swap_profiles": ({"hPSC (Mair 2019 states)": ["state:H1_iCas9_MEF", "state:H1_iCas9_laminin"]}
                                         if bg == "RPE1" else {"RPE1": ["state:DepMap_RPE1"]}),
            "harmonisation_validation": {"arm_qc_pass": R["arm_qc_pass"], "consensus_qc": R["consensus_qc"],
                                         "passed": TR.harmonisation_ok(R["arm_qc_pass"].values(),
                                                                       [v["qc_pass"] for v in R["consensus_qc"].values()])},
            "bootstrap_plan": {"verdict_studies": R["plan"].studies, "n_genes": len(R["plan"].genes),
                               "n_blocks": R["n_blocks"], "max_block_size": R["max_block"]},
            "replicate_ceiling": R["ceilings"],
        }
        if bg == "hPSC":
            bgs[bg]["note"] = ("one verdict-eligible study: the T-skill CI carries gene-block variance only (study "
                               "resampling of one study is the identity); the ceiling CI also carries predictor-pool "
                               "resampling over {Ihry 2019, Yilmaz 2018}, which are predictors only and never evaluated")
    crit = {
        "_status": "PRE-FREEZE (D13). Written by scripts/transfer_criteria.py; binding once committed with tag "
                   "freeze-transfer-criteria. VOID for any set-B prediction unless configs/finalists.json (T, full "
                   "model, H3 no-EC variant, B0/B1/B2 run ids and config hashes) is committed BEFORE that prediction "
                   "(tag freeze-transfer-models). No model has been applied to set A or set B.",
        "decisions": {
            "D7": "skill is computed separately for each of the 5 dev fold models and the 5 skill values are averaged "
                  "(within every bootstrap replicate); the 5-model prediction ensemble is secondary. Anchors: the "
                  "single out-of-fold model (the fold whose training set excludes the anchor line).",
            "D8": "rare-selective >= 1% floor (gene_classes.csv); does not enter the primary statistics",
            "D9": "expression inputs = quantile-normalised to DepMap dev (configs/expression_harmonisation.json)",
            "D10": "harmonised effects on the DepMap-matched affine scale (scaling set -> -1.468, NEG -> 0.007)",
            "D11": "Chronos for TKOv3 not adopted (median guide LFC). The adoption statistic came from a joint fit "
                   "that included set-B count screens; re-enabling Chronos requires a new decision on an "
                   "anchor-only fit.",
            "D12": "no gene-wise correction for Yilmaz 2018 because of a pipeline mismatch: the Wang-library anchors "
                   "are published gene-level scores, Yilmaz is count-scored here. Yilmaz is a predictor-only study.",
            "D13": "this file: study = unit; LOSO replicate ceiling; hierarchical bootstrap; MDE = 2 x SE(ceiling "
                   "skill); transfers needs >= 50% of the ceiling skill; library-matched anchor gate for fails only."},
        "inputs_sha256": {f: sha256(ROOT / f) for f in inputs if (ROOT / f).exists()},
        "models": {"T": "model T as named in configs/finalists.json (E + TP53/RB1 flags, no lineage); "
                        "full = headline model; H3 no-EC = T trained without the 7 germ-cell/EC dev lines",
                   "B0": "per-gene mean over all 1022 dev lines (comparators.parquet column B0)",
                   "B1_nearest": f"shrunk lineage mean (kappa = {B1_KAPPA}, modal per-fold choice of {B1_RUN}) over "
                                 "all 1022 dev lines of the frozen nearest lineage (columns B1_<lineage>)",
                   "B2": "kNN in harmonised expression space, frozen B2 run (secondary comparator only)",
                   "comparator_note": "B0/B1 comparators are the fixed all-dev vectors for every fold model of T "
                                      "(the 80% fold-training versions differ negligibly); for anchors the comparators "
                                      "are refit on the anchor's OOF fold training lines"},
        "harmonisation": {"qc": {"ssmd_max_non_avana": P["qc_ssmd_max"], "nnmd_max_avana": P["qc_nnmd_max"],
                                 "min_eval_controls": P["qc_min_controls"]},
                          "shrinkage_kappa": P["kappa"], "anchor_agreement_improvement_margin": P["improvement_margin"],
                          "affine_targets": [P["affine_target_scaling"], P["affine_target_neg"]],
                          "validation_passed_definition": "every arm of every verdict study has qc_pass in "
                                                          "screen_plan.csv AND every entry consensus passes the same "
                                                          "QC (SSMD of evaluation halves < -1.0, >= 10 of each)"},
        "population": {
            "unit": "study (publication). Entry = study x TP53 status = per-gene NaN-skipping mean of its QC-passing "
                    "genome-wide fitness arms. Study value = mean over its entries; background value = mean over "
                    "verdict studies. Focused libraries and QC-failed arms are excluded from Q3a/Q3b.",
            "gene_universe": "U_e = configs/eval_genes.csv genes observed in the entry consensus minus the scaling "
                             "halves (ceg_scaling U neg_scaling U scaling_set of control_split.json, pinned by the "
                             "affine step); evaluation-half CEGv2/NEGv1 retained",
            "minimum_n": {"universe": TR.MIN_UNIVERSE, "q3b_positives": TR.MIN_Q3B_POS, "q3b_negatives": TR.MIN_Q3B_NEG,
                          "rule": "an entry below a minimum is NaN and dropped; a study with no finite entry is "
                                  "dropped; a background with no finite study is inconclusive; NaN bootstrap "
                                  "replicates are dropped (n_boot_valid reported)"},
            "prediction_sharing": "within RPE1 every entry shares one expression input, so T's prediction is one "
                                  "vector per TP53 flag; the CI describes study choice and gene sampling, not "
                                  "independent replication of the prediction"},
        "Q3a": {"primary": {"label": f"rank label: genes <= the {TR.RANK_Q:.0%} quantile (numpy linear) of the entry "
                                     "consensus over U_e",
                            "statistic": "skill = AP(T) - AP(B0); AP = sklearn step-wise average precision (ties "
                                         "grouped), score = -prediction (essml.transfer.APPrep)",
                            "model": "T only"},
                "secondary_not_verdict_bearing": ["label effect < -0.5", "absolute AUPRC of T, B0 and of the ceiling "
                                                  "predictor", "full model", "B2 comparator", "RPE1 without CEGv2",
                                                  "lost/gained essentials"]},
        "Q3b": {"primary": {"positives": f"consensus < {TR.EFFECT_THR} AND DepMap-dev frac_below <= {TR.FRAC_MAX} "
                                         "(gene_classes.csv), within U_e, excluding NEGv1 evaluation genes",
                            "negatives": "NEGv1 evaluation half (control_split.json neg_evaluation) within U_e",
                            "statistic": f"recall at {TR.FPR:.0%} FPR (threshold = (k+1)-th highest negative score, "
                                         "k = floor(0.05 n_neg), called iff score > threshold; metrics.recall_at_fpr "
                                         "convention) of T minus that of B1-nearest",
                            "model": "T only"},
                "secondary_not_verdict_bearing": ["vs B0", f"frac_below <= {TR.FRAC_MAX_RARE} positives",
                                                  "full model", "B2", "deviation correlation (selective + rare)",
                                                  "Clopper-Pearson recall CIs per entry"]},
        "replicate_ceiling": {"definition": "for entry e of study s: the same statistic with the predictor = mean of "
                                            "the study profiles of the OTHER studies of the background (study profile "
                                            "= mean of its entries; predictor-only studies included), B0 for genes "
                                            "no other study observes; ceiling skill over B0 (Q3a) / B1-nearest (Q3b)",
                              "background_ceiling": "mean over verdict studies of the study mean over entries"},
        "bootstrap": {"n_resamples": TR.N_BOOT, "seed": seed, "seed_source": "configs/seeds.json split_seed",
                      "streams": "np.random.SeedSequence(seed).spawn(3): [0] verdict-study draws (n_studies with "
                                 "replacement), [1] gene-block draws (n_blocks with replacement; gene weight = "
                                 "multiplicity, shared by all entries), [2] ceiling only: each verdict study's LOSO "
                                 "predictor pool resampled with replacement (studies in sorted order)",
                      "pairing": "the T-skill bootstrap and the ceiling bootstrap of a background use identical "
                                 "study and gene draws (same plan genes, blocks and seed)",
                      "blocks": f"paralog families: connected components of cache/priors/paralog_pairs.parquet pairs "
                                f"with min(identity, rev_identity) >= {TR.PARALOG_MIN_IDENTITY:.0f}%; unpaired genes are "
                                "singletons (D13.6; complexes not added)",
                      "plan_genes": "union of U_e over the background's verdict entries",
                      "ci": "percentile 95% (numpy linear)", "se": "SD (ddof=1) of the bootstrap replicates",
                      "MDE": f"{TR.MDE_MULT:g} x SE of the background ceiling skill (no anchor-noise term: struck by D13)"},
        "verdicts": {
            "transfers": f"lower CI(skill) > 0 AND point skill >= {TR.CEILING_FRACTION} x ceiling skill (applied "
                         "literally: with a ceiling skill <= 0 the ratio condition is vacuous and only the lower-CI "
                         "condition binds; flagged in the output)",
            "fails": "upper CI(skill) < MDE (strict) AND harmonisation validation passed AND the library-matched "
                     "anchor gate is informative and positive",
            "both_conditions_met": "inconclusive", "otherwise": "inconclusive",
            "anchor_gate": {"mapping": {k: [a["id"] for a in v] for k, v in ANCHORS.items()},
                            "anchor_details": ANCHORS,
                            "statistic": "the same primary statistic computed on each anchor's harmonised external "
                                         "profile (frozen step: affine; no gene-wise step for TKOv1/TKOv3/Brunello), "
                                         "U = eval genes observed in the anchor minus scaling halves; predictions = the "
                                         "single OOF fold model of T (D7); Q3a comparator = B0 and Q3b comparator = "
                                         "B1 of the anchor's own DepMap lineage, both refit on that fold's training "
                                         "lines (kappa 10); CI = gene-block bootstrap, same seed",
                            "informative_positive": ">= 1 library-matched anchor of the background with lower CI > 0",
                            "TKOv1_usability": "Hart 2015 DLD1/HCT116 are published guide LFCs scored by the same "
                                               "median-guide-LFC step; declared usable (both pass QC); data-level "
                                               "difference recorded",
                            "TKOv2": "no TKOv2 anchor exists (Brown 2019 null entry); the background gate uses the "
                                     "anchors of its other libraries",
                            "role": "used only to attribute 'fails'"},
            "expression_swap_null": "T predicted for this background's entries with each declared swap profile "
                                    "(other verdict background's expression; same TP53 flag as the entry); same "
                                    "statistic, same plan. If any swap has lower CI > 0, 'transfers' is capped at "
                                    "'inconclusive' ('fails' unaffected). The paired matched-minus-swap difference is "
                                    "reported (secondary).",
            "overall": "per question: 'transfers' / 'fails' only if every verdict-eligible background (RPE1, hPSC) has "
                       "that verdict; 'mixed' if both transfers and fails occur; otherwise 'inconclusive'. The "
                       "conjunction is an intersection-union test, so no multiplicity adjustment is applied to the "
                       "2 primary statistics x 2 backgrounds."},
        "backgrounds": {
            **bgs,
            "primary_T": {"verdict_eligible": False, "status": "inconclusive by declaration",
                          "reason": "Shifrut 2018 (SSMD -0.81) and Carnevale 2022 ORCS fail QC; Ting 2018 CD4 is a "
                                    "~2.5k-gene focused chip without anchor or ceiling"},
            "HMEC": {"verdict_eligible": False, "status": "not assessable",
                     "reason": "Tuano 2023 libraries hold < 10 NEGv1 evaluation genes; Sack 2018 score direction unknown"},
            "CRISPRi (iPSC, neurons, microglia)": {"verdict_eligible": False, "status": "descriptive only"},
            "other": {"verdict_eligible": False, "status": "not evaluated (no matched expression profile)",
                      "screens": "fetal NSC (ORCS 572/573/1030/1031), MCF-10A, HUDEP-2, intestinal organoid"}},
        "set_A": {
            "tier_A1": {"lines": SET_A1, "names": "RPE1-ss48, -ss77, -ss6, -ss119, -ss51 (one genetic background)",
                        "lineage": "Eye (DepMap annotation; 6 dev lines): B1 = B1_Eye column (kappa 10); no nearest "
                                   "mapping needed",
                        "statistics": "per line the Q3a/Q3b primary statistics with B1-nearest replaced by B1_Eye; "
                                      "U = eval genes observed in >= 3 of the A1 lines and in the line, minus scaling "
                                      "halves; CI = gene-block bootstrap (same seed); expression = the line's own "
                                      "DepMap expression, TP53 flag WT",
                        "ranking": "within-test percentile = 100 x fraction of cancer-test lines (per-line statistic "
                                   "with B1 of their own lineage, same model) <= the line's value",
                        "descriptive_verdict": "median percentile >= 25: comparable to cancer test lines; 5-25: lower "
                                               "than most; < 5: outside the cancer range (not verdict-bearing)"}},
        "hypotheses": {
            "H1": {"genes": ent(["MDM2", "PPM1D"]),
                   "criterion": "model T (flag WT) - T (flag null), RPE1 expression input, fold-averaged, < 0 for both "
                                "genes AND observed Brown 2019 WT entry - Brown 2019 null entry < 0 for both; sign "
                                "agreement (k/2) reported",
                   "label": "TP53-feature sanity check (single two-point counterfactual); not independent of H2; "
                            "hPSC reported descriptively", "excluded": "USP7 (broadly essential)"},
            "H1_prime": {"type": "diagnostic only", "genes": ent(["TP53", "CDKN1A", "TP53BP1", "USP28", "CHEK2", "ATM"])},
            "H2": {"genes": ent(p53), "gene_source": "DESIGN §7.2 p53/DDR list + p53-target score genes "
                                                     "(harmonisation_frozen.json p53_ddr_genes)",
                   "observed_delta": "per study WT - null consensus delta: Brown 2019 (WT entry - null entry) and "
                                     "Bowden 2020 (WT_D19 - KO_D19, affine-scaled); per gene the mean over the studies "
                                     "observing it (study = unit); each study also reported alone; the same delta on "
                                     "the pre-affine ('raw') scale is reported as a secondary",
                   "criterion": f"n genes >= {TR.H2_MIN_GENES} (else inconclusive) AND Spearman(model delta, observed "
                                f"delta) > {TR.H2_MIN_SPEARMAN} AND sign agreement >= {TR.H2_MIN_SIGN:.0%}, over all "
                                "genes with both deltas"},
            "H3": {"genes": ent(H3_GENES), "not_in_depmap": [g for g in H3_GENES if g not in sym2ent],
                   "criterion": f"among genes observed < {TR.EFFECT_THR} in the Mair 2019 consensus and predicted by T "
                                f"(>= {TR.H3_MIN_OBSERVED} needed, else inconclusive): pass if T predicts < -0.5 for "
                                f">= {TR.H3_PASS_FRACTION:.0%}; the no-EC variant's count is reported; 'narrow EC "
                                "transfer' if T passes and the no-EC variant does not",
                   "prediction": "T for the Mair entry = mean of the two arm predictions (MEF, laminin states)"},
            "H4": {"status": "inconclusive by declaration", "reason": "primary T-cell screens fail QC"}},
        "negative_controls_descriptive": ["oncogene addictions KRAS/BRAF/CTNNB1 predicted ~0", "suppressor tail",
                                          "CN null: SMARCA2, ARID1B, STAG1, ME2/ME3, ENO1 ~0"],
        "superseded": {"results/harmonisation/transfer_ceiling_pairs.csv": "pre-D13 directed-pair ceilings (deleted)",
                       "pre-D13 MDEs": "2 x SE of absolute AUPRC over study representatives (wrong quantity, F1)"},
    }
    (ROOT / "configs" / "transfer_criteria.json").write_text(json.dumps(crit, indent=1, ensure_ascii=False),
                                                             encoding="utf-8")
    print("wrote configs/transfer_criteria.json")


if __name__ == "__main__":
    main()
