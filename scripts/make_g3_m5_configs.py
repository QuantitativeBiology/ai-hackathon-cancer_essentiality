"""Write the G3 close-out run configs to configs/runs/ (RESEARCH_TREE D20c and M5).

    python scripts/make_g3_m5_configs.py <group> [<group> ...]

Groups (each written and committed BEFORE its runs are launched):
    lamcrit   D20c lambda-criterion diagnostic on the headline m1_E_resB0: LOO-MSE and grouped
              inner-CV r instead of the LOO-r criterion of record
    m5        M5 stacking of m1_E_resB0, m2n_E_paralogs and g3_m3n (DESIGN §5), the plain-average
              baseline and the permuted-input sanity control

Every config follows the launcher protocol of the headline config (configs/runs/m1_E_resB0.json:
reference baselines B0/B1 inline, inner grouped folds, min_obs 40) with 6 BLAS threads.
``declared`` / ``expected`` are part of the config hash, so they are fixed before the run.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from make_g2_configs import KAPPA  # noqa: E402

OUT = ROOT / "configs" / "runs"
HEADLINE = "m1_E_resB0"  # D17 current best (run 20261026-222117)
E0 = {"modalities": ["E"], "residual": "B0"}
RULE = ("DESIGN §4 rule vs m1_E_resB0 (paired fold-wise dP1 on dev-wide variable genes > 0.005 "
        "and > 1 fold SE).")


def base(name, model, params, gate="G1", **extra):
    cfg = {"name": name, "model": model, "features": [], "params": params, "grid": {},
           "reference_baselines": {"B0": {}, "B1": {"grid": {"kappa": KAPPA}}},
           "inner_folds": 5, "min_obs": 40, "inner_min_obs": 30, "threads": 6,
           "requires_gate": gate, "cache_reference": True}
    cfg.update(extra)
    return cfg


def lamcrit():
    d = ("D20c diagnostic: m1_E_resB0 re-fitted with lambda chosen by {what} instead of the LOO "
         "Pearson r of record (lambda still shared within the training-fold gene groups "
         "variable / other; same grid 1e2..1e7, 16 points). " + RULE + " If it clears the rule "
         "it is treated like D14 (dev-CV model selection; test still blind); otherwise it is "
         "reported as robustness only. Also reported: lambda per group and fold, and P1 by "
         "training-fold gene-SD decile (all eval genes) vs m1_E_resB0.")
    return [
        base("m1_E_resB0_loomse", "M1", dict(E0, lam_criterion="loo_mse"),
             declared=d.format(what="the closed-form leave-one(-patient)-out MSE of y, each gene "
                               "scaled by its training variance (mean over the group's genes)")),
        base("m1_E_resB0_icvr", "M1", dict(E0, lam_criterion="inner_cv_r", lam_inner_folds=5),
             declared=d.format(what="grouped inner-CV r: 5 patient-grouped, lineage-stratified "
                               "inner folds of the training lines (seed 20261026); per inner fold "
                               "the pattern ridge is refitted on the inner-training rows, per-gene "
                               "r on the inner-validation lines (>= 20 observed), mean over the "
                               "group's genes, then over inner folds")),
    ]


def _params(name):
    return json.loads((OUT / f"{name}.json").read_text(encoding="utf-8"))["params"]


def m5_inputs(m3_inner=None):
    """The three pre-specified inputs with the params of their committed configs (D19: g3_m3n is
    the only admissible M3 ingredient). M2 inner fits: declared simplification D-M5-1."""
    m3 = {"name": "g3_m3n", "model": "M3", "params": _params("g3_m3n")}
    if m3_inner:
        m3["inner_params"] = m3_inner
    return [{"name": HEADLINE, "model": "M1", "params": _params(HEADLINE)},
            {"name": "m2n_E_paralogs", "model": "M2", "params": _params("m2n_E_paralogs"),
             "inner_params": {"nested_lambda": False}},
            m3]


M5_PREREG = (
    "M5 stacking (DESIGN §5) of m1_E_resB0, m2n_E_paralogs and g3_m3n. Per outer fold: 5 "
    "patient-grouped, lineage-stratified inner folds of the training lines (seed 20261026, the "
    "splits of the cached M1 inner OOF); inner out-of-fold predictions of each input on the "
    "training lines; per training-fold gene class (non-variable; variable x common-essential / "
    "selective / other) non-negative least squares of the training targets on the inner OOF "
    "predictions, per gene centred and scaled by 1/SD(y), weights normalised to sum 1 (classes "
    "with < 20 genes use the pooled variable-gene weights); inputs refitted on all training lines "
    "and combined. Baseline: plain average m5_avg. Rung representative: m5_stack if its paired "
    "dP1 vs m5_avg > 0 and > 1 fold SE, otherwise m5_avg (stacking adds nothing over averaging). "
    "The representative replaces m1_E_resB0 only under the " + RULE + " "
    "DEVIATION D-M5-1 (declared before any M5 run): fully nested inner OOF for M2 needs a third "
    "level of M1 fits (25 per inner split, ~650 fits, projected > 5 h with M3), so the M2 inner "
    "fits choose their local lambda by the non-nested inner-fold CV of G2 (nested_lambda=False) "
    "inside the inner-training lines; this is leak-free with respect to the inner-held-out "
    "lines, but the weights are learned on predictions of this slightly different M2. M1 inner "
    "OOF: exact-key cache cache/m2_inner. M3 inner fits: exact (nested_es=True). Time box 5 h "
    "wall from the launch of m5_stack; if its fold 0 takes > 55 min (projected > 5 h incl. the "
    "average and control runs), it is stopped and replaced by m5_stack_f1 (fallback F1: M3 inner "
    "fits with nested_es=False, i.e. the un-nested early stopping, logged as deviation).")


def m5():
    inp = m5_inputs()
    return [
        base("m5_stack", "M5", {"inputs": inp, "mode": "stack", "inner_folds": 5,
                                "seed": 20261026, "min_class_genes": 20}, gate="G2",
             declared=M5_PREREG),
        base("m5_avg", "M5", {"inputs": inp, "mode": "average"}, gate="G2",
             declared="Plain average of m1_E_resB0, m2n_E_paralogs and g3_m3n (equal weights, "
                      "inputs refitted on all training lines): the baseline M5 must beat. " + RULE),
        base("m5_ctl_perm", "M5", {"inputs": inp, "mode": "stack", "inner_folds": 5,
                                   "seed": 20261026, "min_class_genes": 20,
                                   "permute_input": "m2n_E_paralogs", "permute_seed": 0},
             gate="G2",
             expected="Sanity control: m5_stack with the m2n_E_paralogs predictions row-permuted "
                      "across lines (inner OOF among the training lines, outer predictions among "
                      "the held-out lines). Expected: its normalised weight <= 0.01 in every "
                      "class and fold, and P1 variable not below m5_stack - 0.002 (the other two "
                      "inputs carry the signal)."),
        base("m5_stack_f1", "M5", {"inputs": m5_inputs(m3_inner={"nested_es": False}),
                                   "mode": "stack", "inner_folds": 5, "seed": 20261026,
                                   "min_class_genes": 20}, gate="G2",
             declared="FALLBACK F1 of m5_stack (run only if the pre-declared time-box trigger "
                      "fires): as m5_stack but the M3 inner fits use nested_es=False (un-nested "
                      "early stopping inside the inner-training lines; still leak-free with respect "
                      "to the inner-held-out lines); the statement 'M3 inner fits: exact' in the "
                      "m5_stack text below does NOT apply to F1. " + M5_PREREG),
    ]


GROUPS = {"lamcrit": lamcrit, "m5": m5}


def main(argv):
    for g in argv:
        for cfg in GROUPS[g]():
            p = OUT / f"{cfg['name']}.json"
            p.write_text(json.dumps(cfg, indent=1) + "\n", encoding="utf-8")
            print("wrote", p.relative_to(ROOT))


if __name__ == "__main__":
    main(sys.argv[1:] or list(GROUPS))
