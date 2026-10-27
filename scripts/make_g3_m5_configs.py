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


GROUPS = {"lamcrit": lamcrit}


def main(argv):
    for g in argv:
        for cfg in GROUPS[g]():
            p = OUT / f"{cfg['name']}.json"
            p.write_text(json.dumps(cfg, indent=1) + "\n", encoding="utf-8")
            print("wrote", p.relative_to(ROOT))


if __name__ == "__main__":
    main(sys.argv[1:] or list(GROUPS))
