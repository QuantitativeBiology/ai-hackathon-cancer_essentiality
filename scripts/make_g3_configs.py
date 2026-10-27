"""Write the G3 kernel-track run configs to configs/runs/ (RESEARCH_TREE D14-D18).

    python scripts/make_g3_configs.py <group> [<group> ...]

Groups (each written and committed BEFORE its runs are launched):
    controls   stronger negative controls on the B0 path (D17; review finding 11)
    m1b_EC     M1b block-weight tuning, E+C (D15; DESIGN §5)
    m1b_more   M1b tuning, full set and nested E+C+M+... (only if tuned E+C "helps", see M1B_RULE)
    m2         M2 on the D14 base (B0-residual M1), pass-through + nested local-lambda CV
    secondary  leave-one-lineage-out and EC leave-out splits (DESIGN §3)
    noec       H3 no-EC variant of T (configs/transfer_criteria.json; DESIGN §7)

Every config uses the common launcher protocol of scripts/make_g2_configs.py (reference baselines
B0/B1 inline on the same folds, inner grouped folds, min_obs 40), 6 BLAS threads (G3 shares the
machine with the nonlinear track), and requires gate G2. Each config states its expected outcome
(controls) or its declared contrast (candidates) in ``expected`` / ``declared``; those keys are
part of the config hash, so they are fixed before the run.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from make_g2_configs import KAPPA  # noqa: E402

OUT = ROOT / "configs" / "runs"
HEADLINE = "m1_E_resB0"  # D17 current best (run 20261026-222117)
FULL = ["E", "C", "M", "H", "F", "Meta", "composites"]  # never Meta_prov (D15)
# the 7 germ-cell/EC lines of DESIGN §3 that are DEV lines (NCCIT, NTERA-2 cl.D1, PA-1, SUSA);
# GCT27 is a test line and 1411H / SCH are not in the screened population (not in splits.csv)
EC_DEV = ["ACH-001578", "ACH-002291", "ACH-001374", "ACH-001668"]
M1B_RULE = ("tuned E+C 'helps' (and m1b_more is run) iff its paired fold-wise dP1 vs m1_E_resB0 > 0 "
            "with two-sided paired t p < 0.05 (DESIGN §6 'helps' = CI excludes 0; first contrast, "
            "so Holm = unadjusted). Selection is separate: DESIGN §4 rule dP1 > 0.005 and > 1 SE.")


def base(name, model, params, **extra):
    cfg = {"name": name, "model": model, "features": [], "params": params, "grid": {},
           "reference_baselines": {"B0": {}, "B1": {"grid": {"kappa": KAPPA}}},
           "inner_folds": 5, "min_obs": 40, "inner_min_obs": 30, "threads": 6,
           "requires_gate": "G2", "cache_reference": True}
    cfg.update(extra)
    return cfg


def controls():
    e0 = {"modalities": ["E"], "residual": "B0"}
    return [
        base("ctl_resB0_permuted", "M1", dict(e0, permute_seed=0),
             expected="Global shuffle of expression rows on the B0 path (headline path): the LOO "
                      "lambda goes to the upper grid edge and predictions ~ B0, so P1 variable "
                      "~ 0 (|P1| < 0.01) and dP1 vs B1 ~ -0.184."),
        base("ctl_resB0_perm_lineage", "M1", dict(e0, permute_seed=0, permute_within="lineage"),
             expected="Expression rows deranged within lineage: only lineage information "
                      "survives, so P1 variable ~ B1 (within [B1 - 0.04, B1 + 0.01]), "
                      "within-lineage P1 ~ 0 (|.| < 0.01), and P1 far below m1_E_resB0 (> 0.08)."),
        base("ctl_resB0_perm_linlib", "M1", dict(e0, permute_seed=0, permute_within="lineage_library"),
             expected="Expression rows deranged within lineage x library group: lineage and "
                      "library information survive, so P1 variable within [B1 - 0.04, B1 + 0.03] "
                      "(library-predictable structure is worth <= ~0.03, review §4), within-lineage "
                      "P1 < 0.03, within-library P1 <= B1's within-library P1 + 0.01, and P1 far "
                      "below m1_E_resB0 (> 0.08)."),
    ]


def m1b(mods, name=None):
    tag = "_".join("comp" if m == "composites" else m for m in mods)
    return base(name or f"m1bt_{tag}", "M1b",
                {"modalities": list(mods), "residual": "B0", "tune_weights": True},
                declared=f"M1b tuned {'+'.join(mods)} vs {HEADLINE}: paired dP1 (Holm over the "
                         f"declared M1b contrasts); selection DESIGN §4 (dP1 > 0.005 and > 1 SE). "
                         + M1B_RULE)


def m1b_more():
    out = [m1b(FULL, "m1bt_full")]
    for k in range(3, len(FULL)):  # E+C+M, E+C+M+H, ... (nested, DESIGN §6 order)
        out.append(m1b(FULL[:k]))
    return out


def m2():
    out = []
    for pk in ["none", "paralogs", "partners"]:
        out.append(base(f"m2n_E_{pk}", "M2",
                        {"m1_params": {"modalities": ["E"], "residual": "B0"}, "base": "M1",
                         "priors": pk, "local_genes": "variable", "nested_lambda": True},
                        declared=f"M2 (D14 base m1_E_resB0, pass-through for non-variable genes, "
                                 f"nested local-lambda CV) priors={pk} vs {HEADLINE}: DESIGN §4 rule; "
                                 f"prior contrasts paralogs-none and partners-paralogs Holm."))
    return out


def secondary(names):
    """Leave-one-lineage-out and EC leave-out versions of the given base configs."""
    src = {n: json.loads((OUT / f"{n}.json").read_text()) for n in names}
    out = []
    for n, c in src.items():
        lolo = dict(c, name=f"{n}_lolo", threads=6, requires_gate="G2",
                    fold_scheme="lineage", min_obs=20, lolo_min_lines=20,
                    expected="Secondary split (DESIGN §3): each dev lineage held out in turn. "
                             "min_obs 20 = DESIGN §4 'within-lineage per-gene r (lineages with >= 20 "
                             "held-out lines)'; lineages with < 20 dev lines are never held out "
                             "(they stay in every training set; lolo_min_lines). B1 = B0 for an unseen "
                             "lineage, so its P1 is 0 by construction.")
        out.append(lolo)
        if n == "t_E_status":  # the EC leave-out of T IS the H3 no-EC variant (one run)
            out += noec()
            continue
        ec = dict(c, name=f"{n}_noEC", threads=6, requires_gate="G2",
                  exclude_from_training=EC_DEV,
                  expected="Same 5 outer folds; the 4 dev germ-cell/EC lines are removed from every "
                           "training set (they are still predicted in their own fold). EC lines are "
                           "scored per line (P2a) and on the H3 pluripotency genes.")
        out.append(ec)
    return out


def noec():
    c = json.loads((OUT / "t_E_status.json").read_text())
    return [dict(c, name="t_E_status_noEC", threads=6, requires_gate="G2",
                 exclude_from_training=EC_DEV,
                 expected="H3 no-EC variant (configs/transfer_criteria.json models.T): model T trained "
                          "without the dev germ-cell/EC lines (4 of the DESIGN §3 seven are dev lines). "
                          "Dev CV only here; frozen at G4.")]


GROUPS = {"controls": controls, "m1b_EC": lambda: [m1b(["E", "C"])], "m1b_more": m1b_more,
          "m2": m2, "noec": noec,
          "secondary": lambda: secondary(["b1", "b2", "m1_E_resB0", "t_E_status"])}


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    OUT.mkdir(parents=True, exist_ok=True)
    for g in argv:
        for cfg in GROUPS[g]():
            (OUT / f"{cfg['name']}.json").write_text(json.dumps(cfg, indent=1) + "\n")
            print(cfg["name"])


if __name__ == "__main__":
    main()
