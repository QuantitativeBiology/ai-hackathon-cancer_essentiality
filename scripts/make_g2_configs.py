"""Write the G2 run configs (DESIGN §5 M1/M1b/R1/T/M2, §6 ablation) to configs/runs/.

    python scripts/make_g2_configs.py

Every config uses the launcher's common protocol (5 outer dev folds, inner grouped folds,
reference baselines B0 and B1 computed inline on the same folds). The lambda grid
(logspace(2, 7, 16)), leave-patient-out LOO, NaN-pattern merging (1%) and the <50-row B1 fallback
are model defaults (src/essml/models/kernel_ridge.py).
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "configs" / "runs"
KAPPA = [0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 10000]
FULL = ["E", "C", "M", "H", "F", "Meta", "composites"]


def base(name, model, params, **extra):
    cfg = {"name": name, "model": model, "features": [], "params": params, "grid": {},
           "reference_baselines": {"B0": {}, "B1": {"grid": {"kappa": KAPPA}}},
           "inner_folds": 5, "min_obs": 40, "inner_min_obs": 30, "threads": 12,
           "requires_gate": "G1", "cache_reference": True}
    cfg.update(extra)
    return cfg


def tag(mods):
    return "_".join("comp" if m == "composites" else m for m in mods)


def configs():
    c = [base("r1_E", "R1", {"modalities": ["E"]}),
         base("m1_E", "M1", {"modalities": ["E"]})]
    for m in ["C", "M", "H", "F", "Meta", "composites"]:  # add-one to E
        c.append(base(f"m1b_{tag(['E', m])}", "M1b", {"modalities": ["E", m]}))
    c.append(base("m1_E_Psub", "M1", {"modalities": ["E"]}, line_subset="P"))
    c.append(base("m1b_E_P_Psub", "M1b", {"modalities": ["E", "P"]}, line_subset="P"))
    c.append(base("m1b_full", "M1b", {"modalities": FULL}))
    for m in FULL:  # leave-one-out from full
        c.append(base(f"m1b_full_minus_{tag([m])}", "M1b", {"modalities": [x for x in FULL if x != m]}))
    c.append(base("m1b_Efree", "M1b", {"modalities": ["C", "M", "H", "F", "Meta"]}))
    c.append(base("m1b_full_prov", "M1b", {"modalities": FULL + ["Meta_prov"]}))
    c.append(base("t_E_status", "T", {"modalities": ["E", "status"]}))
    # negative control (DESIGN §8): feature rows permuted across lines -> dP1 vs B1 ~ 0
    c.append(base("m1_E_permuted", "M1", {"modalities": ["E"], "permute_seed": 0}))
    # diagnostic (added after R1/T > M1-E was seen): M1-E with the B0 instead of the B1 residual
    c.append(base("m1_E_resB0", "M1", {"modalities": ["E"], "residual": "B0"}))
    return c


def m2_configs(m1_modalities, base_model="M1"):
    out = []
    for pk in ["none", "paralogs", "partners"]:
        out.append(base(f"m2_{tag(m1_modalities)}_{pk}", "M2",
                        {"m1_params": {"modalities": list(m1_modalities)}, "base": base_model,
                         "priors": pk}))
    return out


def main():
    import sys
    OUT.mkdir(parents=True, exist_ok=True)
    cfgs = configs()
    if len(sys.argv) > 1 and sys.argv[1] == "--m2":
        cfgs = m2_configs(sys.argv[2].split(","))
    for cfg in cfgs:
        (OUT / f"{cfg['name']}.json").write_text(json.dumps(cfg, indent=1) + "\n")
        print(cfg["name"])


if __name__ == "__main__":
    main()
