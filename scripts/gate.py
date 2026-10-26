"""Gate evidence (DESIGN §9, §12): write results/gates/G<k>.json with checks run, pass/fail, commit.

    python scripts/gate.py G0     # runs the full test suite (pytest -q tests) + frozen-hash check
    python scripts/gate.py G1     # baseline sanity on the latest ok b0/b1/b2 registry runs

    python scripts/gate.py G2     # M1/M1b/R1/T/M2 + ablation: scripts/g2_report.py output

G1 requires a passed G0. Checks (DESIGN §12): B0 P1 is exactly 0 in every fold and gene class;
B1 P1 (variable) > 0 with every fold > 0. B2 is reported, not gated.
"""
import datetime as dt
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
GATES = ROOT / "results" / "gates"


def _git():
    h = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    d = subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"], cwd=ROOT,
                       capture_output=True, text=True).stdout.strip()
    return h, int(bool(d))


def g0():
    p = subprocess.run([sys.executable, "-m", "pytest", "-q", "tests", "-p", "no:cacheprovider"],
                       cwd=ROOT, capture_output=True, text=True)
    tail = [l for l in p.stdout.strip().splitlines() if l.strip()][-1]
    frozen = json.loads((ROOT / "configs" / "frozen_hashes.json").read_text())["files"]
    checks = [{"check": "pytest -q tests", "result": tail, "pass": p.returncode == 0},
              {"check": "frozen hashes present", "files": frozen, "pass": len(frozen) == 4},
              {"check": "configs/eval_genes.csv and configs/gene_classes.csv exist",
               "pass": (ROOT / "configs/eval_genes.csv").exists()
               and (ROOT / "configs/gene_classes.csv").exists()}]
    return checks


def _latest_ok(name):
    reg = pd.read_csv(ROOT / "results" / "registry.csv", dtype=str)
    r = reg[(reg["name"] == name) & (reg["status"] == "ok")]
    if r.empty:
        raise SystemExit(f"no ok run named {name} in registry")
    rid = r.iloc[-1]["run_id"]
    return rid, json.loads((ROOT / "results" / "runs" / rid / "summary.json").read_text())


def g1():
    g0j = GATES / "G0.json"
    if not g0j.exists() or not json.loads(g0j.read_text())["pass"]:
        raise SystemExit("G0 has not passed")
    b0id, b0 = _latest_ok("b0")
    b1id, b1 = _latest_ok("b1")
    b2id, b2 = _latest_ok("b2")
    b0_folds = {c: v["P1_folds"] for c, v in b0["summary"].items()}
    b0_zero = all(x == 0.0 for v in b0_folds.values() for x in v if x is not None)
    b1v = b1["summary"]["variable"]
    checks = [{"check": "B0 P1 exactly 0 in every fold and class", "run_id": b0id,
               "P1_folds": b0_folds, "pass": b0_zero},
              {"check": "B1 P1 (variable) > 0, every fold > 0", "run_id": b1id,
               "P1": b1v["P1"], "P1_folds": b1v["P1_folds"],
               "pass": b1v["P1"] > 0 and all(x > 0 for x in b1v["P1_folds"])},
              {"check": "B2 reported (not gated)", "run_id": b2id,
               "P1_variable": b2["summary"]["variable"]["P1"], "pass": True}]
    return checks


def g2():
    """DESIGN §12 G2 "choose features". Checks fixed before the G2 results were seen."""
    g1j = GATES / "G1.json"
    if not g1j.exists() or not json.loads(g1j.read_text())["pass"]:
        raise SystemExit("G1 has not passed")
    p = subprocess.run([sys.executable, "-m", "pytest", "-q", "tests", "-p", "no:cacheprovider"],
                       cwd=ROOT, capture_output=True, text=True)
    tail = [l for l in p.stdout.strip().splitlines() if l.strip()][-1]
    summ = json.loads((ROOT / "results" / "g2_summary.json").read_text())
    reg = pd.read_csv(ROOT / "results" / "registry.csv", dtype=str)
    runs = {r["name"]: r for r in summ["runs"]}
    sp = pd.read_csv(ROOT / "configs" / "splits.csv", index_col=0)
    n_dev = int((sp["role"] == "dev").sum())
    bad = [n for n, r in runs.items() if not r["line_subset"] and int(r["n_lines"]) != n_dev]
    leaky = reg[reg["run_id"].isin([r["run_id"] for r in runs.values()]) & (reg["leaky"] != "0")]
    old_b1 = reg[(reg["name"] == "b1") & (reg["status"] == "ok")].iloc[0]
    old_b2 = reg[(reg["name"] == "b2") & (reg["status"] == "ok")].iloc[0]
    perm = runs.get("m1_E_permuted", {})
    checks = [
        {"check": "pytest -q tests (incl. kernel-ridge known answers, frozen hashes)", "result": tail,
         "pass": p.returncode == 0},
        {"check": "every planned G2 run has an ok registry row", "missing": summ["missing_runs"],
         "pass": not summ["missing_runs"]},
        {"check": "full-line runs cover all dev lines (5 folds)", "bad": bad, "pass": not bad},
        {"check": "no leaky runs", "pass": leaky.empty},
        {"check": "B1 re-run reproduces the G1 P1 variable (classes changed only rare_selective)",
         "old": float(old_b1["P1_variable"]), "new": runs["b1"]["P1_var"],
         "pass": abs(float(old_b1["P1_variable"]) - runs["b1"]["P1_var"]) < 1e-4},
        {"check": "vectorised B2 reproduces the G1 P1 variable",
         "old": float(old_b2["P1_variable"]), "new": runs["b2"]["P1_var"],
         "pass": abs(float(old_b2["P1_variable"]) - runs["b2"]["P1_var"]) < 1e-4},
        {"check": "negative control: permuted features do not beat B1 (dP1 <= 0.01)",
         "dP1_vs_B1": perm.get("dP1_B1"), "pass": perm.get("dP1_B1") is not None
         and perm["dP1_B1"] <= 0.01},
        {"check": "features chosen by the DESIGN §4 rule", "selected_m1": summ["selected_m1"],
         "selected_model": summ["selected_model"], "pass": summ["selected_m1"] not in (None, "b1")},
    ]
    return checks


def main():
    gate = sys.argv[1]
    checks = {"G0": g0, "G1": g1, "G2": g2}[gate]()
    commit, dirty = _git()
    out = {"gate": gate, "pass": all(c["pass"] for c in checks), "checks": checks,
           "commit": commit, "dirty_tree": dirty,
           "timestamp": dt.datetime.now().isoformat(timespec="seconds")}
    GATES.mkdir(parents=True, exist_ok=True)
    (GATES / f"{gate}.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({"gate": gate, "pass": out["pass"],
                      "checks": [(c["check"], c["pass"]) for c in checks]}, indent=1))
    sys.exit(0 if out["pass"] else 1)


if __name__ == "__main__":
    main()
