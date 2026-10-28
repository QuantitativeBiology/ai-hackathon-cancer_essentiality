"""Gate evidence (DESIGN §9, §12): write results/gates/G<k>.json with checks run, pass/fail, commit.

    python scripts/gate.py G0     # runs the full test suite (pytest -q tests) + frozen-hash check
    python scripts/gate.py G1     # baseline sanity on the latest ok b0/b1/b2 registry runs

    python scripts/gate.py G2     # M1/M1b/R1/T/M2 + ablation: scripts/g2_report.py output
    python scripts/gate.py G3     # model ladder closed: headline chosen by the rule, finalists' CV runs,
                                  # frozen tiers, dev-only G4 dry run reproduces CV fold 0
    python scripts/gate.py G4     # one-time test evaluation complete + independent reproduction (DESIGN §9)
    python scripts/gate.py G6     # transfer stages A and B complete + independent reproduction

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


def g3():
    """DESIGN §12 G3 "model ladder closed" (D23): evidence that the G4 finalists are ready to freeze."""
    g2j = GATES / "G2.json"
    if not g2j.exists() or not json.loads(g2j.read_text())["pass"]:
        raise SystemExit("G2 has not passed")
    p = subprocess.run([sys.executable, "-m", "pytest", "-q", "tests", "-p", "no:cacheprovider"],
                       cwd=ROOT, capture_output=True, text=True)
    tail = [l for l in p.stdout.strip().splitlines() if l.strip()][-1]
    fin = json.loads((ROOT / "configs" / "finalists.example.json").read_text())
    entries = [fin["headline"], *fin["additional"], fin["T"], fin["H3_noEC"], *fin["baselines"].values()]
    reg = pd.read_csv(ROOT / "results" / "registry.csv", dtype=str).set_index("run_id")
    cv_bad = []
    for e in entries:
        r = reg.loc[e["cv_run_id"]] if e["cv_run_id"] in reg.index else None
        cfg = json.loads((ROOT / e["config"]).read_text())
        h = __import__("hashlib").sha256(json.dumps(cfg, sort_keys=True).encode()).hexdigest()
        if r is None or r["status"] != "ok" or r["leaky"] != "0" or r["config_hash"] != e["config_hash"]                 or h != e["config_hash"]:
            cv_bad.append(e["name"])
    m5 = json.loads((ROOT / "results" / "g3_m5_summary.json").read_text())
    sel = [s for s in m5["selection_log"] if s["candidate"] == fin["headline"]["name"]]
    tiers = pd.read_csv(ROOT / "configs" / "predictable_genes.csv", dtype=str)
    tracked = subprocess.run(["git", "ls-files", "--error-unmatch", "configs/predictable_genes.csv"], cwd=ROOT,
                             capture_output=True).returncode == 0
    dry = json.loads((ROOT / "results" / "dryrun" / "test" / "summary.json").read_text())
    dry_p1 = {t["finalist"]: t.get("P1_variable") for t in dry["table"]}
    cv_p1 = {k: v["cv_fold0"]["variable"]["P1"] for k, v in dry["cv_comparison"].items()}
    repro = {k: (dry_p1.get(k), cv_p1[k]) for k in cv_p1}
    cold = json.loads((ROOT / "results" / "dryrun" / "test" / "m5_cold_refit_check.json").read_text())
    checks = [
        {"check": "pytest -q tests", "result": tail, "pass": p.returncode == 0},
        {"check": "every finalist's CV run is registered ok, not leaky, and its config hashes to config_hash",
         "finalists": [e["name"] for e in entries], "bad": cv_bad, "pass": not cv_bad},
        {"check": "headline chosen by the DESIGN §4 rule (dP1 > delta and > 1 SE vs the incumbent)",
         "headline": fin["headline"]["name"], "selected_model": m5["selected_model"], "selection": sel,
         "delta": m5["delta"], "pass": m5["selected_model"] == fin["headline"]["name"] and bool(sel)
         and sel[0]["replaces"] and sel[0]["dP1"] > m5["delta"] and sel[0]["dP1"] > sel[0]["se"]},
        {"check": "at most 2 additional finalists besides headline, T, H3 no-EC and baselines (DESIGN §4, §5)",
         "additional": [a["name"] for a in fin["additional"]], "pass": len(fin["additional"]) <= 2},
        {"check": "per-gene tiers committed, primary, 5 folds, from the headline's CV run (D22)",
         "tracked": tracked, "source_run_id": sorted(tiers["source_run_id"].unique()),
         "pass": tracked and set(tiers["table_role"]) == {"primary"}
         and set(tiers["folds_used"]) == {"0;1;2;3;4"}
         and set(tiers["source_run_id"]) == {fin["headline"]["cv_run_id"]}},
        {"check": "dev-only G4 dry run (fold 0 stands in) reproduces CV fold-0 P1 (variable) for every finalist",
         "dry_vs_cv": repro, "pass": all(a is not None and abs(a - b) < 1e-5 for a, b in repro.values())},
        {"check": "headline refit from fresh caches reproduces CV fold-0 predictions",
         "max_abs_diff": cold["max_abs_diff_vs_cv_fold0"], "pass": cold["max_abs_diff_vs_cv_fold0"] < 1e-6},
    ]
    return checks


def g4():
    """DESIGN §12 G4 "unblind test once; independent reproduction"."""
    g3j = GATES / "G3.json"
    if not g3j.exists() or not json.loads(g3j.read_text())["pass"]:
        raise SystemExit("G3 has not passed")
    fin = json.loads((ROOT / "configs" / "finalists.json").read_text())
    entries = [fin["headline"], *fin["additional"], fin["T"], fin["H3_noEC"], *fin["baselines"].values()]
    hashes = sorted(e["config_hash"] for e in entries)
    summ = json.loads((ROOT / "results" / "test" / "summary.json").read_text())
    granted = sorted(a["config_hash"] for a in summ["accesses"] if a["granted"])
    status = {t["finalist"]: t["status"] for t in summ["table"]}
    log = pd.read_csv(ROOT / "results" / "test_access_log.csv", dtype=str)
    guarded = log[(log["access_type"] == "targets") & (log["scope"] == "test") & (log["granted"] == "True")]
    comp = pd.read_csv(ROOT / "reviews" / "g4_review" / "comparison.csv")
    maxdiff = float(pd.to_numeric(comp["diff"], errors="coerce").abs().max())
    review = (ROOT / "reviews" / "g4_review" / "REPORT.md").read_text(encoding="utf-8")
    checks = [
        {"check": "real run (not dry), finalists file is configs/finalists.json",
         "dry_run": summ["dry_run"], "finalists_file": summ["finalists_file"],
         "pass": not summ["dry_run"] and summ["finalists_file"].endswith("configs/finalists.json")},
        {"check": "exactly one granted test access per finalist config_hash",
         "n_granted": len(granted), "pass": granted == hashes},
        {"check": "every finalist's results complete", "status": status,
         "pass": all(v == "complete" for v in status.values()) and len(status) == len(entries)},
        {"check": "the access log has exactly one granted test-target row per finalist config_hash",
         "n": int(len(guarded)), "pass": sorted(guarded["config_hash"]) == hashes},
        {"check": "independent reproduction (own metric code): every compared number within 1e-4",
         "n_compared": int(len(comp)), "max_abs_diff": maxdiff, "pass": maxdiff < 1e-4},
        {"check": "independent review verdict: no critical or major finding",
         "pass": "Critical: none. Major: none." in review},
    ]
    return checks


def g6():
    """DESIGN §12 G6 "apply to set A, then set B; verdicts"; independent reproduction (DESIGN §9)."""
    g4j = GATES / "G4.json"
    if not g4j.exists() or not json.loads(g4j.read_text())["pass"]:
        raise SystemExit("G4 has not passed")
    tr = ROOT / "results" / "transfer"
    sa = json.loads((tr / "setA" / "summary.json").read_text())
    sb = json.loads((tr / "setB" / "summary.json").read_text())
    rows = [l.split() for l in (ROOT / "reviews" / "g6_review" / "cmp_boot.txt").read_text().splitlines()]
    diffs = []
    for r in rows:
        try:
            diffs.append(abs(float(r[-1])))
        except (ValueError, IndexError):
            pass
    review = (ROOT / "reviews" / "g6_review" / "REPORT.md").read_text(encoding="utf-8")
    allowed = {"transfers", "fails", "inconclusive"}
    verdicts = {q: v["verdict"] for q, v in sb["overall"].items()}
    checks = [
        {"check": "set A complete, no errors", "status": sa["status"], "pass": sa["status"] == "complete"
         and not sa.get("errors")},
        {"check": "set B complete, no errors, bootstrap plan matches the frozen plan", "status": sb["status"],
         "plan_matches_frozen": sb["plan_matches_frozen"], "pass": sb["status"] == "complete"
         and not sb.get("errors") and all(sb["plan_matches_frozen"].values())},
        {"check": "an overall verdict for each primary question", "verdicts": verdicts,
         "pass": set(verdicts) == {"Q3a", "Q3b"} and set(verdicts.values()) <= allowed},
        {"check": "independent reproduction (own code): every compared value within 1e-4",
         "n_compared": len(diffs), "max_abs_diff": max(diffs), "pass": bool(diffs) and max(diffs) < 1e-4},
        {"check": "independent review: no critical finding", "pass": "### Critical
None." in review},
    ]
    return checks


def main():
    gate = sys.argv[1]
    checks = {"G0": g0, "G1": g1, "G2": g2, "G3": g3, "G4": g4, "G6": g6}[gate]()
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
