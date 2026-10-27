"""Gate G6: apply the frozen models to set A, then set B, and compute the pre-registered verdicts
(DESIGN §7.4-7.6, configs/transfer_criteria.json, D7, D13, D20d).

    # G6 (PI), after G4 passed, on a clean tree with configs/finalists.json committed:
    ESSML_UNBLIND=1 python scripts/evaluate_transfer.py --stage A --threads 4    # set A first (G5a)
    ESSML_UNBLIND=1 python scripts/evaluate_transfer.py --stage B --threads 4    # then set B
    # rehearsal on dev only (stand-in set A / set B built from dev fold-0 lines; no non-dev ID is requested)
    python scripts/evaluate_transfer.py --dry-run-dev --finalists configs/finalists.example.json --threads 4

Every population, threshold, comparator and statistic comes from configs/transfer_criteria.json;
every statistic is computed by the frozen ``essml.transfer`` functions (EntryEval, BootstrapPlan,
skill_bootstrap, ceiling_bootstrap, anchor_gate, background_verdict, overall_verdict,
set_a_descriptive, h1-h4). This script only assembles inputs:

* Models (DESIGN §5 "Frozen models", D7): the 5 dev fold models of each finalist (persisted by G4
  for T, else refit with the CV protocol and checked against the registered CV run's out-of-fold
  predictions). Skill = per fold model, averaged inside every bootstrap replicate; the 5-model
  ensemble is secondary. Anchors (dev lines): the single out-of-fold model; their B0/B1 comparators
  are refit on that model's training lines (kappa 10).
* Predictions: T with the frozen expression input of each entry (``expression_input_T``: one
  profile for RPE1; per-arm states for hPSC, entry prediction = mean of the arm predictions) and the
  entry's TP53 flag (WT -> TP53_mut 0, null -> 1; RB_altered 0). Expression = quantile-normalised
  harmonised profiles (D9; cache/external_expr/harmonised.parquet). Set A uses the lines' own
  DepMap expression and the WT flag.
* H3 no-EC variant: T trained without the 4 EC dev lines (D20d erratum; checked).
* Observed set-B profiles: harmonised 'final' column via ``guard.load_harmonised`` with the
  finalist's config hash, one access per finalist x screen (DESIGN §9); set A via
  ``guard.load_targets(..., screen="A1")``. The tracked access log is committed after each access
  (see scripts/heldout_common.py). Observed data stay in memory.
* Replicate ceiling and MDE: the frozen values in the criteria file (computed pre-freeze with the
  same plan); the model-skill bootstrap uses the same plan (study list, plan genes = union of U_e,
  paralog blocks, seed), so it is paired with the frozen ceiling bootstrap. The plan size is
  checked against the frozen ``bootstrap_plan``.
* Verdicts: ``background_verdict`` per verdict-eligible background x primary question (Q3a, Q3b),
  with the frozen harmonisation validation, the library-matched anchor gate and the expression-swap
  cap; ``overall_verdict`` across backgrounds. Non-eligible backgrounds keep their declared status.
* Secondary (not verdict-bearing; criteria lists): effect label, Q3b vs B0, rare positives, full
  (headline) model, B2, the 5-model ensemble, RPE1 without CEGv2, absolute APs, deviation
  correlation, Clopper-Pearson recall CIs, lost/gained essentials, negative controls; H1-H4, H1'.

Outputs: results/transfer/{setA,setB}/summary.json, results/transfer/summary_<stage>.json,
preflight_<stage>.json, access_ledger_<stage>.json, registry_record_<stage>.csv (per-stage names, so
committing stage A's outputs never makes stage B's tree dirty). Dry run: results/dryrun/transfer/.
"""
from __future__ import annotations

import os
import sys

if __name__ == "__main__":  # thread env before numpy is imported
    _n = int(sys.argv[sys.argv.index("--threads") + 1]) if "--threads" in sys.argv else 4
    for _v in ["OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS",
               "VECLIB_MAXIMUM_THREADS"]:
        os.environ[_v] = str(_n)

import argparse  # noqa: E402
import datetime as dt  # noqa: E402
import json  # noqa: E402
import traceback  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import beta  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
import heldout_common as HC  # noqa: E402
from essml import guard  # noqa: E402
from essml import harmonise as H  # noqa: E402
from essml import metrics as M  # noqa: E402
from essml import transfer as TR  # noqa: E402
from essml.baselines import GeneMean, LineageMean  # noqa: E402

B1_KAPPA = 10
STATS = {  # key -> (question, label, comparator, frac_max); as scripts/transfer_criteria.py
    "Q3a_primary": ("Q3a", "rank10", "B0", TR.FRAC_MAX),
    "Q3a_effect_label": ("Q3a", "effect", "B0", TR.FRAC_MAX),
    "Q3b_primary": ("Q3b", "rank10", "B1n", TR.FRAC_MAX),
    "Q3b_vs_B0": ("Q3b", "rank10", "B0", TR.FRAC_MAX),
    "Q3b_rare_vs_B1n": ("Q3b", "rank10", "B1n", TR.FRAC_MAX_RARE),
}
PRIMARY = {"Q3a": "Q3a_primary", "Q3b": "Q3b_primary"}
SECONDARY_MODELS = ("headline", "B2")
NEG_CONTROL_SYMBOLS = {"oncogene_addiction": ["KRAS", "BRAF", "CTNNB1"],
                       "cn_null": ["SMARCA2", "ARID1B", "STAG1", "ME2", "ME3", "ENO1"]}


def _now():
    return dt.datetime.now().isoformat(timespec="seconds")


def r6(x):
    return None if x is None or not np.isfinite(x) else round(float(x), 6)


# ============================================================================ world (fixed inputs)
class World:
    """Frozen non-target inputs: gene lists, control split, comparators, paralog blocks, symbols."""

    def __init__(self, paths: HC.Paths, crit: dict, comparators: pd.DataFrame):
        self.paths, self.crit = paths, crit
        self.eval_genes = pd.read_csv(paths.eval_genes, dtype={"entrez": str})["entrez"].tolist()
        self.gc = pd.read_csv(paths.gene_classes, dtype={"entrez": str}).set_index("entrez")
        self.frac = self.gc["frac_below"].astype(float)
        self.cs = json.loads(Path(paths.control_split).read_text())
        self.excl = TR.scaling_halves(self.cs)
        self.neg_eval = self.cs["neg_evaluation"]
        self.ceg_all = set(self.cs["ceg_scaling"]) | set(self.cs["ceg_evaluation"])
        self.pairs = (pd.read_parquet(paths.paralog_pairs) if Path(paths.paralog_pairs).exists()
                      else pd.DataFrame(columns=["a", "b", "identity", "rev_identity"]))
        self.comp = comparators
        self.seed = int(crit["bootstrap"]["seed"])
        self.n_boot = int(crit["bootstrap"]["n_resamples"])
        self.sym2ent = {s: e for e, s in self.gc["symbol"].items()} if "symbol" in self.gc else {}


# ============================================================================ predictions
class ModelBank:
    """Fold models per finalist role and cached profile predictions {fold: Series}."""

    def __init__(self, fin, pop, paths, dry, x_builder, expr, feature_loader=None, log=print):
        self.fin, self.pop, self.paths, self.dry = fin, pop, paths, dry
        self.x_builder, self.expr, self.feature_loader, self.log = x_builder, expr, feature_loader, log
        self.models, self.info, self._cache = {}, {}, {}

    def get(self, role):
        if role not in self.models:
            f = self.fin[role]
            models, info = HC.load_or_fit_fold_models(f, self.pop, self.paths, self.dry, self.feature_loader,
                                                      log=self.log)
            self.models[role] = models
            self.info[role] = {k: v for k, v in info.items() if k != "oof_pred"}
            if "oof_pred" in info:
                self.info[role]["cv_reproduction"] = cv_reproduction(f, info["oof_pred"], self.paths)
        return self.models[role]

    def predict(self, role, profile_id, tp53="WT", models=None) -> dict:
        key = (role, profile_id, tp53, None if models is None else tuple(k for k, _, _ in models))
        if key not in self._cache:
            ms = self.get(role) if models is None else models
            out = {}
            for k, m, _ in ms:
                pid, X = self.x_builder(self.fin[role], m, profile_id, self.expr(profile_id), tp53)
                p = m.predict([pid], lineage=None, X=X).iloc[0]
                p.index = p.index.astype(str)
                out[k] = p.astype(np.float64)
            self._cache[key] = out
        return self._cache[key]


def cv_reproduction(f, oof, paths):
    """Refit fold models' OOF predictions vs the registered CV run (same folds and seeds => equal)."""
    if not f.cv_run_id:
        return None
    p = Path(paths.cv_cache) / f.cv_run_id / "pred.parquet"
    if not p.exists():
        return {"cv_run_id": f.cv_run_id, "status": "cached CV predictions not found"}
    P = pd.read_parquet(p)
    P.columns = P.columns.astype(str)
    ids, g = oof.index.intersection(P.index), oof.columns.astype(str).intersection(P.columns)
    o = oof.copy()
    o.columns = o.columns.astype(str)
    return {"cv_run_id": f.cv_run_id, "n_lines": int(len(ids)),
            "max_abs_diff": float(np.abs(o.loc[ids, g].to_numpy(np.float64) - P.loc[ids, g].to_numpy(np.float64)).max())}


def real_x_builder(f, model, profile_id, expr_row: pd.Series, tp53):
    """X for one external/DepMap profile: kernel ridge -> raw blocks {'E', 'status'}; B2 -> 'expression:' frame."""
    pid = f"profile::{profile_id}::{tp53}"
    row = expr_row.astype(np.float64)
    if f.cfg["model"] == "B2":
        return pid, pd.DataFrame([row.add_prefix("expression:").to_numpy()], index=[pid],
                                 columns=[f"expression:{c}" for c in row.index])
    E = pd.DataFrame([row.to_numpy()], index=[pid], columns=row.index.astype(str))
    if tp53 == "depmap":  # anchors (DepMap dev lines): their own DepMap status flags
        from essml.models.kernel_ridge import load_raw_blocks
        st = load_raw_blocks(["status"])["status"].loc[str(profile_id)[5:]]
        status = pd.DataFrame([st.to_numpy(np.float64)], index=[pid], columns=st.index)
    else:
        status = pd.DataFrame({"TP53_mut": [1.0 if tp53 == "null" else 0.0], "RB_altered": [0.0]}, index=[pid])
    return pid, {"E": E, "status": status}


class Expression:
    """Expression profiles: 'line:<ModelID>' = DepMap E row (cache/features); else harmonised external profile."""

    def __init__(self, paths, depmap_E=None):
        self.paths = paths
        self._ext = None
        self._E = depmap_E

    def __call__(self, pid) -> pd.Series:
        if str(pid).startswith("line:"):
            if self._E is None:
                from essml.models.kernel_ridge import load_raw_blocks
                self._E = load_raw_blocks(["E"])["E"]
            return self._E.loc[str(pid)[5:]]
        if self._ext is None:
            self._ext = pd.read_parquet(self.paths.external_expr)
            self._ext.columns = self._ext.columns.astype(str)
        return self._ext.loc[pid]


def mean_folds(preds: list) -> dict:
    """Entry prediction = mean over its arm predictions, per fold (hPSC: per-arm states)."""
    return {k: sum(p[k] for p in preds) / len(preds) for k in preds[0]}


def ens(pred_folds: dict) -> pd.Series:
    return sum(pred_folds.values()) / len(pred_folds)


# ============================================================================ observed data (guarded)
class Observed:
    """One guarded access per (finalist, screen); values cached in memory only."""

    def __init__(self, paths, ledger, dry, standin=None, log=print, fake=None):
        self.paths, self.ledger, self.dry, self.standin, self.log = paths, ledger, dry, standin, log
        self.fake = fake  # rehearsal: (kind, finalist, name) -> DataFrame, no access at all
        self.mem = {}

    def arm(self, f: HC.Finalist, arm, column="final") -> pd.Series:
        key = (f.config_hash, arm)
        if key not in self.mem:
            scope = f"setB:{arm}"
            if self.fake is not None:
                df = self.fake("arm", f, arm)
            elif self.dry:
                df = self.standin(f, arm)
                self.ledger.record(f, scope + "(dry stand-in)", 1, True)
            else:
                self.ledger.check(f.config_hash, scope)
                try:
                    df = guard.load_harmonised(arm, config_hash=f.config_hash, harmonised_dir=self.paths.harmonised_dir,
                                               splits_path=self.paths.splits, log_path=self.paths.access_log,
                                               repo_root=self.paths.root)
                except guard.BlindingError as e:
                    self.ledger.record(f, scope, 1, False, note=str(e)[:300])
                    raise
                c = HC.commit_access_log(self.paths, f"G6 access log: {f.name} ({f.config_hash[:8]}) {scope}")
                self.ledger.record(f, scope, 1, True, commit=c)
            df.index = df.index.astype(str)
            self.mem[key] = df
        df = self.mem[key]
        col = column if column in df.columns else "final"
        return df[col].astype(np.float64).rename(arm)

    def set_a(self, f: HC.Finalist, lines) -> pd.DataFrame:
        key = (f.config_hash, "setA:A1")
        if key not in self.mem and self.fake is not None:
            self.mem[key] = self.fake("setA", f, lines)
        if key not in self.mem:
            Y = HC.guarded_targets(f, lines, self.paths, self.ledger, self.dry, screen="A1",
                                   scope_label="setA:A1(dry stand-in)" if self.dry else None)
            Y.columns = Y.columns.astype(str)
            self.mem[key] = Y.astype(np.float64)
        return self.mem[key]


# ============================================================================ background evaluation
def entries_for(bg, spec, obs: Observed, f, verdict_only=True):
    """TR.build_entries on the verdict studies only (real run: predictor-only studies serve the replicate
    ceiling, which is frozen, so their screens are never accessed); all studies when the ceiling is
    recomputed (dry run)."""
    if verdict_only:
        spec = dict(spec, studies={k: v for k, v in spec["studies"].items() if v["role"] == "verdict"})
    return TR.build_entries(spec, lambda a: obs.arm(f, a))


def arm_profile(spec, arm):
    ei = spec.get("expression_input_T", {})
    return ei.get("per_arm", {}).get(arm) or ei.get("profile_id")


def entry_preds(bank, role, spec, entries, profile_override=None, models=None):
    """{entry name: {fold: Series}}; entry prediction = mean of its arms' profile predictions."""
    out = {}
    for e in entries:
        arms = spec["studies"][e.study]["entries"][e.tp53]["arms"]
        prof = [profile_override or arm_profile(spec, a) for a in arms]
        out[e.name] = mean_folds([bank.predict(role, p, e.tp53, models) for p in prof])
    return out


def make_evals(entries, key, W: World, plan, universe=None):
    q, label, _, fm = STATS[key]
    return {e.name: TR.EntryEval(e, q, W.eval_genes, W.excl, frac_below=W.frac, neg_eval=W.neg_eval, label=label,
                                 plan_genes=plan.genes, frac_max=fm,
                                 universe=None if universe is None else universe(e)) for e in entries}


def comparator(W: World, spec, key, b0=None, b1n=None):
    c = STATS[key][2]
    if c == "B0":
        return W.comp["B0"] if b0 is None else b0
    return W.comp[spec["B1_nearest_column"]] if b1n is None else b1n


def skill_summary(res):
    keep = ("point", "ci_lo", "ci_hi", "se", "n_boot_valid", "n_entries_valid", "entries_invalid")
    out = {k: res.get(k) for k in keep}
    out["per_study"] = {k: r6(v) for k, v in res.get("per_study", {}).items()}
    out["per_entry"] = {k: r6(v) for k, v in res.get("per_entry", {}).items()}
    out["model_abs_per_entry"] = {k: r6(v) for k, v in res.get("model_abs_per_entry", {}).items()}
    out["comparator_abs_per_entry"] = {k: r6(v) for k, v in res.get("comparator_abs_per_entry", {}).items()}
    return out


def clopper_pearson(k, n, level=0.95):
    a = (1 - level) / 2
    lo = 0.0 if k == 0 else float(beta.ppf(a, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(1 - a, k + 1, n - k))
    return lo, hi


def anchor_results(bank, spec, W: World, anchors: dict, key, log=print):
    """D13.7 anchor gate inputs: per library-matched anchor, the statistic with T's single OOF fold model
    (D7), comparators refit on that model's training lines, gene-block bootstrap with the same seed."""
    out = {}
    for aid in spec.get("library_matched_anchors", []):
        if aid not in anchors:
            out[aid] = {"error": "anchor profile unavailable"}
            continue
        prof, depmap_id = anchors[aid]
        models = [(k, m, meta) for k, m, meta in bank.get("T") if depmap_id not in set(meta["train_ids"])]
        if not models:
            out[aid] = {"error": "every fold model trained on the anchor"}
            continue
        tr = sorted(set().union(*[set(meta["train_ids"]) for _, _, meta in models]))
        Y = guard.load_targets(tr, **W.paths.guard_kw()).astype(np.float64)
        lin = bank.pop.lineage.reindex(Y.index)
        b0 = pd.Series(GeneMean().fit(Y).mu_, index=Y.columns.astype(str))
        lm = LineageMean(B1_KAPPA).fit(Y, lin).lineage_means()
        lm.columns = lm.columns.astype(str)
        own = bank.pop.lineage.get(depmap_id)
        b1 = lm.loc[own] if own in lm.index else b0
        del Y
        ent = TR.Entry(name=aid, study=aid, observed=prof)
        U = TR.gene_universe(prof, W.eval_genes, W.excl)
        plan = TR.BootstrapPlan([aid], U, TR.paralog_blocks(U, W.pairs), W.n_boot, W.seed)
        ev = {aid: make_evals([ent], key, W, plan)[aid]}
        pred = bank.predict("T", f"line:{depmap_id}", "depmap", models)
        res = TR.skill_bootstrap(ev, {aid: pred}, {aid: b0 if STATS[key][2] == "B0" else b1}, plan)
        out[aid] = dict(skill_summary(res), depmap_id=depmap_id, oof_folds=[k for k, _, _ in models],
                        n_universe=int(len(U)))
    return out


def evaluate_background(bg, spec, W: World, bank: ModelBank, obs: Observed, fin, anchors, frozen_ceiling=True,
                        secondary=True, log=print):
    out = {"background": bg, "errors": {}}
    T = fin["T"]
    eligible, everyone = entries_for(bg, spec, obs, T, verdict_only=frozen_ceiling)
    studies = sorted({e.study for e in eligible})
    U_all = sorted(set().union(*[set(TR.gene_universe(e.observed, W.eval_genes, W.excl)) for e in eligible]))
    blocks = TR.paralog_blocks(U_all, W.pairs)
    plan = TR.BootstrapPlan(studies, U_all, blocks, W.n_boot, W.seed)
    fp = spec.get("bootstrap_plan", {})
    out["plan"] = {"verdict_studies": plan.studies, "n_genes": len(plan.genes), "n_blocks": plan.n_blocks,
                   "matches_frozen_plan": (fp.get("n_genes") == len(plan.genes) and fp.get("n_blocks") == plan.n_blocks
                                           and fp.get("verdict_studies") == plan.studies) if fp else None}
    # harmonisation validation: frozen value (criteria), recomputed on the consensus for the record
    hv = spec.get("harmonisation_validation", {})
    cqc = {e.name: H.qc_pass(H.qc_metrics(e.observed, W.cs["ceg_evaluation"], W.cs["neg_evaluation"])) for e in eligible}
    harm_ok = bool(hv["passed"]) if "passed" in hv else TR.harmonisation_ok(
        hv.get("arm_qc_pass", {}).values(), cqc.values())
    out["harmonisation"] = {"frozen_passed": hv.get("passed"), "recomputed_consensus_qc": cqc, "used": harm_ok}

    preds_T = entry_preds(bank, "T", spec, eligible)
    stats = {}
    for key in STATS:
        evals = make_evals(eligible, key, W, plan)
        cvec = comparator(W, spec, key)
        res = TR.skill_bootstrap(evals, preds_T, {n: cvec for n in evals}, plan)
        s = {"T": skill_summary(res), "verdict_bearing": key in PRIMARY.values()}
        s["T_ensemble"] = skill_summary(TR.skill_bootstrap(evals, {n: ens(p) for n, p in preds_T.items()},
                                                           {n: cvec for n in evals}, plan))
        rc = spec.get("replicate_ceiling", {}).get(key)
        if frozen_ceiling:  # real run: the frozen ceiling (never recomputed after the transfer freeze)
            s["ceiling"] = ({"ceiling_skill": rc["ceiling_skill"], "MDE": rc.get("MDE"), "source": "frozen criteria"}
                            if rc else {"ceiling_skill": None, "MDE": None, "source": "not in criteria"})
        else:
            pools = TR.loso_pools(spec)
            profiles = TR.study_profiles(everyone)
            cr = TR.ceiling_bootstrap(evals, pools, profiles, {n: cvec for n in evals}, W.comp["B0"], plan)
            s["ceiling"] = {"ceiling_skill": cr["point"], "MDE": cr["mde"], "se": cr["se"],
                            "ci": [cr["ci_lo"], cr["ci_hi"]], "source": "computed (dry run stand-in)"}
        cs_ = s["ceiling"]["ceiling_skill"] if s["ceiling"]["ceiling_skill"] is not None else np.nan
        s["skill_over_ceiling"] = (res["point"] / cs_) if (cs_ and np.isfinite(res["point"]) and cs_ > 0) else None
        if key == "Q3b_primary":  # Clopper-Pearson recall CIs per entry (ensemble prediction; secondary)
            cp = {}
            for n, e in evals.items():
                if e.valid:
                    rec = e.prep(ens(preds_T[n]))()
                    k = int(round(rec * e.n_pos))
                    cp[n] = {"recall": r6(rec), "n_pos": e.n_pos, "ci": clopper_pearson(k, e.n_pos)}
            s["recall_clopper_pearson_T_ensemble"] = cp
        stats[key] = s
        stats[key]["_evals"], stats[key]["_res"] = evals, res
    # expression-swap null (primary statistics; same plan)
    for q, key in PRIMARY.items():
        evals = stats[key]["_evals"]
        cvec = comparator(W, spec, key)
        sw = {}
        for label, profs in spec.get("expression_swap_profiles", {}).items():
            for pid in profs:
                sp = entry_preds(bank, "T", spec, eligible, profile_override=pid)
                r = TR.skill_bootstrap(evals, sp, {n: cvec for n in evals}, plan)
                sw[f"{label}|{pid}"] = skill_summary(r)
        stats[key]["swap"] = sw
        stats[key]["matched_minus_swap_point"] = {k: (stats[key]["T"]["point"] - v["point"])
                                                  if v["point"] is not None and stats[key]["T"]["point"] is not None
                                                  else None for k, v in sw.items()}
    # anchor gate, per primary statistic
    gates = {}
    for q, key in PRIMARY.items():
        try:
            ar = anchor_results(bank, spec, W, anchors, key, log)
            gates[key] = {"anchors": ar, **TR.anchor_gate({a: r for a, r in ar.items() if "error" not in r})}
        except Exception as e:  # noqa: BLE001
            out["errors"][f"anchor_gate_{key}"] = f"{type(e).__name__}: {e}"
            gates[key] = {"informative_positive": False, "error": str(e)}
    out["anchor_gate"] = gates
    # verdicts
    verd = {}
    for q, key in PRIMARY.items():
        s = stats[key]
        cz = s["ceiling"]["ceiling_skill"]
        v = TR.background_verdict(s["_res"], np.nan if cz is None else cz, s["ceiling"]["MDE"] or np.nan, harm_ok,
                                  gates[key].get("informative_positive", False), list(s["swap"].values()))
        verd[q] = v
    out["verdicts"] = verd
    # secondary: RPE1 without CEGv2, deviation correlation, lost/gained essentials, full model, B2
    sec = {}
    try:
        if bg == "RPE1" or spec.get("report_without_cegv2"):
            for key in PRIMARY.values():
                ev2 = make_evals(eligible, key, W, plan,
                                 universe=lambda e: TR.gene_universe(e.observed, W.eval_genes, W.excl | W.ceg_all))
                r = TR.skill_bootstrap(ev2, preds_T, {n: comparator(W, spec, key) for n in ev2}, plan)
                sec[f"{key}_without_CEGv2"] = skill_summary(r)
        sec["deviation_corr_T_ensemble"] = deviation_corrs(eligible, preds_T, W)
        sec["lost_gained_essentials"] = lost_gained(eligible, preds_T, W)
    except Exception as e:  # noqa: BLE001
        out["errors"]["secondary_T"] = f"{type(e).__name__}: {e}"
    if secondary:
        for role in SECONDARY_MODELS:
            try:
                f = fin[role]
                el2, _ = entries_for(bg, spec, obs, f)  # the finalist's own access to each screen
                p2 = entry_preds(bank, role, spec, el2)
                for key in PRIMARY.values():
                    ev2 = make_evals(el2, key, W, plan)
                    r = TR.skill_bootstrap(ev2, p2, {n: comparator(W, spec, key) for n in ev2}, plan)
                    sec[f"{key}_{role}"] = skill_summary(r)
                sec[f"deviation_corr_{role}"] = deviation_corrs(el2, p2, W)
            except Exception as e:  # noqa: BLE001
                out["errors"][f"secondary_{role}"] = f"{type(e).__name__}: {e}"
    out["secondary"] = sec
    out["statistics"] = {k: {kk: vv for kk, vv in v.items() if not kk.startswith("_")} for k, v in stats.items()}
    out["_entries"] = {e.name: e for e in everyone}
    out["_preds_T"] = preds_T
    return out


def deviation_corrs(entries, preds, W: World):
    sel = set(W.gc.index[W.gc["selective"].astype(bool) | W.gc["rare_selective"].astype(bool)])
    b0 = W.comp["B0"]
    out = {}
    for e in entries:
        U = [g for g in TR.gene_universe(e.observed, W.eval_genes, W.excl) if g in sel]
        p = ens(preds[e.name]).reindex(U).to_numpy(np.float64)
        o = e.observed.reindex(U).to_numpy(np.float64)
        out[e.name] = {"n_genes": len(U), "r": r6(M.deviation_corr(p, o, b0.reindex(U).to_numpy(np.float64)))}
    return out


def lost_gained(entries, preds, W: World):
    """Cancer common essentials not essential in the screen (lost) and screen essentials that are not cancer
    common essentials (gained); did T (ensemble) move each gene in the observed direction relative to B0?"""
    ce = set(W.gc.index[W.gc["common_essential"].astype(bool)])
    b0 = W.comp["B0"]
    out = {}
    for e in entries:
        U = TR.gene_universe(e.observed, W.eval_genes, W.excl)
        o = e.observed.reindex(U)
        d = (ens(preds[e.name]).reindex(U) - b0.reindex(U))
        lost = [g for g in U if g in ce and o[g] >= TR.EFFECT_THR]
        gained = [g for g in U if g not in ce and o[g] < TR.EFFECT_THR and W.frac.get(g, np.nan) <= TR.FRAC_MAX]
        out[e.name] = {"n_lost": len(lost), "frac_lost_predicted_up": r6(float(np.mean(d[lost] > 0))) if lost else None,
                       "n_gained": len(gained),
                       "frac_gained_predicted_down": r6(float(np.mean(d[gained] < 0))) if gained else None}
    return out


# ============================================================================ hypotheses
def hypotheses(crit, bg_results, bank, obs, fin, W: World, h2_arms, log=print):
    out = {}
    hyp = crit["hypotheses"]
    rp = bg_results.get("RPE1")
    errs = {}
    # H1 / H2: T(flag WT) - T(flag null) on the RPE1 expression input; observed Brown WT - null (+ Bowden for H2)
    try:
        spec = crit["backgrounds"]["RPE1"]
        pid = spec["expression_input_T"]["profile_id"]
        dWT, dN = ens(bank.predict("T", pid, "WT")), ens(bank.predict("T", pid, "null"))
        md = dWT - dN
        ent = rp["_entries"]
        bw, bn = ent.get(h2_arms["brown_wt_entry"]), ent.get(h2_arms["brown_null_entry"])
        od = (bw.observed - bn.observed) if (bw is not None and bn is not None) else pd.Series(dtype=float)
        g1 = hyp["H1"]["genes"]
        out["H1"] = TR.h1({s: md.get(e, np.nan) for s, e in g1.items()}, {s: od.get(e, np.nan) for s, e in g1.items()},
                          genes=tuple(g1))
        g2 = [e for e in hyp["H2"]["genes"].values() if e]
        deltas = {"brown2019": od.reindex(g2)}
        raw = {}
        if h2_arms.get("bowden_wt") and h2_arms.get("bowden_ko"):
            T = fin["T"]
            bwd = obs.arm(T, h2_arms["bowden_wt"], "affine") - obs.arm(T, h2_arms["bowden_ko"], "affine")
            deltas["bowden2020"] = bwd.reindex(g2)
            raw["bowden2020"] = (obs.arm(T, h2_arms["bowden_wt"], "raw") - obs.arm(T, h2_arms["bowden_ko"], "raw")).reindex(g2)
        out["H2"] = TR.h2(md.reindex(g2), {k: v.dropna() for k, v in deltas.items()})
        out["H2"]["per_study"] = {k: TR.h2(md.reindex(g2), {k: v.dropna()}) for k, v in deltas.items()}
        if raw:
            out["H2"]["bowden2020_raw_scale"] = TR.h2(md.reindex(g2), {k: v.dropna() for k, v in raw.items()})
        # H1': suppressor tail (screen diagnostic), observed per RPE1 entry
        g1p = hyp["H1_prime"]["genes"]
        out["H1_prime"] = {n: {s: r6(e_.observed.get(g, np.nan)) for s, g in g1p.items()} for n, e_ in ent.items()}
    except Exception as e:  # noqa: BLE001
        errs["H1_H2"] = f"{type(e).__name__}: {e}\n{traceback.format_exc()[-800:]}"
    # H3: pluripotency genes in the Mair consensus; T vs the no-EC variant (D20d: 4 EC dev lines excluded)
    try:
        hp = bg_results.get("hPSC")
        spec = crit["backgrounds"]["hPSC"]
        mair = hp["_entries"][h2_arms["mair_entry"]]
        arms = spec["studies"][mair.study]["entries"][mair.tp53]["arms"]
        obs.arm(fin["H3_noEC"], arms[0])  # the no-EC finalist's own access (its prediction is scored on them)
        for a in arms[1:]:
            obs.arm(fin["H3_noEC"], a)
        pT = ens(mean_folds([bank.predict("T", arm_profile(spec, a), mair.tp53) for a in arms]))
        pN = ens(mean_folds([bank.predict("H3_noEC", arm_profile(spec, a), mair.tp53) for a in arms]))
        genes = [e for e in hyp["H3"]["genes"].values() if e]
        out["H3"] = TR.h3(mair.observed, pT, pN, genes)
        out["H3"]["per_gene"] = {g: {"observed": r6(mair.observed.get(g, np.nan)), "T": r6(pT.get(g, np.nan)),
                                     "noEC": r6(pN.get(g, np.nan))} for g in genes}
        out["H3"]["noEC_excluded_lines"] = HC.EC_DEV
    except Exception as e:  # noqa: BLE001
        errs["H3"] = f"{type(e).__name__}: {e}"
    out["H4"] = TR.h4()
    # negative controls (descriptive): predictions and observations near 0
    nc = {}
    for bg, r in bg_results.items():
        for grp, syms in NEG_CONTROL_SYMBOLS.items():
            for s in syms:
                g = W.sym2ent.get(s)
                if not g:
                    continue
                for n, p in r["_preds_T"].items():
                    o = r["_entries"][n].observed.get(g, np.nan)
                    nc.setdefault(bg, {}).setdefault(grp, {}).setdefault(s, {})[n] = {
                        "T": r6(ens(p).get(g, np.nan)), "B0": r6(W.comp["B0"].get(g, np.nan)), "observed": r6(o)}
    out["negative_controls_descriptive"] = nc
    out["errors"] = errs
    return out


# ============================================================================ set A
def evaluate_set_a(crit, W: World, bank: ModelBank, obs: Observed, fin, test_ref: pd.DataFrame | None, b1_col,
                   secondary=True):
    spec = crit["set_A"]["tier_A1"]
    lines = list(spec["lines"])
    out = {"lines": lines, "errors": {}}
    T = fin["T"]
    Y = obs.set_a(T, lines)
    obs_mask = Y.reindex(columns=W.eval_genes).notna()
    in3 = obs_mask.sum(0) >= 3
    res = {}
    for q, key in PRIMARY.items():
        rows = {}
        for ln in lines:
            o = Y.loc[ln].dropna()
            U = pd.Index(sorted(set(TR.gene_universe(o, W.eval_genes, W.excl)) & set(in3.index[in3])))
            plan = TR.BootstrapPlan([ln], U, TR.paralog_blocks(U, W.pairs), W.n_boot, W.seed)
            ent = TR.Entry(name=ln, study=ln, observed=o, tp53="WT")
            ev = {ln: make_evals([ent], key, W, plan, universe=lambda e, U=U: U)[ln]}
            pred = bank.predict("T", f"line:{ln}", "WT")
            cvec = W.comp["B0"] if STATS[key][2] == "B0" else W.comp[b1_col]
            r = TR.skill_bootstrap(ev, {ln: pred}, {ln: cvec}, plan)
            rows[ln] = dict(skill_summary(r), n_universe=int(len(U)))
        dist = None if test_ref is None else test_ref[f"{q}_fold"].to_numpy(np.float64)
        desc = TR.set_a_descriptive({k: {"point": v["point"] if v["point"] is not None else np.nan,
                                         "ci_lo": v["ci_lo"] if v["ci_lo"] is not None else np.nan,
                                         "ci_hi": v["ci_hi"] if v["ci_hi"] is not None else np.nan}
                                     for k, v in rows.items()}, [] if dist is None else dist)
        res[key] = {"per_line": rows, "descriptive": desc,
                    "test_reference_n": None if dist is None else int(np.isfinite(dist).sum())}
    out["statistics"] = res
    # functional p53 check (descriptive diagnostic; DESIGN §3)
    tp53, mdm2 = W.sym2ent.get("TP53", "7157"), W.sym2ent.get("MDM2", "4193")
    out["p53_functional_check"] = {ln: {"TP53": r6(Y.loc[ln].get(tp53, np.nan)), "MDM2": r6(Y.loc[ln].get(mdm2, np.nan)),
                                        "TP53_KO_gt_0.5": bool(Y.loc[ln].get(tp53, np.nan) > 0.5),
                                        "MDM2_lt_-1": bool(Y.loc[ln].get(mdm2, np.nan) < -1)} for ln in lines}
    if secondary:
        try:
            f = fin["headline"]
            Y2 = obs.set_a(f, lines)
            sec = {}
            for q, key in PRIMARY.items():
                rows = {}
                for ln in lines:
                    o = Y2.loc[ln].dropna()
                    U = pd.Index(sorted(set(TR.gene_universe(o, W.eval_genes, W.excl)) & set(in3.index[in3])))
                    plan = TR.BootstrapPlan([ln], U, TR.paralog_blocks(U, W.pairs), W.n_boot, W.seed)
                    ent = TR.Entry(name=ln, study=ln, observed=o, tp53="WT")
                    ev = {ln: make_evals([ent], key, W, plan, universe=lambda e, U=U: U)[ln]}
                    cvec = W.comp["B0"] if STATS[key][2] == "B0" else W.comp[b1_col]
                    rows[ln] = skill_summary(TR.skill_bootstrap(ev, {ln: bank.predict("headline", f"line:{ln}", "WT")},
                                                                {ln: cvec}, plan))
                sec[key] = rows
            out["secondary_headline"] = sec
        except Exception as e:  # noqa: BLE001
            out["errors"]["secondary_headline"] = f"{type(e).__name__}: {e}"
    return out


# ============================================================================ anchors (dev-line external profiles)
def real_anchor_profiles(crit, paths) -> dict:
    """{anchor id: (affine-harmonised external profile, depmap id)} recomputed with the frozen harmonisation
    parameters (scripts/harmonise_screens.py scoring; anchors are DepMap dev lines, not set-B data)."""
    import harmonise_screens as HS
    fr = json.loads((paths.root / "configs" / "harmonisation_frozen.json").read_text())
    H.DEFAULTS.update(fr["parameters"])
    cs = json.loads(Path(paths.control_split).read_text())
    entries = {e["id"]: e for e in json.loads((paths.root / "configs" / "harmonisation_screens.json")
                                              .read_text(encoding="utf-8"))["entries"]}
    det = crit["verdicts"]["anchor_gate"]["anchor_details"]
    want = {a["id"]: a["depmap_id"] for lib in det.values() for a in lib}
    roles = guard.roles_of(list(want.values()), paths.splits)
    if (roles != "dev").any():
        raise guard.BlindingError(f"anchor lines must be dev lines: {roles[roles != 'dev'].to_dict()}")
    sym = HS.SymbolMap()
    out = {}
    for aid, mid in want.items():
        raw, _ = HS.raw_arm_scores(entries[aid], sym)
        ms, mn = H.affine_params(raw, cs["scaling_set"], cs["neg_scaling"])
        out[aid] = (H.affine_apply(raw, ms, mn).astype(np.float64), mid)
    return out


# ============================================================================ dry-run stand-in criteria
def dry_spec(crit, paths, pop, W_comp_lineages):
    """Stand-in backgrounds / set A from dev fold-0 lines (their DepMap CRISPR profiles play the observed
    screens; their DepMap expression the expression inputs); same structure as the criteria file."""
    sp = pd.read_csv(paths.splits, index_col=0)
    sp.index = sp.index.astype(str)
    ho = sp.loc[pop.heldout_ids]
    ho = ho[~ho.index.isin(HC.EC_DEV)]
    comp_status = pd.read_parquet(paths.root / "cache" / "features" / "composite.parquet")["TP53_mut"]
    by = ho.groupby("OncotreeLineage").groups
    lins = sorted(by, key=lambda l: -len(by[l]))
    la, lb, lc = lins[0], lins[1], lins[2]
    A = sorted(by[la])
    wt = [i for i in A if comp_status.get(i, 1) == 0]
    mu = [i for i in A if comp_status.get(i, 0) == 1]
    if len(wt) < 4 or len(mu) < 3:
        raise RuntimeError("dry run needs >= 4 TP53-WT and >= 3 TP53-mutant fold-0 lines in the largest lineage")
    Bl = sorted(by[lb])
    C = sorted(by[lc])[:5]
    lp = lambda ids: {i: f"line:{i}" for i in ids}  # noqa: E731
    c = json.loads(json.dumps(crit))
    rpe = c["backgrounds"]["RPE1"]
    rpe.update({
        "studies": {"s_hart": {"role": "verdict", "entries": {"WT": {"arms": [wt[0]], "library": "TKOv1"}}},
                    "s_brown": {"role": "verdict", "entries": {"WT": {"arms": wt[1:3], "library": "TKOv3"},
                                                               "null": {"arms": [mu[0]], "library": "TKOv2"}}},
                    "s_olivieri": {"role": "verdict", "entries": {"null": {"arms": mu[1:3], "library": "TKOv3"}}},
                    "s_haap": {"role": "verdict", "entries": {"null": {"arms": [wt[3]], "library": "Brunello"}}}},
        "expression_input_T": {"per_arm": lp(wt[:4] + mu[:3])},
        "B1_nearest_column": f"B1_{la}", "replicate_ceiling": {}, "bootstrap_plan": {},
        "expression_swap_profiles": {"hPSC stand-in": [f"line:{Bl[0]}"]},
        "harmonisation_validation": {"passed": True}})
    hp = c["backgrounds"]["hPSC"]
    hp.update({
        "studies": {"s_mair": {"role": "verdict", "entries": {"WT": {"arms": Bl[0:2], "library": "TKOv3"}}},
                    "s_ihry": {"role": "predictor_only", "entries": {"WT": {"arms": [Bl[2]], "library": "x"}}},
                    "s_yilmaz": {"role": "predictor_only", "entries": {"WT": {"arms": [Bl[3]], "library": "x"}}}},
        "expression_input_T": {"per_arm": lp(Bl[:4])},
        "B1_nearest_column": f"B1_{lb}", "replicate_ceiling": {}, "bootstrap_plan": {},
        "expression_swap_profiles": {"RPE1 stand-in": [f"line:{wt[0]}"]},
        "harmonisation_validation": {"passed": True}})
    # H1 needs one RPE1 expression input for the WT-vs-null counterfactual
    rpe["expression_input_T"]["profile_id"] = f"line:{wt[0]}"
    c["set_A"]["tier_A1"]["lines"] = C
    c["set_A"]["tier_A1"]["lineage"] = lc
    h2 = {"brown_wt_entry": "s_brown:WT", "brown_null_entry": "s_brown:null", "bowden_wt": wt[4] if len(wt) > 4 else None,
          "bowden_ko": mu[3] if len(mu) > 3 else None, "mair_entry": "s_mair:WT"}
    used = set(wt[:5] + mu[:4] + Bl[:4] + C)
    return c, h2, {"lineages": {"RPE1": la, "hPSC": lb, "setA": lc}, "lines": sorted(used)}, f"B1_{lc}"


def real_h2_arms():
    return {"brown_wt_entry": "brown2019:WT", "brown_null_entry": "brown2019:null", "bowden_wt": "bowden2020_WT_D19",
            "bowden_ko": "bowden2020_KO_D19", "mair_entry": "mair2019:WT"}


def dry_comparators(paths, pop, lineages):
    Y = guard.load_targets(list(pop.train_ids), **paths.guard_kw()).astype(np.float64)
    b0 = pd.Series(GeneMean().fit(Y).mu_, index=Y.columns.astype(str))
    lm = LineageMean(B1_KAPPA).fit(Y, pop.lineage.reindex(Y.index)).lineage_means()
    lm.columns = lm.columns.astype(str)
    comp = pd.DataFrame({"B0": b0})
    for l in lineages:
        comp[f"B1_{l}"] = lm.loc[l]
    return comp


# ============================================================================ preconditions
def preflight(fin, crit, paths, dry, stage, finalists_rel, secondary, h2_arms, resume=False):
    checks = {}
    allow = [paths.rel(paths.out_dir) + "/"] if resume else []
    dirty = HC.strict_dirty(paths.root, allow)
    checks["tree_clean"] = {"pass": not dirty, "dirty": dirty[:20]}
    checks["finalists_committed"] = {"pass": HC.is_tracked(paths.root, finalists_rel)}
    checks["transfer_freeze_committed"] = {"pass": all(HC.is_tracked(paths.root, f) for f in guard.TRANSFER_FREEZE)}
    checks["gate_G4"] = {"pass": HC.gate_passed(paths, "G4")}
    bad = [f for f, h in crit.get("inputs_sha256", {}).items()
           if (paths.root / f).exists() and HC.sha256_lf(paths.root / f) != h
           and HC.sha256_file(paths.root / f) != h]
    checks["criteria_inputs_unchanged"] = {"pass": not bad, "changed": bad}
    HC.check_noec(fin["H3_noEC"])
    checks["noEC_excludes_4_EC_dev_lines"] = {"pass": True}
    planned = planned_accesses(crit, stage, secondary, h2_arms)
    grant = {}
    for role, scope in planned:
        f = fin[role]
        ok, reason = guard.unblinding_check(f.config_hash, scope, paths.root, paths.access_log)
        grant[f"{f.name}|{scope}"] = reason
    checks["guard_would_grant"] = {"pass": all(v == "ok" for v in grant.values()), "n_planned": len(planned),
                                   "per_access": grant}
    if stage == "B":
        p = Path(paths.out_dir) / "setA" / "summary.json"
        checks["setA_done_first"] = {"pass": p.exists() and json.loads(p.read_text()).get("status") == "complete"}
    if not dry:
        failed = [k for k, v in checks.items() if not v["pass"]]
        if failed:
            raise SystemExit(f"G6 refused to start: {failed}\n{json.dumps(checks, indent=1)[:4000]}")
    return checks


def verdict_arms(crit):
    arms = []
    for bg, spec in crit["backgrounds"].items():
        if spec.get("verdict_eligible"):
            for st, s in spec["studies"].items():
                if s["role"] == "verdict":
                    for e in s["entries"].values():
                        arms += e["arms"]
    return arms


def planned_accesses(crit, stage, secondary, h2_arms):
    out = []
    if stage in ("A", "all"):
        out.append(("T", "setA:A1"))
        if secondary:
            out.append(("headline", "setA:A1"))
    if stage in ("B", "all"):
        arms = verdict_arms(crit)
        out += [("T", f"setB:{a}") for a in arms]
        out += [("T", f"setB:{a}") for a in (h2_arms.get("bowden_wt"), h2_arms.get("bowden_ko")) if a]
        if secondary:
            out += [(r, f"setB:{a}") for r in SECONDARY_MODELS for a in arms]
        st, tp = h2_arms["mair_entry"].split(":")
        m = crit["backgrounds"]["hPSC"]["studies"][st]["entries"][tp]["arms"]
        out += [("H3_noEC", f"setB:{a}") for a in m]
    return out


# ============================================================================ driver
def strip_private(o):
    if isinstance(o, dict):
        return {k: strip_private(v) for k, v in o.items() if not str(k).startswith("_")}
    if isinstance(o, list):
        return [strip_private(v) for v in o]
    return o


def rehearse(crit, W, bank, fin, paths, stage, secondary, h2_arms, b1_setA, anchors, log=print):
    """Run the whole post-access computation once on stand-in observed data (B0 + noise; no access)."""
    rng = np.random.default_rng(1)
    b0 = W.comp["B0"]

    def fake(kind, f, name):
        if kind == "setA":
            return pd.DataFrame(np.add.outer(np.zeros(len(name)), b0.to_numpy()) + rng.normal(0, 0.3, (len(name), len(b0))),
                                index=list(name), columns=b0.index)
        v = b0 + rng.normal(0, 0.3, len(b0))
        return pd.DataFrame({"final": v, "affine": v, "raw": v})

    crit_r = json.loads(json.dumps(crit))
    crit_r["bootstrap"]["n_resamples"] = 20
    Wr = World(paths, crit_r, W.comp)
    led = HC.AccessLedger()
    obs = Observed(paths, led, True, log=log, fake=fake)
    if stage in ("A", "all"):
        sa = evaluate_set_a(crit_r, Wr, bank, obs, fin, None, b1_setA, secondary)
        if sa["errors"]:
            raise RuntimeError(f"rehearsal (set A) failed: {sa['errors']}")
    if stage in ("B", "all"):
        bgr = {bg: evaluate_background(bg, spec, Wr, bank, obs, fin, anchors or {}, True, secondary, log)
               for bg, spec in crit_r["backgrounds"].items() if spec.get("verdict_eligible")}
        errs = {bg: r["errors"] for bg, r in bgr.items() if r["errors"]}
        hyp = hypotheses(crit_r, bgr, bank, obs, fin, Wr, h2_arms, log)
        if errs or hyp["errors"]:
            raise RuntimeError(f"rehearsal (set B) failed: {errs} {hyp['errors']}")
    log("rehearsal on stand-in observed data: every step ran (no access)")


def run_g6(finalists_path, paths: HC.Paths, dry=False, stage="all", threads=4, secondary=True, x_builder=None,
           expression=None, anchors=None, standin=None, spec_override=None, comparators=None, feature_loader=None,
           n_boot=None, log=print, resume=False, rehearsal=True):
    stamp = dt.datetime.now()
    rel = paths.rel(finalists_path)
    fin = HC.load_finalists(finalists_path, paths.root)
    crit = json.loads(Path(paths.criteria).read_text(encoding="utf-8"))
    pop = HC.population(paths, dry)
    h2_arms = real_h2_arms()
    b1_setA = "B1_Eye"
    dry_info = None
    if spec_override is not None:  # tests: a toy criteria-like spec
        crit, h2_arms, b1_setA = spec_override["crit"], spec_override["h2_arms"], spec_override["b1_setA"]
    elif dry:
        crit, h2_arms, dry_info, b1_setA = dry_spec(crit, paths, pop, None)
    if n_boot:
        crit["bootstrap"]["n_resamples"] = int(n_boot)
    checks = preflight(fin, crit, paths, dry, stage, rel, secondary, h2_arms, resume)
    out_dir = Path(paths.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    # per-stage file names: stage A's committed outputs are never modified before stage B's accesses
    HC.write_json(out_dir / f"preflight_{stage}.json", HC.clean_floats(checks))
    if comparators is None:
        if dry:
            lins = {crit["backgrounds"][b]["B1_nearest_column"][3:] for b in ("RPE1", "hPSC")} | {b1_setA[3:]}
            comparators = dry_comparators(paths, pop, sorted(lins))
        else:
            comparators = pd.read_parquet(paths.comparators)
            comparators.index = comparators.index.astype(str)
    W = World(paths, crit, comparators)
    expr = expression or Expression(paths)
    bank = ModelBank(fin, pop, paths, dry, x_builder or real_x_builder, expr, feature_loader, log)
    ledger = HC.AccessLedger(out_dir / f"access_ledger_{stage}.json")
    if dry and standin is None:
        def standin(f, arm):  # a dev line's DepMap profile stands in for a harmonised screen (dev IDs only)
            Y = HC.guarded_targets(f, [arm], paths, ledger, True, scope_label=f"dry:{arm}:{f.config_hash[:8]}")
            return Y.iloc[0].to_frame("final")
    if stage in ("B", "all") and anchors is None:
        try:
            anchors = real_anchor_profiles(crit, paths)
        except Exception as e:  # noqa: BLE001 - the gate then cannot be informative ('fails' blocked)
            log(f"anchor profiles unavailable: {e}")
            anchors = {}
    if not dry and rehearsal:
        rehearse(crit, W, bank, fin, paths, stage, secondary, h2_arms, b1_setA, anchors, log)
    obs = Observed(paths, ledger, dry, standin, log)
    summary = {"gate": "G6", "dry_run": dry, "stage": stage, "timestamp": stamp.isoformat(timespec="seconds"),
               "git_commit": HC.git_head(paths.root), "population": pop.label, "preflight": checks,
               "dry_run_standins": dry_info, "n_boot": W.n_boot, "seed": W.seed}
    # models first (no target access): fold models + CV reproduction check
    for role in ("T", "H3_noEC") + (SECONDARY_MODELS if secondary else ()):
        bank.get(role)
    summary["models"] = bank.info

    if stage in ("A", "all"):
        log("set A (tier A1) ...")
        ref = None
        tref = Path(paths.extra.get("test_dir", out_dir.parent / "test")) / fin["T"].name / "setA_reference.csv"
        if tref.exists():
            ref = pd.read_csv(tref, dtype={"ModelID": str})
            if dry:
                ref = ref[~ref["ModelID"].isin(crit["set_A"]["tier_A1"]["lines"])]
        sa = evaluate_set_a(crit, W, bank, obs, fin, ref, b1_setA, secondary)
        sa["test_reference_file"] = paths.rel(tref) if ref is not None else None
        sa["status"] = "complete" if not sa["errors"] else "complete_with_errors"
        HC.write_json(out_dir / "setA" / "summary.json", HC.clean_floats(sa))
        summary["setA"] = {k: v["descriptive"] for k, v in sa["statistics"].items()}
    if stage in ("B", "all"):
        bgr = {}
        for bg, spec in crit["backgrounds"].items():
            if not spec.get("verdict_eligible"):
                continue
            log(f"set B background {bg} ...")
            bgr[bg] = evaluate_background(bg, spec, W, bank, obs, fin, anchors, frozen_ceiling=not dry,
                                          secondary=secondary, log=log)
            HC.write_json(out_dir / "setB" / f"{bg}.json", HC.clean_floats(strip_private(bgr[bg])))
        hyp = hypotheses(crit, bgr, bank, obs, fin, W, h2_arms, log)
        overall = {q: TR.overall_verdict({bg: r["verdicts"][q]["verdict"] for bg, r in bgr.items()}) for q in PRIMARY}
        declared = {bg: {k: spec.get(k) for k in ("status", "reason")} for bg, spec in crit["backgrounds"].items()
                    if not spec.get("verdict_eligible")}
        sb = {"verdicts": {bg: {q: r["verdicts"][q] for q in PRIMARY} for bg, r in bgr.items()},
              "overall": overall, "declared_backgrounds": declared, "hypotheses": hyp,
              "ceiling_comparison": {bg: {k: {"T_skill": v["T"]["point"], "T_ci": [v["T"]["ci_lo"], v["T"]["ci_hi"]],
                                              "ceiling_skill": v["ceiling"]["ceiling_skill"], "MDE": v["ceiling"]["MDE"],
                                              "skill_over_ceiling": v["skill_over_ceiling"]}
                                          for k, v in r["statistics"].items()} for bg, r in bgr.items()},
              "status": "complete"}
        HC.write_json(out_dir / "setB" / "summary.json", HC.clean_floats(strip_private(sb)))
        summary["setB"] = strip_private({k: sb[k] for k in ("verdicts", "overall", "ceiling_comparison")})
        summary["hypotheses"] = {k: (v.get("result") if isinstance(v, dict) else v) for k, v in hyp.items()
                                 if k in ("H1", "H2", "H3", "H4")}
    summary["accesses"] = ledger.rows
    summary["n_accesses_granted"] = ledger.n_granted()
    HC.write_json(out_dir / f"summary_{stage}.json", HC.clean_floats(strip_private(summary)))
    pd.DataFrame([{"run_id": f"{stamp:%Y%m%d-%H%M%S}_{'dryG6' if dry else 'G6'}_{stage}", "timestamp": summary["timestamp"],
                   "git_commit": summary["git_commit"], "stage": stage, "n_accesses_granted": ledger.n_granted(),
                   "finalists": ";".join(f"{f.name}:{f.config_hash[:8]}" for f in fin.values()),
                   "threads": threads, "peak_rss_mb": HC.peak_rss()[0]}]).to_csv(
        out_dir / f"registry_record_{stage}.csv", index=False)
    log(f"wrote {paths.rel(out_dir / f'summary_{stage}.json')}")
    return summary


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--finalists", default=None)
    ap.add_argument("--dry-run-dev", action="store_true")
    ap.add_argument("--stage", choices=["A", "B", "all"], default=None, help="real run: A then B (default all in dry run)")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--no-secondary", action="store_true", help="skip headline/B2 secondary models (fewer accesses)")
    ap.add_argument("--n-boot", type=int, default=None, help="dry run only: fewer resamples")
    ap.add_argument("--resume", action="store_true")
    a = ap.parse_args(argv)
    try:
        from threadpoolctl import threadpool_limits
        threadpool_limits(a.threads)
    except ImportError:
        pass
    dry = a.dry_run_dev
    paths = HC.Paths()
    if dry:
        paths = paths.with_(out_dir=ROOT / "results" / "dryrun" / "transfer",
                            extra={"test_dir": ROOT / "results" / "dryrun" / "test"})
    else:
        paths = paths.with_(out_dir=ROOT / "results" / "transfer", extra={"test_dir": ROOT / "results" / "test"})
        if a.n_boot:
            raise SystemExit("--n-boot is for the dry run only (the criteria fix 2000 resamples)")
        if a.stage is None:
            raise SystemExit("real run: give --stage A (set A first), then --stage B")
    fpath = Path(a.finalists) if a.finalists else ROOT / "configs" / ("finalists.example.json" if dry else "finalists.json")
    if not fpath.is_absolute():
        fpath = ROOT / fpath
    if not dry and fpath.name != "finalists.json":
        raise SystemExit("the real G6 run reads configs/finalists.json")
    run_g6(fpath, paths, dry, a.stage or "all", a.threads, secondary=not a.no_secondary, n_boot=a.n_boot,
           resume=a.resume)


if __name__ == "__main__":
    main()
