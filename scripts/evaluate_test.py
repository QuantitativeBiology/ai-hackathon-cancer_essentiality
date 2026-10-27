"""Gate G4: the one-time cancer test-set evaluation of the committed finalists (DESIGN §4, §5, §9).

    # G4 (PI): commit configs/finalists.json first (schema: configs/finalists.example.json); clean tree
    ESSML_UNBLIND=1 python scripts/evaluate_test.py --threads 4
    # rehearsal on dev only: dev fold-0 validation lines stand in for the test set
    python scripts/evaluate_test.py --dry-run-dev --finalists configs/finalists.example.json --threads 4

Protocol
--------
0. Refuses to start unless: the work tree is clean (no modified tracked file and no untracked
   non-ignored file), configs/finalists.json is committed, every config hashes to its recorded
   config_hash, results/gates/G3.json passed, the frozen evaluation hashes match, and the guard would
   grant every finalist's test access (``guard.unblinding_check``: ESSML_UNBLIND=1, finalist hash,
   not accessed before). Nothing is fitted if any check fails.
1. Phase A (no test targets): each finalist is refit on all 1022 dev lines (DESIGN §5 "Cancer test
   set: the model refit on all dev lines"; scripts/run.py tuning protocol, seed = split_seed, leakage
   sentinel) and predicts the test lines from features only. Model T additionally gets its 5 dev
   fold models (D7; persisted for G6) and their test predictions, for the set-A reference
   distribution (transfer_criteria.json set_A.tier_A1.ranking). Each finalist runs in its own
   subprocess so fit/predict time and peak RSS are its own. Predictions are written to
   cache/test/preds/ and their sha256 to results/test/predictions_manifest.json BEFORE any target
   request.
2. Rehearsal: the whole post-access pipeline runs once on stand-in targets (B1 predictions + noise,
   dev NaN patterns), so a code error cannot surface after a real access.
3. Phase B: per finalist, in the order headline, T, H3 no-EC, B1, B2, B0: exactly ONE
   ``guard.load_targets(test_ids, config_hash)`` (an in-process ledger refuses a second request
   before the guard sees it; the guard refuses repeats across processes). The tracked access log is
   committed right after each granted access (scripts/heldout_common.py explains why). Targets live
   in memory only, never on disk. Post-access steps run one by one, each guarded by try/except, and
   results are written after each step, core metrics first.
4. Results: results/test/<finalist>/summary.json (+ per_gene.csv.gz, per_line.csv), the combined
   results/test/summary.json, results/test/registry_record.csv and (real run only, after the last
   access) rows appended to results/registry.csv.

Metrics (DESIGN §4; frozen ``essml.evaluate``, single fold, min_obs 50, eval genes, dev-wide gene
classes): P1 by class (+ median, n r > 0.3, frac beating B1), dP1 vs B0/B1/B2, P2a and dP2a, P2b nAP
and AUROC (selective, prevalence <= 0.5). Uncertainty (``essml.heldout``): 95% percentile CIs from
2000 patient resamples stratified by lineage (seed = configs/seeds.json split_seed), identical draws
for every model, so dP1/dP2a vs B1/B2 are paired. Report columns: within-library, within-lineage and
within lineage x library P1 (D16/D20b; scripts/g3_kernel_report.group_centre), per-lineage P1 for
lineages with >= 20 test lines and per-library-group P1 (min_obs 20 inside those strata), and the
published-comparable block (loosely comparable). Test-line library groups = the dev NaN-pattern
groups (configs/library_groups.csv definition) matched on the test line's own NaN pattern.

Crash after an access: predictions are on disk, so ``--resume`` refits nothing and skips finalists
whose results are complete; a finalist whose access was granted but whose core results are missing
cannot be re-accessed (guard) and is reported for a PI decision.
"""
from __future__ import annotations

import os
import sys

if __name__ == "__main__":  # thread env before numpy is imported (DESIGN §10)
    _n = int(sys.argv[sys.argv.index("--threads") + 1]) if "--threads" in sys.argv else 4
    for _v in ["OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS",
               "VECLIB_MAXIMUM_THREADS"]:
        os.environ[_v] = str(_n)

import argparse  # noqa: E402
import datetime as dt  # noqa: E402
import json  # noqa: E402
import subprocess  # noqa: E402
import time  # noqa: E402
import traceback  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
import heldout_common as HC  # noqa: E402
from essml import genesets, guard, heldout  # noqa: E402
from essml import metrics as M  # noqa: E402
from essml import transfer as TR  # noqa: E402
from essml.baselines import GeneMean, LineageMean  # noqa: E402
from essml.evaluate import evaluate  # noqa: E402

MIN_OBS_TEST = 50       # DESIGN §2/§4: >= 50 observed held-out lines
MIN_OBS_CV = 40
MIN_STRATUM_LINES = 20  # DESIGN §4: within-lineage r for lineages with >= 20 held-out lines
B1_KAPPA = 10           # transfer_criteria.json models.B1_nearest (set-A reference comparator)


# ============================================================================ helpers
def _now():
    return dt.datetime.now().isoformat(timespec="seconds")


def finalist_paths(paths: HC.Paths, f: HC.Finalist):
    return Path(paths.pred_dir) / f"{f.name}.parquet", Path(paths.out_dir) / f.name


def manifest_path(paths):
    return Path(paths.out_dir) / "predictions_manifest.json"


def read_manifest(paths):
    p = manifest_path(paths)
    return json.loads(p.read_text()) if p.exists() else {}


def write_manifest_entry(paths, name, entry):
    m = read_manifest(paths)
    m[name] = entry
    HC.write_json(manifest_path(paths), HC.clean_floats(m))


def load_pred(path, sha):
    if HC.sha256_file(path) != sha:
        raise RuntimeError(f"{path}: sha256 differs from the manifest (prediction changed after phase A)")
    df = pd.read_parquet(path)
    df.columns = df.columns.astype(str)
    df.index = df.index.astype(str)
    return df.astype(np.float64)


# ============================================================================ phase 0: preconditions
def preflight(fin: dict, paths: HC.Paths, dry: bool, finalists_rel: str, resume=False) -> dict:
    """Checks before anything is fitted. Real run: every failure is fatal. Dry run: reported."""
    checks = {}
    allow = [paths.rel(paths.out_dir) + "/"] if resume else []
    dirty = HC.strict_dirty(paths.root, allow)
    checks["tree_clean"] = {"pass": not dirty, "dirty": dirty[:20]}
    checks["finalists_committed"] = {"pass": HC.is_tracked(paths.root, finalists_rel), "file": finalists_rel}
    checks["gate_G3"] = {"pass": HC.gate_passed(paths, "G3")}
    bad = HC.frozen_ok(paths)
    checks["frozen_hashes"] = {"pass": not bad, "changed": bad}
    HC.check_noec(fin["H3_noEC"])
    checks["noEC_excludes_4_EC_dev_lines"] = {"pass": True, "lines": HC.EC_DEV}
    done = completed(paths, fin)
    unb = {}
    for role, f in fin.items():
        if f.name in done:
            unb[f.name] = "complete (resume)"
            continue
        ok, reason = guard.unblinding_check(f.config_hash, "test", paths.root, paths.access_log)
        unb[f.name] = reason
    checks["guard_would_grant"] = {"pass": all(v in ("ok", "complete (resume)") for v in unb.values()),
                                   "per_finalist": unb}
    if not dry:
        failed = [k for k, v in checks.items() if not v["pass"]]
        if failed:
            raise SystemExit(f"G4 refused to start: {failed}\n{json.dumps(checks, indent=1)}")
    return checks


def completed(paths, fin):
    out = set()
    for f in fin.values():
        p = Path(paths.out_dir) / f.name / "summary.json"
        if p.exists() and json.loads(p.read_text()).get("status") in ("complete", "complete_with_errors"):
            out.add(f.name)
    return out


# ============================================================================ phase A: fit + predict (no test targets)
def phase_a_one(f: HC.Finalist, paths: HC.Paths, dry: bool, threads: int, feature_loader=None, log=print):
    pop = HC.population(paths, dry)
    pred_path, _ = finalist_paths(paths, f)
    t0 = time.perf_counter()
    r = HC.refit_predict(f, pop, paths, feature_loader)
    pred = r["pred"]
    if not np.isfinite(pred.to_numpy()).all():
        raise ValueError(f"{f.name}: non-finite predictions")
    pred_path.parent.mkdir(parents=True, exist_ok=True)
    pred.astype(np.float32).to_parquet(pred_path)
    entry = {"file": paths.rel(pred_path), "sha256": HC.sha256_file(pred_path), "n_lines": int(pred.shape[0]),
             "n_genes": int(pred.shape[1]), "fit_s": r["fit_s"], "predict_s": r["predict_s"],
             "n_train_ids": r["n_train_ids"], "chosen_params": r["chosen"], "threads": threads,
             "train_ids_sha256": HC.hashlib.sha256(",".join(sorted(r["train_ids"])).encode()).hexdigest(),
             "excluded_from_training": list(map(str, f.cfg.get("exclude_from_training", []))),
             "config_hash": f.config_hash, "timestamp": _now(), "git_commit": HC.git_head(paths.root)}
    if f.role == "T":  # D7 fold models + their held-out predictions (set-A reference distribution)
        models, info = HC.load_or_fit_fold_models(f, pop, paths, dry, feature_loader, log=log)
        X = (feature_loader or HC.default_feature_loader)(f.cfg.get("features", []), pop.heldout_ids)
        fp = HC.predict_folds(models, pop.heldout_ids, pop.lineage, X)
        entry["fold_preds"] = {}
        for k, p in fp.items():
            fpath = Path(paths.pred_dir) / f"{f.name}_fold{k}.parquet"
            p.astype(np.float32).to_parquet(fpath)
            entry["fold_preds"][str(k)] = {"file": paths.rel(fpath), "sha256": HC.sha256_file(fpath)}
        entry["fold_models"] = {"tag": info["tag"], "source": info["source"]}
    rss, how = HC.peak_rss()
    entry.update(peak_rss_mb=rss, rss_method=how, total_s=time.perf_counter() - t0)
    write_manifest_entry(paths, f.name, entry)
    log(f"phase A {f.name}: fit {r['fit_s']:.1f}s predict {r['predict_s']:.2f}s, peak RSS {rss:.0f} MB")
    return entry


def phase_a(fin, paths, dry, threads, isolate, resume, feature_loader=None, log=print, argv_extra=()):
    man = read_manifest(paths)
    for f in fin.values():
        e = man.get(f.name)
        if resume and e and Path(paths.root, e["file"]).exists() and HC.sha256_file(Path(paths.root, e["file"])) == e["sha256"]:
            log(f"phase A {f.name}: predictions present (resume)")
            continue
        if isolate:
            cmd = [sys.executable, str(Path(__file__).resolve()), "--phase-a", f.role, "--threads", str(threads),
                   *argv_extra]
            r = subprocess.run(cmd, cwd=paths.root)
            if r.returncode != 0:
                raise RuntimeError(f"phase A subprocess failed for {f.name}")
        else:
            phase_a_one(f, paths, dry, threads, feature_loader, log)
    return read_manifest(paths)


# ============================================================================ library groups (dev NaN patterns)
def dev_library_patterns(paths, dev_ids):
    """Top-4 dev row NaN patterns (configs/library_groups.csv definition; dev targets via the guard).
    Returns (list of pattern bytes in group order, dev labels Series)."""
    Y = guard.load_targets(list(dev_ids), **paths.guard_kw())
    cols = Y.columns
    mask = Y.isna().to_numpy()
    del Y
    keys = pd.Series([m.tobytes() for m in mask], index=pd.Index(dev_ids))
    top = list(keys.value_counts().index[:4])
    lab = keys.map({k: i for i, k in enumerate(top)}).fillna(4).astype(int)
    return top, lab, cols


def library_labels(truth_all: pd.DataFrame, top, cols):
    T = truth_all.reindex(columns=cols)
    keys = [m.tobytes() for m in T.isna().to_numpy()]
    pos = {k: i for i, k in enumerate(top)}
    return pd.Series([pos.get(k, 4) for k in keys], index=truth_all.index, name="lib_group")


# ============================================================================ post-access steps
def step_core(ctx, f, pred, truth, refs):
    ev = evaluate(pred, truth, ctx["gene_cls"], None, baseline_preds=refs, genes=ctx["eval_genes"],
                  min_obs=ctx["min_obs"])
    ctx["ev"] = ev
    od = ctx["fdir"]
    od.mkdir(parents=True, exist_ok=True)
    pg = ev["per_gene"].drop(columns=["fold"]).round(5)
    pg.to_csv(od / "per_gene.csv.gz", index=False)
    ev["per_line"].drop(columns=["fold"]).round(5).to_csv(od / "per_line.csv", index=False)
    return {"headline": ev["headline"], "summary": {c: {k: v for k, v in d.items() if not k.endswith(("_sd", "_folds"))}
                                                    for c, d in ev["summary"].items()}}


def step_bootstrap(ctx, f, pred, truth, refs):
    b = heldout.bootstrap_test_metrics(pred, truth, ctx["gene_cls"], ctx["W"], ctx["ev"],
                                       baselines={k: v for k, v in refs.items()}, min_obs=ctx["min_obs"],
                                       cache=ctx["boot_cache"])
    return {"n_boot": int(ctx["W"].shape[0]), "seed": ctx["seed"],
            "scheme": "patients resampled with replacement within lineage strata (stratum = modal lineage of the "
                      "patient's lines); line weight = draws of its patient; same draws for every model (paired)",
            "fixed_at_point_estimate": "evaluable genes, class membership, P2b gene subset",
            "metrics": b}


def step_within(ctx, f, pred, truth, refs):
    from g3_kernel_report import group_centre
    cls = ctx["gene_cls"].reindex(truth.columns).fillna(False)["variable"].to_numpy()
    T = truth.to_numpy(np.float64)
    ev_mask = cls & ((~np.isnan(T)).sum(0) >= ctx["min_obs"])
    lin = ctx["lineage"].reindex(truth.index).to_numpy()
    lib = ctx["lib"].reindex(truth.index).astype(int).astype(str).to_numpy()
    groupings = {"within_library": lib, "within_lineage": lin,
                 "within_lin_lib": np.array([a + "|" + b for a, b in zip(lin, lib)])}
    out = {}
    for who, P in [("model", pred)] + [(k, v) for k, v in refs.items()]:
        Pm = P.loc[truth.index, truth.columns].to_numpy(np.float64)
        row = {"raw": float(np.nanmean(M.per_gene_pearson(Pm, T, ctx["min_obs"])[ev_mask]))}
        for name, lab in groupings.items():
            Pc, Tc = group_centre(Pm, T, lab)
            row[name] = float(np.nanmean(M.per_gene_pearson(Pc, Tc, ctx["min_obs"])[ev_mask]))
        out[who] = row
    for b in refs:
        out[f"delta_vs_{b}"] = {k: out["model"][k] - out[b][k] for k in out["model"]}
    out["_definition"] = ("variable genes with >= min_obs observed test lines; truth and prediction centred per gene "
                          "within library group / lineage / lineage x library group of the test lines "
                          "(scripts/g3_kernel_report.group_centre, D16/D20b), then mean per-gene Pearson r")
    out["library_group_counts"] = ctx["lib"].reindex(truth.index).value_counts().sort_index().to_dict()
    return out


def step_strata(ctx, f, pred, truth, refs):
    """DESIGN §4 'within-lineage per-gene r (lineages with >= 20 held-out lines)' and 'results by NaN/library
    group': P1 (variable) computed inside each stratum with min_obs = 20."""
    cls = ctx["gene_cls"].reindex(truth.columns).fillna(False)["variable"].to_numpy()
    out = {}
    for sname, lab in (("lineage", ctx["lineage"].reindex(truth.index)),
                       ("library_group", ctx["lib"].reindex(truth.index).astype(int).astype(str))):
        rows = {}
        for g, ids in lab.groupby(lab).groups.items():
            if len(ids) < MIN_STRATUM_LINES:
                continue
            T = truth.loc[ids].to_numpy(np.float64)
            m = cls & ((~np.isnan(T)).sum(0) >= MIN_STRATUM_LINES)
            r = M.per_gene_pearson(pred.loc[ids, truth.columns].to_numpy(np.float64), T, MIN_STRATUM_LINES)
            row = {"n_lines": int(len(ids)), "n_genes": int(np.isfinite(r[m]).sum()), "P1_variable": float(np.nanmean(r[m]))}
            for b, P in refs.items():
                rb = M.per_gene_pearson(P.loc[ids, truth.columns].to_numpy(np.float64), T, MIN_STRATUM_LINES)
                row[f"dP1_vs_{b}"] = float(np.nanmean((r - rb)[m]))
            rows[str(g)] = row
        out[sname] = rows
    return out


def step_published(ctx, f, pred, truth, refs):
    pg = ctx["ev"]["per_gene"]
    pg = pg.assign(gene=pg["gene"].astype(str))
    sets = pd.read_csv(ctx["paths"].published_sets, dtype=str).dropna(subset=["entrez"])
    out = {"label": "loosely comparable (DESIGN §4); not used for selection"}
    for s, d in sets.groupby("set"):
        sub = pg[pg["gene"].isin(set(d["entrez"]))]
        out[f"rho_{s}"] = float(sub["rho"].mean())
        out[f"n_{s}"] = int(sub["rho"].notna().sum())
    out["n_r_ge_0.2"] = int((pg["r"] >= 0.2).sum())
    shi = pg[(pg["n_pos"] >= 5) & (pg["n_obs"] - pg["n_pos"] >= 5)]
    out["shi_eligible"] = int(len(shi))
    out["shi_n_r_gt_0.2"] = int((shi["r"] > 0.2).sum())
    return out


def step_setA_reference(ctx, f, pred, truth, refs):
    """transfer_criteria.json set_A.tier_A1.ranking: per cancer-test line, the Q3a/Q3b primary statistics of
    model T (Q3a vs B0; Q3b vs B1 of the line's own lineage, kappa 10, all training lines); U = eval genes
    observed in the line minus the scaling halves. 'fold' = D7 (skill per fold model, averaged; primary,
    matches set A), 'refit' = the all-dev refit model (secondary)."""
    if f.role != "T":
        return None
    c = ctx["setA"]
    rows = []
    for line in truth.index:
        obs = truth.loc[line].dropna()
        ent = TR.Entry(name=line, study=line, observed=obs)
        row = {"ModelID": line, "lineage": ctx["lineage"].get(line)}
        b1 = c["b1_means"].loc[row["lineage"]] if row["lineage"] in c["b1_means"].index else c["b0"]
        for q, comp in (("Q3a", c["b0"]), ("Q3b", b1)):
            e = TR.EntryEval(ent, q, c["eval_genes"], c["exclude"], frac_below=c["frac"], neg_eval=c["neg_eval"])
            row[f"{q}_n_universe"], row[f"{q}_valid"] = int(len(e.U)), bool(e.valid)
            if not e.valid:
                row[f"{q}_fold"] = row[f"{q}_refit"] = np.nan
                continue
            cv = e.prep(comp)()
            row[f"{q}_fold"] = float(np.mean([e.prep(p.loc[line])() for p in ctx["fold_preds"].values()])) - cv
            row[f"{q}_refit"] = e.prep(pred.loc[line])() - cv
        rows.append(row)
    df = pd.DataFrame(rows)
    df.round(6).to_csv(ctx["fdir"] / "setA_reference.csv", index=False)
    return {"file": ctx["paths"].rel(ctx["fdir"] / "setA_reference.csv"), "n_lines": int(len(df)),
            "median_Q3a_fold": float(df["Q3a_fold"].median()), "median_Q3b_fold": float(df["Q3b_fold"].median()),
            "definition": step_setA_reference.__doc__.split("\n", 1)[1].strip()}


STEPS = [("core", step_core), ("bootstrap", step_bootstrap), ("within_group", step_within),
         ("strata", step_strata), ("published_comparable", step_published), ("setA_reference", step_setA_reference)]


def evaluate_finalist(ctx, f, pred, truth_all, refs, steps=STEPS, log=print):
    """Post-access pipeline for one finalist; writes <out>/<name>/summary.json after every step."""
    truth = truth_all.loc[:, ctx["eval_genes"]].astype(np.float64)
    fdir = ctx["fdir"] = Path(ctx["out_dir"]) / f.name
    fdir.mkdir(parents=True, exist_ok=True)
    summ = {"finalist": f.name, "role": f.role, "config": f.config_path, "config_hash": f.config_hash,
            "population": ctx["population_label"], "n_lines": int(truth.shape[0]), "min_obs": ctx["min_obs"],
            "baselines_for_deltas": sorted(refs), "status": "running", "errors": {}, "timestamp": _now()}
    for name, fn in steps:
        try:
            res = fn(ctx, f, pred, truth, refs)
            if res is not None:
                summ[name] = res
        except Exception as e:  # noqa: BLE001 - a later step must not lose the earlier results
            summ["errors"][name] = f"{type(e).__name__}: {e}"
            (fdir / f"error_{name}.txt").write_text(traceback.format_exc())
            if name == "core":
                summ["status"] = "core_failed"
                HC.write_json(fdir / "summary.json", HC.clean_floats(summ))
                return summ
        HC.write_json(fdir / "summary.json", HC.clean_floats(summ))
    summ["status"] = "complete" if not summ["errors"] else "complete_with_errors"
    HC.write_json(fdir / "summary.json", HC.clean_floats(summ))
    log(f"  {f.name}: P1 var {summ['core']['headline'].get('P1_variable', np.nan):.4f}"
        + (f"; errors {list(summ['errors'])}" if summ["errors"] else ""))
    return summ


# ============================================================================ phase B
def build_context(fin, paths, pop, man, n_boot, min_obs, dry):
    eval_genes = pd.Index(pd.read_csv(paths.eval_genes, dtype={"entrez": str})["entrez"])
    gene_cls = genesets.load_gene_classes(paths.gene_classes)
    gc_raw = pd.read_csv(paths.gene_classes, dtype={"entrez": str}).set_index("entrez")
    seed = int(json.loads(Path(paths.seeds).read_text())["split_seed"])
    ho = pop.heldout_ids
    W = heldout.patient_bootstrap_weights(pop.lineage.reindex(ho), pop.patient.reindex(ho), n_boot, seed)
    sp = pd.read_csv(paths.splits, index_col=0)
    dev_all = sp.index[sp["role"] == "dev"].astype(str)
    top, dev_lab, cols = dev_library_patterns(paths, dev_all)
    lib_check = None
    if Path(paths.library_groups).exists():
        ref = pd.read_csv(paths.library_groups, index_col=0)["lib_group"].reindex(dev_lab.index)
        lib_check = bool((ref.astype(int) == dev_lab).all())
    # set-A reference comparators: B0 and B1 (kappa 10) on the training population (all dev / dry: other folds)
    Ytr = guard.load_targets(list(pop.train_ids), **paths.guard_kw()).astype(np.float64)
    b0 = pd.Series(GeneMean().fit(Ytr).mu_, index=Ytr.columns.astype(str))
    b1 = LineageMean(B1_KAPPA).fit(Ytr, pop.lineage.reindex(Ytr.index)).lineage_means()
    b1.columns = b1.columns.astype(str)
    del Ytr
    cs = json.loads(Path(paths.control_split).read_text()) if Path(paths.control_split).exists() else None
    setA = {"b0": b0, "b1_means": b1, "eval_genes": list(eval_genes),
            "exclude": TR.scaling_halves(cs) if cs else set(), "neg_eval": cs["neg_evaluation"] if cs else [],
            "frac": gc_raw["frac_below"].astype(float) if "frac_below" in gc_raw else pd.Series(dtype=float)}
    fold_preds = {}
    t = fin["T"]
    for k, e in man.get(t.name, {}).get("fold_preds", {}).items():
        fold_preds[int(k)] = load_pred(Path(paths.root, e["file"]), e["sha256"])
    return {"eval_genes": eval_genes, "gene_cls": gene_cls, "seed": seed, "W": W, "boot_cache": {},
            "lineage": pop.lineage, "lib_top": top, "lib_cols": cols, "library_groups_reproduced": lib_check,
            "setA": setA, "fold_preds": fold_preds, "min_obs": min_obs, "paths": paths,
            "population_label": pop.label, "out_dir": paths.out_dir, "dry": dry}


def rehearse(fin, preds, ctx, pop, paths, log=print):
    """Run every post-access step on stand-in targets (no target access) so code errors surface now."""
    rng = np.random.default_rng(0)
    base = preds["B1"]
    stand = base + rng.normal(0, 0.3, base.shape)
    stand = stand.reindex(columns=ctx["lib_cols"]).fillna(0.0)
    mask = rng.random(stand.shape) < 0.02
    stand = stand.mask(mask)
    rctx = dict(ctx, out_dir=Path(paths.pred_dir).parent / "rehearsal", W=ctx["W"][:20], boot_cache={})
    rctx["lib"] = library_labels(stand, ctx["lib_top"], ctx["lib_cols"])
    for f in fin.values():
        s = evaluate_finalist(rctx, f, preds[f.role], stand, refs_for(f, preds), log=lambda *_: None)
        if s["status"] != "complete":
            raise RuntimeError(f"rehearsal failed for {f.name}: {s['errors']} (see {rctx['out_dir']})")
    log("rehearsal on stand-in targets: every post-access step ran")


def refs_for(f, preds):
    return {b: preds[b] for b in ("B0", "B1", "B2") if b != f.role}


def cv_comparison(f, pred, truth, ctx, paths):
    """Dry run only: compare with the registered CV run's fold-0 numbers and OOF predictions."""
    if not f.cv_run_id:
        return None
    out = {"cv_run_id": f.cv_run_id}
    pf = Path(paths.root) / "results" / "runs" / f.cv_run_id / "per_fold.csv"
    if pf.exists():
        cv = pd.read_csv(pf)
        cv0 = cv[(cv["fold"] == 0)].set_index("class")
        out["cv_fold0"] = {c: {k: float(cv0.loc[c, k]) for k in ("P1", "P2b_nap_prev05") if k in cv0}
                           for c in cv0.index}
        if "dP1_vs_B1" in cv0:
            out["cv_fold0"]["variable"]["dP1_vs_B1"] = float(cv0.loc["variable", "dP1_vs_B1"])
    ev40 = evaluate(pred, truth.loc[:, ctx["eval_genes"]], ctx["gene_cls"], None,
                    baseline_preds=refs_for(f, ctx["preds"]), genes=ctx["eval_genes"], min_obs=MIN_OBS_CV)
    out["dryrun_min_obs40"] = {c: {"P1": ev40["summary"][c]["P1"], "P2b_nap_prev05": ev40["summary"][c]["P2b_nap_prev05"]}
                               for c in ev40["summary"]}
    out["dryrun_min_obs40"]["variable"]["dP1_vs_B1"] = ev40["summary"]["variable"].get("dP1_vs_B1")
    oof = Path(paths.cv_cache) / f.cv_run_id / "pred.parquet"
    if oof.exists():
        P = pd.read_parquet(oof)
        P.columns = P.columns.astype(str)
        ids = pred.index.intersection(P.index)
        g = pred.columns.intersection(P.columns)
        out["max_abs_diff_vs_cv_oof_pred"] = float(np.abs(pred.loc[ids, g].to_numpy() - P.loc[ids, g].to_numpy()).max())
        out["n_lines_compared"] = int(len(ids))
    return out


def phase_b(fin, paths, pop, man, n_boot, dry, log=print, min_obs=MIN_OBS_TEST, steps=STEPS, rehearsal=True):
    preds = {f.name: load_pred(Path(paths.root, man[f.name]["file"]), man[f.name]["sha256"]) for f in fin.values()}
    preds_role = {f.role: preds[f.name] for f in fin.values()}
    ctx = build_context(fin, paths, pop, man, n_boot, min_obs, dry)
    ctx["preds"] = preds_role
    if rehearsal and len(completed(paths, fin)) < len(fin):
        rehearse(fin, preds_role, ctx, pop, paths, log)
    ledger = HC.AccessLedger(Path(paths.out_dir) / "access_ledger.json")
    if (Path(paths.out_dir) / "access_ledger.json").exists():
        ledger.rows = json.loads((Path(paths.out_dir) / "access_ledger.json").read_text())
    done = completed(paths, fin)
    results, previous = {}, None
    for role in HC.ROLE_ORDER:
        f = fin[role]
        if f.name in done:
            log(f"{f.name}: complete (resume), not accessed again")
            results[f.name] = json.loads((Path(paths.out_dir) / f.name / "summary.json").read_text())
            continue
        if ledger.n_granted(f.config_hash) and not dry:
            raise SystemExit(f"{f.name}: test access already granted but results incomplete; the guard refuses a "
                             "second access. PI decision needed (see results/test/access_ledger.json).")
        truth_all = HC.guarded_targets(f, pop.heldout_ids, paths, ledger, dry,
                                       scope_label=None if not dry else f"test-standin:dev-fold0")
        truth_all.columns = truth_all.columns.astype(str)
        if previous is not None:  # every access must return the same matrix
            same = np.array_equal(previous.to_numpy(), truth_all.to_numpy(), equal_nan=True)
            if not same:
                raise RuntimeError("test targets differ between accesses")
        previous = truth_all
        ctx["lib"] = library_labels(truth_all, ctx["lib_top"], ctx["lib_cols"])
        res = evaluate_finalist(ctx, f, preds_role[role], truth_all, refs_for(f, preds_role), steps, log)
        if dry:
            try:
                res["cv_comparison"] = cv_comparison(f, preds_role[role], truth_all, ctx, paths)
            except Exception as e:  # noqa: BLE001
                res["cv_comparison"] = {"error": f"{type(e).__name__}: {e}"}
            HC.write_json(Path(paths.out_dir) / f.name / "summary.json", HC.clean_floats(res))
        results[f.name] = res
    del previous
    return results, ledger, ctx


# ============================================================================ outputs
def headline_table(results, fin, man):
    rows = []
    for f in fin.values():
        r = results.get(f.name, {})
        h = (r.get("core") or {}).get("headline", {})
        bm = (r.get("bootstrap") or {}).get("metrics", {})
        e = man.get(f.name, {})
        row = {"finalist": f.name, "role": f.role, "config_hash": f.config_hash, "status": r.get("status"),
               "fit_s": e.get("fit_s"), "predict_s": e.get("predict_s"), "threads": e.get("threads"),
               "peak_rss_mb": e.get("peak_rss_mb"), "n_train_ids": e.get("n_train_ids")}
        for k in ("P1_variable", "P1_all", "P1_selective", "P1_rare_selective", "dP1_vs_B1_variable",
                  "dP1_vs_B2_variable", "frac_beat_B1_variable", "n_r_gt_0.3_variable", "P2a", "dP2a_vs_B0",
                  "P2b_nap_prev05", "P2b_auroc_prev05"):
            row[k] = h.get(k)
            b = bm.get(k)
            if b:
                row[k + "_ci"] = [b["ci_lo"], b["ci_hi"]]
        w = r.get("within_group", {}).get("model", {})
        row.update({f"P1_{k}": w.get(k) for k in ("within_library", "within_lineage", "within_lin_lib")})
        rows.append(row)
    return rows


def registry_rows(results, fin, man, paths, dry, stamp, threads, ctx):
    from run import REG_COLS
    rows = []
    for f in fin.values():
        r, e = results.get(f.name, {}), man.get(f.name, {})
        h = (r.get("core") or {}).get("headline", {})
        row = {"run_id": f"{stamp:%Y%m%d-%H%M%S}_{'dryG4' if dry else 'G4test'}_{f.name}_{f.config_hash[:8]}",
               "timestamp": stamp.isoformat(timespec="seconds"), "name": f"{'dryG4' if dry else 'G4test'}_{f.name}",
               "model": f.cfg["model"], "config": f.config_path, "config_hash": f.config_hash,
               "git_commit": HC.git_head(paths.root), "dirty_tree": 0,
               "splits_sha256": HC.sha256_file(paths.splits), "roles_used": ctx.get("roles_used", ""),
               "n_train_ids": e.get("n_train_ids"), "leaky": 0, "status": r.get("status"),
               "n_lines": r.get("n_lines"), "n_genes_eval": len(ctx["eval_genes"]),
               "chosen_params": json.dumps(e.get("chosen_params", {})), "fit_s": e.get("fit_s"),
               "predict_s": e.get("predict_s"), "total_s": e.get("total_s"), "threads": e.get("threads", threads),
               "peak_rss_mb": e.get("peak_rss_mb"), "rss_method": e.get("rss_method"),
               "eval_genes_sha256": HC.sha256_file(paths.eval_genes),
               "gene_classes_sha256": HC.sha256_file(paths.gene_classes), "tree_clean": int(not dry)}
        row.update({k: v for k, v in h.items() if k in REG_COLS})
        rows.append({k: row.get(k, "") for k in REG_COLS})
    return rows


def run_g4(finalists_path, paths: HC.Paths, dry=False, threads=4, n_boot=2000, isolate=True, resume=False,
           feature_loader=None, log=print, steps=STEPS, argv_extra=(), min_obs=MIN_OBS_TEST):
    stamp = dt.datetime.now()
    rel = paths.rel(finalists_path)
    fin = HC.load_finalists(finalists_path, paths.root)
    checks = preflight(fin, paths, dry, rel, resume)
    Path(paths.out_dir).mkdir(parents=True, exist_ok=True)
    HC.write_json(Path(paths.out_dir) / "preflight.json", HC.clean_floats(checks))
    pop = HC.population(paths, dry)
    log(f"G4 {'DRY RUN ' if dry else ''}{pop.label}: {len(pop.train_ids)} training lines, "
        f"{len(pop.heldout_ids)} held-out lines")
    man = phase_a(fin, paths, dry, threads, isolate, resume, feature_loader, log, argv_extra)
    results, ledger, ctx = phase_b(fin, paths, pop, man, n_boot, dry, log, min_obs=min_obs, steps=steps)
    roles = guard.roles_of(list(pop.train_ids) + list(pop.heldout_ids), paths.splits)
    ctx["roles_used"] = ";".join(f"{k}:{v}" for k, v in sorted(roles.value_counts().items()))
    summary = {"gate": "G4", "dry_run": dry, "population": pop.label, "timestamp": stamp.isoformat(timespec="seconds"),
               "git_commit": HC.git_head(paths.root), "finalists_file": rel, "preflight": checks,
               "n_boot": n_boot, "min_obs": min_obs, "library_groups_reproduced": ctx["library_groups_reproduced"],
               "accesses": ledger.rows, "n_accesses_granted": ledger.n_granted(),
               "table": headline_table(results, fin, man),
               "per_finalist": {n: paths.rel(Path(paths.out_dir) / n / "summary.json") for n in results}}
    if dry:
        summary["cv_comparison"] = {n: r.get("cv_comparison") for n, r in results.items()}
    HC.write_json(Path(paths.out_dir) / "summary.json", HC.clean_floats(summary))
    reg = registry_rows(results, fin, man, paths, dry, stamp, threads, ctx)
    pd.DataFrame(reg).to_csv(Path(paths.out_dir) / "registry_record.csv", index=False)
    if not dry:  # after the last access: the registry is a tracked file
        for row in reg:
            HC.RUN.write_registry(row, paths.registry)
    log(f"wrote {paths.rel(Path(paths.out_dir) / 'summary.json')}")
    return summary


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--finalists", default=None, help="default configs/finalists.json (dry run: the example)")
    ap.add_argument("--dry-run-dev", action="store_true", help="dev fold 0 stands in for the test set (dev IDs only)")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--no-isolate", action="store_true", help="fit finalists in this process (timing less clean)")
    ap.add_argument("--phase-a", default=None, help=argparse.SUPPRESS)
    a = ap.parse_args(argv)
    try:
        from threadpoolctl import threadpool_limits
        threadpool_limits(a.threads)
    except ImportError:
        pass
    dry = a.dry_run_dev
    paths = HC.Paths()
    if dry:
        paths = paths.with_(out_dir=ROOT / "results" / "dryrun" / "test", pred_dir=ROOT / "cache" / "dryrun" / "test" / "preds")
    fpath = Path(a.finalists) if a.finalists else ROOT / "configs" / ("finalists.example.json" if dry else "finalists.json")
    if not fpath.is_absolute():
        fpath = ROOT / fpath
    if not dry and fpath.name != "finalists.json":
        raise SystemExit("the real G4 run reads configs/finalists.json (the guard only honours that file)")
    extra = (["--dry-run-dev"] if dry else []) + ["--finalists", str(fpath)]
    if a.phase_a:
        fin = HC.load_finalists(fpath, paths.root)
        phase_a_one(fin[a.phase_a], paths, dry, a.threads)
        return
    run_g4(fpath, paths, dry, a.threads, a.n_boot, isolate=not a.no_isolate, resume=a.resume, argv_extra=extra)


if __name__ == "__main__":
    main()
