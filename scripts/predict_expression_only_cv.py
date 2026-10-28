"""Post-hoc dev-CV diagnostic (G7, after the REPORT red team): m5_avg applied with EXPRESSION ONLY.

    python scripts/predict_expression_only_cv.py            # -> results/predict_expression_only_cv.json

scripts/predict.py reads expression only. With ``--model m5_avg`` its genomic inputs (CN, damaging / hotspot
mutations, fusions) are absent = diploid wild type and its metadata (lineage, sex, age) unknown, through
``StackM5.predict_from_blocks`` (the G6 external path). The G4 number of m5_avg used real genomics. This script
scores that shipped configuration on the dev CV folds:

1. The persisted G6 fold models cache/models/m5_avg_6dbebcf4 (fold k trained on the other dev folds) are loaded and
   checked against the CV run 20261027-161041_m5_avg_6dbebcf4 (config hash, fold signature, training lines).
   m5_avg ("average" mode) stores its training data only; its three inputs (M1, M2, M3) are refitted on the fold's
   training lines inside the fold model at first prediction, exactly as in CV, G6 and predict.py
   (``StackM5.fitted_inputs``; KEEP_INPUTS). No hyperparameter, weight or input changes.
2. Self-check: with ALL of a held-out line's own raw blocks (E, C, M, H, meta, gene-level CN / damaging; as the G6
   self-check, evaluate_transfer.m5_blocks_for), the block path must reproduce the registered CV out-of-fold
   predictions (cache/runs/<cv run>/pred.parquet) for every held-out line of the fold.
3. Expression only: the held-out lines are predicted exactly as predict.py predicts a new line: IDs
   ``input::<line>``, X = {"E": the line's DepMap log2(TPM+1) row} (= cache/features/E.parquet, what a user
   supplying the DepMap profile gives; --harmonise none), lineage None.
4. Scoring with the frozen essml.evaluate.evaluate (as scripts/run.py: dev-wide gene classes, eval genes, min_obs
   40); paired fold differences vs m1_E_resB0 and m5_avg (registered per_fold.csv; scripts/g2_report.paired);
   NRAS / KRAS / BRAF / PIK3CA hotspot AUPRC with essml.mutdeps.fold_metrics and the declared metric_params of
   configs/mutation_dependencies.json (as docs/mutation_dependencies.md).

Dev lines only: targets via essml.guard (dev IDs). No test or external target is read.
"""
from __future__ import annotations

import os
import sys

for _v in ["OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"]:
    os.environ.setdefault(_v, "6")

import argparse  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
from essml import genesets, guard, modelstore  # noqa: E402
from essml import mutdeps as MD  # noqa: E402
from essml.evaluate import evaluate  # noqa: E402

TAG = "m5_avg_6dbebcf4"
M5_RUN = "20261027-161041_m5_avg_6dbebcf4"
M1_RUN = "20261026-222117_m1_E_resB0_91e52f9d"
CONFIG = ROOT / "configs" / "runs" / "m5_avg.json"
OUT = ROOT / "results" / "predict_expression_only_cv.json"
PRED_CACHE = ROOT / "cache" / "predict_expression_only_cv"
HOTSPOT_PAIRS = ["NRAS_hotspot", "KRAS_hotspot", "BRAF_hotspot", "PIK3CA_hotspot"]
CV_REPRO_TOL = 1e-3          # heldout_common.CV_REPRO_TOL (pred.parquet is float32)
DELTA = 0.005                # DESIGN §4
CLASSES = ["variable", "all", "selective", "rare_selective", "common_essential", "nonessential"]


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def git_commit():
    import subprocess
    h = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"], cwd=ROOT,
                           capture_output=True, text=True).stdout.strip()
    return h, bool(dirty)


def dev_frame():
    sp = pd.read_csv(ROOT / "configs" / "splits.csv", index_col=0)
    sp.index = sp.index.astype(str)
    return sp[sp["role"] == "dev"]


def check_fold_models(models, dev, log):
    """The persisted fold models are those of the registered m5_avg CV run (config, folds, training lines)."""
    folds = dev["cv_fold"].astype(int)
    meta = json.loads((modelstore.MODEL_DIR / TAG / "meta.json").read_text())
    cfg = json.loads(CONFIG.read_text())
    summ = json.loads((ROOT / "results" / "runs" / M5_RUN / "summary.json").read_text())
    ch = hashlib.sha256(json.dumps(cfg, sort_keys=True).encode()).hexdigest()
    sig = hashlib.sha256(json.dumps(sorted(map(str, folds.index)) + [int(x) for x in folds.sort_index()]).encode())
    out = {"config_hash_file": ch, "config_hash_meta": meta["config_hash"], "config_hash_cv_run": summ["config_hash"],
           "folds_sha256_meta": meta["folds_sha256"], "folds_sha256_splits": sig.hexdigest(), "dry": meta["dry"],
           "per_fold": {}}
    ok = ch == meta["config_hash"] == summ["config_hash"] and meta["folds_sha256"] == sig.hexdigest() and not meta["dry"]
    for k, m, mt in models:
        exp = set(folds.index[folds != k])
        tr = set(map(str, mt["train_ids"]))
        same = tr == exp and set(map(str, m.Y_.index)) == exp and m.mode == "average"
        out["per_fold"][str(k)] = {"n_train": len(tr), "train_ids_equal_other_folds": same,
                                   "fold_pkl_sha256": sha256(modelstore.MODEL_DIR / TAG / f"fold{k}.pkl")}
        ok &= same
    out["ok"] = bool(ok)
    if not ok:
        raise RuntimeError(f"persisted fold models do not match {M5_RUN}: {json.dumps(out, indent=1)}")
    log(f"fold models {TAG} match {M5_RUN} (config, fold signature, training lines)")
    return out


def own_blocks(ids):
    """Every raw block an M5 input reads for DepMap lines (evaluate_transfer.m5_blocks_for, vectorised)."""
    from essml.io import load_matrix
    from essml.models.g3_common import load_blocks
    blocks = {b: v.reindex(pd.Index(ids).intersection(v.index)) for b, v in
              load_blocks(["E", "C", "M", "H", "meta"]).items()}
    for key, name in (("_expression", "expression"), ("_cn", "cn"), ("_damaging", "damaging")):
        mat = load_matrix(name)
        blocks[key] = mat.loc[mat.index.intersection(pd.Index(ids))]
    return blocks


def expression_only_X(E, ids):
    """Exactly predict.py's input for a new line: IDs 'input::<id>', X = {'E': expression} (X_BUILDER)."""
    new = pd.Index([f"input::{i}" for i in ids])
    return new, {"E": E.loc[ids].astype(np.float64).set_axis(new, axis=0)}


def predict_folds(models, dev, log):
    from essml.io import load_matrix
    from essml.models.m5_stack import StackM5
    StackM5.KEEP_INPUTS = True
    folds = dev["cv_fold"].astype(int)
    E = load_matrix("expression")
    E.columns, E.index = E.columns.astype(str), E.index.astype(str)
    reg = pd.read_parquet(ROOT / "cache" / "runs" / M5_RUN / "pred.parquet")
    reg.columns, reg.index = reg.columns.astype(str), reg.index.astype(str)
    parts, checks, timing = [], {}, {}
    PRED_CACHE.mkdir(parents=True, exist_ok=True)
    for k, m, mt in models:
        ho = pd.Index(sorted(folds.index[folds == k]))
        cp = PRED_CACHE / f"fold{k}.parquet"
        ck = PRED_CACHE / f"fold{k}_check.json"
        if cp.exists() and ck.exists():
            parts.append(pd.read_parquet(cp))
            checks[str(k)] = json.loads(ck.read_text())
            log(f"fold {k}: loaded {cp.relative_to(ROOT)}")
            continue
        t0 = time.perf_counter()
        m.fitted_inputs()  # the fold model's inputs, refitted on its training lines as in CV (deterministic)
        t_fit = time.perf_counter() - t0
        for e in m.input_fit_log_:
            if set(e["ids"]) & set(ho) or (guard.roles_of(e["ids"]) != "dev").any():
                raise RuntimeError(f"fold {k}: input {e['input']} saw held-out / non-dev lines")
        # (2) self-check: all own blocks -> registered CV OOF
        full = m.predict_from_blocks(ho, own_blocks(ho))
        full.columns = full.columns.astype(str)
        g = full.columns.intersection(reg.columns)
        d_full = np.abs(full.loc[ho, g].to_numpy(np.float64) - reg.loc[ho, g].to_numpy(np.float64))
        chk = {"n_lines": int(len(ho)), "n_genes": int(len(g)), "max_abs_diff_vs_registered_oof": float(d_full.max()),
               "lines_checked": "all held-out lines of the fold", "tol": CV_REPRO_TOL,
               "input_refit_s": round(t_fit, 1), "input_fit_ids_n": {e["input"]: len(e["ids"]) for e in m.input_fit_log_}}
        chk["ok"] = bool(len(g) == reg.shape[1] and chk["max_abs_diff_vs_registered_oof"] <= CV_REPRO_TOL)
        log(f"fold {k}: own-blocks path vs registered OOF max |diff| {chk['max_abs_diff_vs_registered_oof']:.3g} "
            f"({len(ho)} lines; inputs refitted in {t_fit:.0f} s)")
        if not chk["ok"]:
            raise RuntimeError(f"fold {k}: self-check failed {chk}")
        # (3) expression only, exactly as predict.py
        new, X = expression_only_X(E, ho)
        P = m.predict(new, lineage=None, X=X)  # StackM5.predict -> predict_from_blocks (dict X)
        if not np.isfinite(P.to_numpy()).all():
            raise RuntimeError(f"fold {k}: non-finite expression-only predictions")
        P.index = ho
        P.columns = P.columns.astype(str)
        d_eo = (P.loc[ho, g] - full.loc[ho, g]).abs()
        chk["expression_only_vs_own_blocks_mean_abs_diff"] = float(d_eo.to_numpy().mean())
        P.astype(np.float32).to_parquet(cp)
        ck.write_text(json.dumps(chk, indent=1))
        parts.append(P.astype(np.float32))
        checks[str(k)] = chk
        timing[str(k)] = round(time.perf_counter() - t0, 1)
        m.models_ = None  # free the kept inputs
    pred = pd.concat(parts).astype(np.float64)
    return pred, checks


def score(pred, dev, log):
    """Frozen evaluation exactly as scripts/run.py (no baseline predictions: P1 / P2 are unaffected by them)."""
    folds = dev["cv_fold"].astype(int)
    ids = pd.Index(sorted(dev.index))
    Y = guard.load_targets(list(ids)).astype(np.float64)  # dev IDs only
    Y.columns = Y.columns.astype(str)
    eval_genes = pd.Index(pd.read_csv(ROOT / "configs" / "eval_genes.csv", dtype={"entrez": str})["entrez"])
    gene_cls = genesets.load_gene_classes()
    min_obs = int(json.loads(CONFIG.read_text()).get("min_obs", 40))
    ev = evaluate(pred.loc[ids], Y.loc[ids], gene_cls, folds.loc[ids], genes=eval_genes, min_obs=min_obs)
    # consistency: the same call on the registered m5_avg OOF reproduces its registered per_fold.csv
    reg = pd.read_parquet(ROOT / "cache" / "runs" / M5_RUN / "pred.parquet").astype(np.float64)
    reg.columns, reg.index = reg.columns.astype(str), reg.index.astype(str)
    ev_reg = evaluate(reg.loc[ids], Y.loc[ids], gene_cls, folds.loc[ids], genes=eval_genes, min_obs=min_obs)
    pf_reg = pd.read_csv(ROOT / "results" / "runs" / M5_RUN / "per_fold.csv")
    a = ev_reg["per_fold"].set_index(["fold", "class"])["P1"].sort_index()
    b = pf_reg.set_index(["fold", "class"])["P1"].sort_index()
    repro = float(np.abs(a.to_numpy() - b.reindex(a.index).to_numpy()).max())
    log(f"evaluate() on the registered m5_avg OOF reproduces its per_fold.csv P1 (max |diff| {repro:.2g})")
    if repro > 1e-5:
        raise RuntimeError("evaluation call does not reproduce the registered per_fold.csv")
    return ev, Y, repro


def paired(a, b):
    from g2_report import paired as _p
    return _p(a, b)


def per_fold_table(path_or_df):
    pf = pd.read_csv(path_or_df) if not isinstance(path_or_df, pd.DataFrame) else path_or_df
    return pf.set_index(["class", "fold"]).sort_index()


def hotspot_auprc(pred, Y, dev):
    """docs/mutation_dependencies.md AUPRC (all-line AP, effect < -0.5) for the four hotspot pairs."""
    from mutation_dependencies import raw_status_mats
    decl = json.loads((ROOT / "configs" / "mutation_dependencies.json").read_text())
    pairs = pd.read_csv(ROOT / "configs" / "oncogene_pairs.csv", dtype={"dependency_entrez": str,
                                                                         "biomarker_entrez": str}).set_index("pair_id")
    md = json.loads((ROOT / "results" / "mutation_dependencies.json").read_text())
    folds = dev["cv_fold"].astype(int)
    hot, dam, cn = raw_status_mats()
    ids = Y.index
    m1 = pd.read_parquet(ROOT / "cache" / "runs" / M1_RUN / "pred.parquet")
    m5 = pd.read_parquet(ROOT / "cache" / "runs" / M5_RUN / "pred.parquet")
    for p in (m1, m5):
        p.columns, p.index = p.columns.astype(str), p.index.astype(str)
    preds = {"m5_avg_expression_only": pred, "m1_E_resB0": m1, "m5_avg_with_genomics": m5}
    out = {}
    for pid in HOTSPOT_PAIRS:
        pr = pairs.loc[pid]
        g = pr["dependency_entrez"]
        s = MD.biomarker_status(pr["biomarker_type"], pr["biomarker_entrez"], ids, hot, dam, cn)
        res = {}
        for name, P in preds.items():
            rows = []
            for k in sorted(folds.unique()):
                te = folds.index[folds == k]
                rows.append(MD.fold_metrics(P.loc[te, g].to_numpy(np.float64), Y.loc[te, g].to_numpy(),
                                            s.loc[te].to_numpy(), **decl["metric_params"]))
            ap = np.array([r["ap_all"] for r in rows], np.float64)
            res[name] = {"ap_all_mean": float(np.nanmean(ap)), "ap_all_folds": np.round(ap, 5).tolist(),
                         "r_all_mean": float(np.nanmean([r["r_all"] for r in rows])),
                         "auroc_bm_mean": float(np.nanmean([r["auroc_bm"] for r in rows]))}
        for name in ("m1_E_resB0", "m5_avg_with_genomics"):
            ci = MD.paired_ci(np.asarray(res["m5_avg_expression_only"]["ap_all_folds"])
                              - np.asarray(res[name]["ap_all_folds"]))
            res["m5_avg_expression_only"][f"delta_ap_vs_{name}"] = ci
        # the recomputation for the registered runs must equal results/mutation_dependencies.json
        for name, key in (("m1_E_resB0", "m1_E_resB0"), ("m5_avg_with_genomics", "m5_avg")):
            ref = md["per_pair"][pid][key]["ap_all"]["mean"]
            if abs(ref - res[name]["ap_all_mean"]) > 1e-4:
                raise RuntimeError(f"{pid} {name}: AUPRC {res[name]['ap_all_mean']} != committed {ref}")
        res["dependency_gene"], res["n_mutant_dev"] = pr["dependency_gene"], int(pr["n_mutant_dev"])
        out[pid] = res
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--threads", type=int, default=6)
    a = ap.parse_args(argv)
    try:
        from threadpoolctl import threadpool_limits
        threadpool_limits(a.threads)
    except ImportError:
        pass
    np.random.seed(20261026)
    t0 = time.perf_counter()
    log = lambda s: print(s, flush=True)  # noqa: E731
    commit, dirty = git_commit()
    dev = dev_frame()
    models = modelstore.load_fold_models(TAG)
    fm = check_fold_models(models, dev, log)
    pred, checks = predict_folds(models, dev, log)
    ev, Y, repro = score(pred, dev, log)
    pf = ev["per_fold"].set_index(["class", "fold"]).sort_index()
    m1 = per_fold_table(ROOT / "results" / "runs" / M1_RUN / "per_fold.csv")
    m5 = per_fold_table(ROOT / "results" / "runs" / M5_RUN / "per_fold.csv")
    folds = sorted(pf.loc["variable"].index)
    by_class = {}
    for c in CLASSES:
        eo = pf.loc[c].loc[folds]
        by_class[c] = {"P1_folds": eo["P1"].round(6).tolist(), "P1_mean": float(eo["P1"].mean()),
                       "P2a_mean": float(eo["P2a"].mean()), "P2b_nap_prev05_mean": float(eo["P2b_nap_prev05"].mean()),
                       "P2b_auroc_prev05_mean": float(eo["P2b_auroc_prev05"].mean()),
                       "m1_E_resB0_P1_mean": float(m1.loc[c].loc[folds, "P1"].mean()),
                       "m5_avg_with_genomics_P1_mean": float(m5.loc[c].loc[folds, "P1"].mean()),
                       "delta_P1_vs_m1_E_resB0": paired(eo["P1"].to_numpy(), m1.loc[c].loc[folds, "P1"].to_numpy()),
                       "delta_P1_vs_m5_avg_with_genomics": paired(eo["P1"].to_numpy(),
                                                                  m5.loc[c].loc[folds, "P1"].to_numpy())}
    v = by_class["variable"]["delta_P1_vs_m1_E_resB0"]
    beats = bool(v["mean"] > DELTA and v["mean"] > v["se"])
    hs = hotspot_auprc(pred, Y, dev)
    res = {
        "post_hoc": True,
        "note": ("Post-hoc dev-CV diagnostic added at G7 after the REPORT red team (not pre-registered; dev lines only; "
                 "no test-set or external target read). It scores the configuration scripts/predict.py actually runs: "
                 "m5_avg trained with genomics (its persisted G6 fold models) applied to held-out dev lines with "
                 "expression only (genomic inputs absent = diploid wild type, metadata unknown; "
                 "StackM5.predict_from_blocks with X = {'E': expression}, IDs 'input::<line>'). The G4 test number of "
                 "m5_avg (0.319) used real genomics. The predict.py default rule (DESIGN section 4: dP1 variable vs "
                 "m1_E_resB0 > 0.005 and > 1 SE of the paired fold differences) was fixed by the PI before this "
                 "number was seen."),
        "git_commit": commit, "dirty_tree": dirty,
        "fold_models": {"tag": TAG, "dir": f"cache/models/{TAG}", "cv_run_id": M5_RUN, **fm},
        "method": {
            "inputs_refit": ("m5_avg (average mode) persists its training data only; each fold model's three inputs "
                             "(m1_E_resB0, m2n_E_paralogs, g3_m3n) are refitted on that fold's training lines inside "
                             "the fold model (StackM5.fitted_inputs), as in the CV run, G6 and predict.py; configs and "
                             "seeds unchanged."),
            "self_check": ("all own raw blocks (E, C, M, H, meta, gene-level expression / CN / damaging) -> "
                           "predict_from_blocks for every held-out line; must equal the registered OOF "
                           f"cache/runs/{M5_RUN}/pred.parquet within {CV_REPRO_TOL}"),
            "expression_only": "X = {'E': DepMap log2(TPM+1) row}, IDs 'input::<line>', lineage None (= predict.py)",
            "scoring": ("essml.evaluate.evaluate as scripts/run.py (configs/gene_classes.csv dev-wide classes, "
                        "configs/eval_genes.csv, min_obs 40); paired fold differences scripts/g2_report.paired "
                        "(SE = sd(ddof=1)/sqrt(5)) vs the registered per_fold.csv of each comparator"),
            "hotspot_auprc": ("essml.mutdeps.fold_metrics ap_all (AUPRC of effect < -0.5, all held-out lines), "
                              "configs/mutation_dependencies.json metric_params, mean over the 5 folds; registered "
                              "runs recomputed and checked against results/mutation_dependencies.json"),
        },
        "comparators": {"m1_E_resB0": M1_RUN, "m5_avg_with_genomics": M5_RUN},
        "self_check": checks,
        "evaluate_reproduces_registered_per_fold_max_abs_diff": repro,
        "P1_by_class": by_class,
        "hotspot_auprc": hs,
        "default_rule": {"rule": "dP1 variable (m5_avg expression-only - m1_E_resB0) > 0.005 and > 1 SE",
                         "delta": DELTA, "mean": v["mean"], "se": v["se"], "met": beats,
                         "predict_py_default": "m5_avg" if beats else "m1_E_resB0"},
        "runtime_s": round(time.perf_counter() - t0, 1),
    }
    OUT.write_text(json.dumps(res, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o)))
    log(json.dumps({"P1_variable_folds": by_class["variable"]["P1_folds"],
                    "P1_variable_mean": by_class["variable"]["P1_mean"],
                    "vs_m1": by_class["variable"]["delta_P1_vs_m1_E_resB0"],
                    "vs_m5": by_class["variable"]["delta_P1_vs_m5_avg_with_genomics"],
                    "hotspot": {p: {k: round(h[k]["ap_all_mean"], 4) for k in
                                    ("m5_avg_expression_only", "m1_E_resB0", "m5_avg_with_genomics")}
                                for p, h in hs.items()},
                    "default": res["default_rule"]}, indent=1, default=str))
    log(f"wrote {OUT.relative_to(ROOT)}")
    return res


if __name__ == "__main__":
    main()
