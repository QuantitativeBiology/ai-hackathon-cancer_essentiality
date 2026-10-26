"""Single launcher for dev-CV runs (DESIGN.md §4, §9).

    python scripts/run.py --config configs/runs/<name>.json [--threads N] [--folds 0,1]

Config keys:
    name                 run name (part of run_id)
    model                key in essml.baselines.MODELS (B0, B1, B2) or essml.models.MODELS
                         (M1, M1b, R1, T, M2)
    params               fixed constructor params, e.g. {"n_top": 2000}
    grid                 hyper-parameters tuned by INNER grouped CV, e.g. {"kappa": [0, 1, 10]}
    features             feature matrices passed as X (e.g. ["expression"]); [] for none
    reference_baselines  {"B0": {}, "B1": {"grid": {"kappa": [...]}}} computed inline on the same
                         folds for dP1 / dP2a (their runtime is not counted in the run's runtime)
    inner_folds, min_obs, inner_min_obs, threads, seed (default: configs/seeds.json split_seed)
    requires_gate        gate whose results/gates/<gate>.json must have pass=true (default "G0")
    leaky                1 for deliberately leaky diagnostic runs (flagged in the registry)
    cache_reference      true: reuse B0/B1 reference predictions from cache/ref_preds (exact
                         cache keyed by a hash of the training targets, lineages, held-out IDs
                         and the baseline config; reference baselines use no features)
    save_models          true (or CLI --save-models): persist the fitted fold models to
                         cache/models/<run_id>/ (essml.modelstore). Default false: ~100-300 MB per
                         fold, regenerable in ~100 s/fold; nothing in the default path reads them.
                         Needed only for finalists that will predict new lines (set A / set B).
    line_subset          optional raw feature block (cache/features/<block>.parquet, e.g. "P"):
                         restrict the run (training AND evaluation) to dev lines present in it
    fold_scheme          secondary splits (DESIGN §3): "cv" (default: configs/splits.csv cv_fold) or
                         "lineage" (leave-one-lineage-out: outer fold k = k-th dev lineage in sorted
                         order; summary.json records the fold -> lineage map)
    exclude_from_training  dev ModelIDs never used for training (every outer and inner fit; the
                         sentinel raises if one reaches a fit) but still predicted in their own
                         held-out fold (H3 no-EC variant, EC leave-out split)

Protocol per outer fold k (dev lines only; targets come from essml.guard.load_targets):
    training lines = dev folds != k; held-out = fold k.
    hyper-parameter tuning (inner StratifiedGroupKFold by PatientID, stratified by lineage; the
    criterion uses variable genes of the INNER training part), feature selection and fitting see
    TRAINING lines only. Every fit is checked by a
    sentinel that raises if any held-out ID reaches it, and recorded in FIT_LOG (tests).
    Held-out targets are used only by essml.evaluate after all predictions are made, with the
    dev-wide gene classes configs/gene_classes.csv (identical for every model; DESIGN v3 §3) and
    the fixed evaluable genes configs/eval_genes.csv.
Outputs: results/runs/<run_id>/{summary.json, per_fold.csv, per_gene.csv.gz};
    cache/runs/<run_id>/{pred.parquet, per_gene_fold.parquet, per_line.parquet};
    cache/models/<run_id>/fold<k>.pkl ONLY with save_models / --save-models
    (essml.modelstore.predict_new averages them for new lines);
    results/registry.csv: a row written at start (status=running) and replaced at the end
    (ok / partial_folds / failed: ...), with config hash, commit, dirty-tree flag, splits.csv
    sha256, roles of the target rows loaded, number of distinct IDs that reached any fit.
"""
from __future__ import annotations

import json
import os
import sys

THREAD_VARS = ["OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
               "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"]


def _peek_threads(argv):
    """Thread count from --threads or the config, read with the stdlib before numpy is imported."""
    if "--threads" in argv:
        return int(argv[argv.index("--threads") + 1])
    if "--config" in argv:
        try:
            with open(argv[argv.index("--config") + 1]) as fh:
                return int(json.load(fh).get("threads", 12))
        except (OSError, ValueError):
            pass
    return 12


if __name__ == "__main__":
    _T = _peek_threads(sys.argv)
    for _v in THREAD_VARS:
        os.environ[_v] = str(_T)

import argparse  # noqa: E402
import csv  # noqa: E402
import datetime as dt  # noqa: E402
import hashlib  # noqa: E402
import subprocess  # noqa: E402
import time  # noqa: E402
import traceback  # noqa: E402
import warnings  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from sklearn.model_selection import StratifiedGroupKFold  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from essml import baselines, genesets, guard, io, metrics, modelstore  # noqa: E402
from essml import models as _models  # noqa: E402
from essml.modelstore import predict_anchor, predict_new  # noqa: E402,F401  (re-exported)
from essml.evaluate import evaluate  # noqa: E402

REGISTRY = ROOT / "results" / "registry.csv"
from essml.models.g3 import MODELS as _G3_MODELS  # noqa: E402  (G3 nonlinear track)
MODELS = {**baselines.MODELS, **_models.MODELS, **_G3_MODELS}
REF_CACHE = ROOT / "cache" / "ref_preds"
REG_COLS = ["run_id", "timestamp", "name", "model", "config", "config_hash", "git_commit",
            "dirty_tree", "splits_sha256", "roles_used", "n_train_ids", "leaky", "status",
            "n_lines", "n_genes_eval", "chosen_params",
            "P1_variable", "P1_variable_sd", "P1_all", "P1_all_sd", "P1_selective", "P1_selective_sd",
            "P1_rare_selective", "P1_median_variable", "n_r_gt_0.3_variable",
            "dP1_vs_B1_variable", "dP1_vs_B1_variable_sd", "frac_beat_B1_variable",
            "P2a", "P2a_sd", "dP2a_vs_B0", "dP2a_vs_B0_sd",
            "P2b_nap_prev05", "P2b_nap_prev05_sd", "P2b_auroc_prev05", "P2b_auroc_prev05_sd",
            "fit_s", "predict_s", "total_s", "threads", "peak_rss_mb", "rss_method",
            # provenance columns added after G2 (older rows: empty)
            "eval_genes_sha256", "gene_classes_sha256", "tree_clean"]

# ----------------------------------------------------------------------------------------------
# Leakage instrumentation
FIT_LOG: list[dict] = []  # one entry per fit call: {"tag", "outer_fold", "ids"}


class LeakageError(AssertionError):
    pass


class Sentinel:
    """Holds the held-out IDs of the current outer fold; any fit receiving one raises."""

    def __init__(self, forbidden, outer_fold):
        self.forbidden = set(map(str, forbidden))
        self.outer_fold = outer_fold

    def check(self, tag, *frames):
        for fr in frames:
            if fr is None:
                continue
            bad = self.forbidden.intersection(map(str, fr.index))
            if bad:
                raise LeakageError(f"{tag}: {len(bad)} held-out IDs passed to fit, e.g. {sorted(bad)[:3]}")


def guarded_fit(model, sentinel: Sentinel, tag, Y, lineage, X, groups=None):
    """The only place models are fitted: sentinel check, FIT_LOG record, then model.fit.

    Models with ``uses_groups = True`` (M1/M2: leave-patient-out LOO, grouped inner folds) also
    receive the training lines' PatientIDs."""
    lin = None if lineage is None else lineage.loc[Y.index]
    Xf = None if X is None else X.loc[Y.index]
    sentinel.check(tag, Y, None if lin is None else lin.to_frame(), Xf)
    FIT_LOG.append({"tag": tag, "outer_fold": sentinel.outer_fold, "ids": list(map(str, Y.index))})
    if getattr(model, "uses_groups", False) and groups is not None:
        return model.fit(Y, lineage=lin, X=Xf, groups=groups.loc[Y.index].to_numpy())
    return model.fit(Y, lineage=lin, X=Xf)


# ----------------------------------------------------------------------------------------------
# Inner CV tuning (training lines only)
def _grid_points(grid):
    keys = sorted(grid)
    pts = [{}]
    for k in keys:
        pts = [dict(p, **{k: v}) for p in pts for v in grid[k]]
    return pts


def _predict(model, ids, lineage, X, **kw):
    lin = None if lineage is None else lineage.loc[ids]
    Xp = None if X is None else X.loc[ids]
    return model.predict(ids, lineage=lin, X=Xp, **kw)


def inner_splits(ids, groups, strata, n_splits, seed):
    sgk = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)  # small lineages < n_splits
        return [(ids[a], ids[b]) for a, b in sgk.split(np.zeros(len(ids)), strata.loc[ids].to_numpy(),
                                                        groups.loc[ids].to_numpy())]


def tune(model_key, params, grid, Y_tr, lineage, X, groups, sentinel, eval_genes, neg, n_splits,
         seed, inner_min_obs):
    """Choose the grid point maximising inner-CV P1 (variable genes of the INNER training part,
    within each inner held-out fold, averaged). Returns (best params, table of scores)."""
    cls_ = MODELS[model_key]
    pts = _grid_points(grid)
    if len(pts) == 1 and not pts[0]:
        return {}, []
    predict_keys = set(getattr(cls_, "predict_params", ()))
    scores = np.zeros((len(pts), n_splits))
    for j, (itr, iva) in enumerate(inner_splits(Y_tr.index, groups, lineage, n_splits, seed)):
        Yi = Y_tr.loc[itr]
        var = genesets.gene_classes(Yi, neg)["variable"]
        g = var.index[var.to_numpy()].intersection(eval_genes)
        Tva = Y_tr.loc[iva, g].to_numpy()  # inner held-out targets: part of the TRAINING fold
        fitted = {}
        for i, p in enumerate(pts):
            fit_p = {k: v for k, v in p.items() if k not in predict_keys}
            key = json.dumps(fit_p, sort_keys=True)
            if key not in fitted:
                fitted[key] = guarded_fit(cls_(**params, **fit_p), sentinel, f"inner{j}", Yi,
                                          lineage, X, groups)
            pred_p = {k: v for k, v in p.items() if k in predict_keys}
            Pva = _predict(fitted[key], iva, lineage, X, **pred_p).loc[:, g].to_numpy()
            scores[i, j] = np.nanmean(metrics.per_gene_pearson(Pva, Tva, inner_min_obs))
    mean = scores.mean(1)
    best = int(np.argmax(mean))  # first maximum = earliest grid point
    table = [{"params": p, "inner_P1": float(m), "inner_folds": scores[i].round(5).tolist()}
             for i, (p, m) in enumerate(zip(pts, mean))]
    return pts[best], table


def fit_predict_fold(model_key, params, grid, Y_tr, ho_ids, lineage, X, groups, sentinel,
                     eval_genes, neg, inner_folds, seed, inner_min_obs):
    """Tune on training lines, refit on all training lines, predict held-out lines."""
    t0 = time.perf_counter()
    best, table = tune(model_key, params, grid, Y_tr, lineage, X, groups, sentinel, eval_genes,
                       neg, inner_folds, seed, inner_min_obs)
    cls_ = MODELS[model_key]
    pk = set(getattr(cls_, "predict_params", ()))
    model = guarded_fit(cls_(**params, **best), sentinel, "outer", Y_tr, lineage, X, groups)
    t_fit = time.perf_counter() - t0
    t0 = time.perf_counter()
    pred = _predict(model, ho_ids, lineage, X, **{k: v for k, v in best.items() if k in pk})
    t_pred = time.perf_counter() - t0
    for k, v in best.items():  # make the stored model predict with the chosen params
        setattr(model, k, v)
    return model, pred, best, table, t_fit, t_pred


def run_cv(Y, folds, groups, lineage, X, cfg, eval_genes, neg, seed, folds_to_run=None,
           log=print):
    """Outer grouped dev CV. Y: dev targets (DataFrame). folds: Series ModelID -> outer fold.

    Returns dict with OOF predictions, reference baseline predictions,
    fitted models, chosen params, inner tables and timings. Held-out targets are never read here.
    """
    refs = cfg.get("reference_baselines", {})
    fl = [int(f) for f in (sorted(folds.unique()) if folds_to_run is None else folds_to_run)]
    preds, ref_preds = [], {b: [] for b in refs}
    out = {"models": {}, "chosen": {}, "inner": {}, "fit_s": {}, "predict_s": {}}
    excl = pd.Index(map(str, cfg.get("exclude_from_training", [])))
    for f in fl:
        ho = folds.index[folds == f]
        tr = folds.index[folds != f]
        tr = tr[~tr.astype(str).isin(excl)]  # excluded lines never train (still predicted)
        Y_tr = Y.loc[tr]  # training targets only from here on
        sent = Sentinel(ho.union(excl), f)
        model, pred, best, table, tf, tp = fit_predict_fold(
            cfg["model"], cfg.get("params", {}), cfg.get("grid", {}), Y_tr, ho, lineage, X, groups,
            sent, eval_genes, neg, cfg.get("inner_folds", 5), seed + int(f),
            cfg.get("inner_min_obs", 30))
        preds.append(pred)
        out["models"][f], out["chosen"][f], out["inner"][f] = model, best, table
        out["fit_s"][f], out["predict_s"][f] = tf, tp
        log(f"  fold {f}: {best} fit {tf:.1f}s predict {tp:.2f}s")
        if hasattr(model, "report_"):
            out.setdefault("model_report", {})[f] = model.report_
            log(f"    {json.dumps(model.report_, default=str)[:600]}")
        for b, bcfg in refs.items():
            key = None
            if cfg.get("cache_reference"):
                # exact cache: deterministic given (baseline config, inner settings, seed, training
                # targets and lineages, held-out IDs); the reference sees no features (B0/B1 only)
                h = hashlib.sha256(json.dumps([b, bcfg, cfg.get("inner_folds", 5),
                                               cfg.get("inner_min_obs", 30), seed + int(f),
                                               list(map(str, ho))], sort_keys=True).encode())
                h.update(pd.util.hash_pandas_object(Y_tr, index=True).to_numpy().tobytes())
                h.update(Y_tr.columns.str.cat(sep=",").encode())
                h.update(lineage.loc[tr].str.cat(sep=",").encode())
                key = REF_CACHE / f"{b}_{h.hexdigest()[:24]}.pkl"
            if key is not None and key.exists() and b in ("B0", "B1"):
                bp, bb = pd.read_pickle(key)
            else:
                _, bp, bb, _, _, _ = fit_predict_fold(
                    b, bcfg.get("params", {}), bcfg.get("grid", {}), Y_tr, ho, lineage, X, groups,
                    sent, eval_genes, neg, cfg.get("inner_folds", 5), seed + int(f),
                    cfg.get("inner_min_obs", 30))
                if key is not None and b in ("B0", "B1"):
                    key.parent.mkdir(parents=True, exist_ok=True)
                    pd.to_pickle((bp, bb), key)
            ref_preds[b].append(bp)
            out["chosen"].setdefault(f"ref_{b}", {})[f] = bb
    out["pred"] = pd.concat(preds)
    out["ref_pred"] = {b: pd.concat(v) for b, v in ref_preds.items()}
    return out


# ----------------------------------------------------------------------------------------------
# Provenance helpers
def config_hash(cfg):
    return hashlib.sha256(json.dumps(cfg, sort_keys=True).encode()).hexdigest()


def git_commit():
    """(HEAD commit, dirty flag 0/1 for modified tracked files)."""
    try:
        h = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True,
                           text=True).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain", "--untracked-files=no", "--", "."]
                               + [f":(exclude){p}" for p in OUTPUT_PATHS], cwd=ROOT,
                               capture_output=True, text=True).stdout.strip()
        return h, int(bool(dirty))
    except OSError:
        return "unknown", 1


# paths written by the run itself; they must not make the tree "unclean" for the next launch
OUTPUT_PATHS = ("results/registry.csv", "results/runs")


def tree_clean():
    """1 if the work tree is clean, else 0: no modified or untracked (non-ignored) files, ignoring only
    the run's own outputs (OUTPUT_PATHS). Stricter than ``dirty_tree``, which ignores untracked files."""
    try:
        r = subprocess.run(["git", "status", "--porcelain", "--", "."]
                           + [f":(exclude){p}" for p in OUTPUT_PATHS],
                           cwd=ROOT, capture_output=True, text=True)
        return int(r.returncode == 0 and not r.stdout.strip())
    except OSError:
        return 0


def file_sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require_gate(gate, gates_dir=ROOT / "results" / "gates"):
    """Refuse to start unless results/gates/<gate>.json exists with pass == true (DESIGN §9)."""
    p = Path(gates_dir) / f"{gate}.json"
    if not p.exists() or not json.loads(p.read_text()).get("pass"):
        raise RuntimeError(f"gate {gate} has not passed ({p} missing or pass != true)")


def peak_rss_mb():
    try:
        import psutil
        mi = psutil.Process().memory_info()
        if hasattr(mi, "peak_wset"):  # Windows: true peak working set
            return mi.peak_wset / 2 ** 20, "psutil.peak_wset"
        try:
            import resource
            return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, "ru_maxrss"
        except ImportError:
            return mi.rss / 2 ** 20, "psutil.rss_at_end(not peak)"
    except ImportError:
        return np.nan, "psutil unavailable"


def write_registry(row, path=REGISTRY):
    """Insert the row, or replace the existing row with the same run_id (start -> end update)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lock = path.with_suffix(".lock")
    for _ in range(600):  # concurrent launchers: mkdir is atomic
        try:
            os.mkdir(lock)
            break
        except FileExistsError:
            time.sleep(0.1)
    try:
        _write_registry_locked(row, path)
    finally:
        try:
            os.rmdir(lock)
        except OSError:
            pass


def _write_registry_locked(row, path):
    rows = []
    if path.exists():
        with open(path, newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
    new = {k: row.get(k, "") for k in REG_COLS}
    rows = [r for r in rows if r.get("run_id") != row["run_id"]] + [new]
    tmp = path.with_suffix(".tmp")
    with open(tmp, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=REG_COLS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in REG_COLS})
    os.replace(tmp, path)


def _fmt(x):
    return round(x, 5) if isinstance(x, float) and np.isfinite(x) else x


def load_eval_genes(path=ROOT / "configs" / "eval_genes.csv"):
    return pd.Index(pd.read_csv(path, dtype={"entrez": str})["entrez"])


def per_gene_aggregate(pg: pd.DataFrame) -> pd.DataFrame:
    """Gene-level summary across folds (committed); the long per-fold table goes to cache/."""
    num = ["r", "rho", "ap", "nap", "auroc", "prev"] + [c for c in pg.columns if c.startswith("r_")]
    g = pg.groupby("gene", sort=False)
    out = g[num].mean().add_suffix("_mean")
    out.insert(0, "r_sd", g["r"].std(ddof=1))
    out.insert(0, "n_folds", g["r"].count())
    for c in [c for c in pg.columns if c.startswith("is_")]:
        out["frac_folds_" + c[3:]] = g[c].mean()
    return out.round(4).reset_index()


# ----------------------------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", required=True)
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--save-models", action="store_true",
                    help="persist fitted fold models to cache/models (default off; or cfg save_models)")
    ap.add_argument("--folds", default=None, help="comma-separated outer folds (debug; not headline)")
    a = ap.parse_args(argv)
    cfg = json.loads(Path(a.config).read_text())
    threads = a.threads or int(cfg.get("threads", 12))
    try:
        from threadpoolctl import threadpool_limits
        threadpool_limits(threads)
    except ImportError:
        pass
    seed = int(cfg.get("seed", json.loads((ROOT / "configs" / "seeds.json").read_text())["split_seed"]))
    ch = config_hash(cfg)
    stamp = dt.datetime.now()
    run_id = f"{stamp:%Y%m%d-%H%M%S}_{cfg['name']}_{ch[:8]}"
    commit, dirty = git_commit()
    splits_path = ROOT / "configs" / "splits.csv"
    row = {"run_id": run_id, "timestamp": stamp.isoformat(timespec="seconds"), "name": cfg["name"],
           "model": cfg["model"], "config": str(Path(a.config).as_posix()), "config_hash": ch,
           "git_commit": commit, "dirty_tree": dirty, "splits_sha256": file_sha256(splits_path),
           "eval_genes_sha256": file_sha256(ROOT / "configs" / "eval_genes.csv"),
           "gene_classes_sha256": file_sha256(ROOT / "configs" / "gene_classes.csv"),
           "tree_clean": tree_clean(), "leaky": int(bool(cfg.get("leaky", 0))), "threads": threads, "status": "running"}
    t_start = time.perf_counter()
    write_registry(row)  # row at start; replaced at the end (also on failure)
    try:
        require_gate(cfg.get("requires_gate", "G0"))
        splits = pd.read_csv(splits_path, index_col=0)
        dev = splits[splits["role"] == "dev"]
        if cfg.get("line_subset"):
            keep = pd.read_parquet(ROOT / "cache" / "features" / f"{cfg['line_subset']}.parquet",
                                   columns=[]).index
            dev = dev[dev.index.isin(keep)]
        folds = dev["cv_fold"].astype(int)
        lineage = dev["OncotreeLineage"].astype(str)
        scheme = cfg.get("fold_scheme", "cv")
        fold_names = None
        if scheme == "lineage":  # leave-one-lineage-out (DESIGN §3 secondary split)
            levels = sorted(lineage.unique())
            folds = lineage.map({l: i for i, l in enumerate(levels)}).astype(int)
            fold_names = dict(enumerate(levels))
        elif scheme != "cv":
            raise ValueError(f"unknown fold_scheme {scheme}")
        excl = pd.Index(map(str, cfg.get("exclude_from_training", [])))
        if not excl.isin(dev.index).all():
            raise ValueError("exclude_from_training must list dev lines of this run")
        groups = dev["PatientID"].astype(str)
        Y = guard.load_targets(dev.index, config_hash=ch).astype(np.float64)
        roles = guard.roles_of(Y.index).value_counts()
        row["roles_used"] = ";".join(f"{k}:{v}" for k, v in sorted(roles.items()))
        eval_genes = load_eval_genes()
        gene_cls = genesets.load_gene_classes()  # dev-wide, committed (evaluation reporting)
        X = None
        feats = cfg.get("features", [])
        if feats:
            mats = [io.load_matrix(fn).reindex(dev.index).add_prefix(f"{fn}:") for fn in feats]
            X = pd.concat(mats, axis=1)
            if X.isna().any().any():
                raise ValueError("features have NaN for some dev lines")
        neg, neg_rep = genesets.load_hart_neg()
        print(f"run {run_id}: {cfg['model']} dev lines {len(Y)}, eval genes {len(eval_genes)}, "
              f"NEGv1 mapped {neg_rep['n_mapped']}/{neg_rep['n_neg']}")
        fl = None if a.folds is None else [int(x) for x in a.folds.split(",")]
        n_log = len(FIT_LOG)
        res = run_cv(Y, folds, groups, lineage, X, cfg, eval_genes, neg, seed, fl)
        if cfg.get("oracle_select"):  # deliberately leaky diagnostic (e.g. non-nested R4)
            if not cfg.get("leaky"):
                raise LeakageError("oracle_select reads held-out targets: the config must set leaky=1")
            res["pred"] = MODELS[cfg["model"]].oracle_select(res, Y, folds)
        train_ids = set().union(*[set(e["ids"]) for e in FIT_LOG[n_log:]])
        if (guard.roles_of(sorted(train_ids)) != "dev").any():
            raise LeakageError("non-dev IDs reached a fit")
        row["n_train_ids"] = len(train_ids)

        ho_ids = res["pred"].index
        if fl is not None:
            row["status"] = "partial_folds"
        ev = evaluate(res["pred"], Y.loc[ho_ids], gene_cls, folds.loc[ho_ids],
                      baseline_preds={b: p for b, p in res["ref_pred"].items()},
                      genes=eval_genes, min_obs=int(cfg.get("min_obs", 40)))
        hl = ev["headline"]

        out_dir = ROOT / "results" / "runs" / run_id
        cache_dir = ROOT / "cache" / "runs" / run_id
        out_dir.mkdir(parents=True, exist_ok=True)
        cache_dir.mkdir(parents=True, exist_ok=True)
        res["pred"].astype(np.float32).to_parquet(cache_dir / "pred.parquet")
        ev["per_gene"].to_parquet(cache_dir / "per_gene_fold.parquet")
        ev["per_line"].to_parquet(cache_dir / "per_line.parquet")
        ev["per_fold"].round(6).to_csv(out_dir / "per_fold.csv", index=False)
        per_gene_aggregate(ev["per_gene"]).to_csv(out_dir / "per_gene.csv.gz", index=False)
        if a.save_models or cfg.get("save_models"):
            for f, m in res["models"].items():
                if hasattr(m, "compact"):
                    m.compact()  # drop re-buildable arrays before pickling
                modelstore.save_fold_model(run_id, f, m, folds.index[folds != f],
                                           [] if X is None else X.columns)
        fit_s, pred_s = sum(res["fit_s"].values()), sum(res["predict_s"].values())
        rss, rss_m = peak_rss_mb()
        summary = {"run_id": run_id, "config": cfg, "config_hash": ch, "git_commit": commit,
                   "dirty_tree": dirty, "splits_sha256": row["splits_sha256"],
                   "eval_genes_sha256": row["eval_genes_sha256"],
                   "gene_classes_sha256": row["gene_classes_sha256"], "tree_clean": row["tree_clean"],
                   "roles_used": row["roles_used"], "n_train_ids": row["n_train_ids"],
                   "threads": threads, "seed": seed, "n_lines": len(ho_ids),
                   "fold_scheme": scheme, "fold_names": fold_names,
                   "excluded_from_training": list(excl),
                   "n_eval_genes": len(eval_genes), "neg_mapping": neg_rep,
                   "chosen_params": {str(k): v for k, v in res["chosen"].items()},
                   "inner_cv": {str(k): v for k, v in res["inner"].items()},
                   "model_report": {str(k): v for k, v in res.get("model_report", {}).items()},
                   "fit_s_per_fold": res["fit_s"], "predict_s_per_fold": res["predict_s"],
                   "fit_s": fit_s, "predict_s": pred_s, "peak_rss_mb": rss, "rss_method": rss_m,
                   "headline": hl, "summary": ev["summary"]}
        summary["fit_s_per_fold"] = {str(k): v for k, v in res["fit_s"].items()}
        summary["predict_s_per_fold"] = {str(k): v for k, v in res["predict_s"].items()}
        (out_dir / "summary.json").write_text(json.dumps(summary, indent=1, default=float))
        row.update({k: _fmt(v) for k, v in hl.items() if k in REG_COLS})
        row.update({"status": row["status"] if row["status"] != "running" else "ok",
                    "n_lines": len(ho_ids), "n_genes_eval": len(eval_genes),
                    "chosen_params": json.dumps({str(k): v for k, v in res["chosen"].items()
                                                 if not str(k).startswith("ref_")}),
                    "fit_s": round(fit_s, 2), "predict_s": round(pred_s, 2),
                    "total_s": round(time.perf_counter() - t_start, 1),
                    "peak_rss_mb": round(rss, 1), "rss_method": rss_m})
        write_registry(row)
        print(json.dumps({k: _fmt(v) for k, v in hl.items() if not k.endswith("_sd")}, indent=1))
        print(f"saved results/runs/{run_id}")
        return run_id
    except BaseException as e:
        row.update({"status": f"failed: {type(e).__name__}: {str(e)[:200]}",
                    "total_s": round(time.perf_counter() - t_start, 1)})
        write_registry(row)
        traceback.print_exc()
        raise


if __name__ == "__main__":
    main()
