"""Shared plumbing for the one-time held-out evaluations: scripts/evaluate_test.py (G4) and
scripts/evaluate_transfer.py (G6). No metric code lives here (frozen: essml.evaluate / metrics /
transfer; resampling: essml.heldout).

* ``Paths``: every file location, injectable so tests run on toy repositories (tests/test_guard.py
  pattern) and never touch real data.
* ``load_finalists``: parses ``configs/finalists.json`` (schema: configs/finalists.example.json) and
  checks that each config file hashes to its recorded ``config_hash`` (scripts/run.py convention).
* Tree checks: ``strict_dirty`` (modified OR untracked non-ignored files) and ``commit_access_log``.
  ``results/test_access_log.csv`` is TRACKED, and the guard appends a row on every non-dev access,
  so after the first access the tree is dirty and the guard would refuse every later finalist
  ("dirty tree"). The scripts therefore commit the access log (that path only, after verifying it
  is the only modified tracked file) immediately BEFORE each further access and once at the end
  (``ensure_log_committed``, with retries). The access itself is recorded in the ledger and
  evaluated first, so a failing commit can only stop the NEXT access, never waste the current one.
  If the PI untracks the log instead, nothing is committed.
* ``AccessLedger``: in-process record of every target access; refuses a second request for the same
  (config_hash, scope) before the guard is even called, and is written next to the results.
* Model fitting reuses scripts/run.py (``fit_predict_fold``, ``run_cv``, ``Sentinel``, FIT_LOG):
  the same tuning protocol, seeds and leakage sentinel as dev CV.
"""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import os
import subprocess
import sys
import time
from dataclasses import dataclass, field, replace
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
THREAD_VARS = ["OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS",
               "VECLIB_MAXIMUM_THREADS"]
EC_DEV = ["ACH-001578", "ACH-002291", "ACH-001374", "ACH-001668"]  # D20d: the 4 EC dev lines
ROLE_ORDER = ["headline", "T", "H3_noEC", "B1", "B2", "B0"]  # access order (most important first)
ADDITIONAL_PREFIX = "add_"  # role of an additional finalist (D23f: m1_E_resB0) = add_<name>


def access_order(fin: dict) -> list:
    """headline, additional finalists (file order), T, H3 no-EC, B1, B2, B0."""
    adds = [r for r in fin if r.startswith(ADDITIONAL_PREFIX)]
    return ["headline"] + adds + [r for r in ROLE_ORDER if r != "headline" and r in fin]


def set_threads_early(argv, default=4):
    """Thread env vars from --threads (stdlib only; call before numpy is imported)."""
    n = default
    if "--threads" in argv:
        n = int(argv[argv.index("--threads") + 1])
    for v in THREAD_VARS:
        os.environ[v] = str(n)
    return n


import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
from essml import genesets, guard, io, modelstore  # noqa: E402
import run as RUN  # noqa: E402  (scripts/run.py: launcher protocol, sentinel, registry)


# ============================================================================ paths
@dataclass
class Paths:
    root: Path = ROOT
    splits: Path = ROOT / "configs" / "splits.csv"
    crispr: Path | None = None    # None = the guard's own target file (only toy tests substitute one)
    access_log: Path = guard.ACCESS_LOG
    eval_genes: Path = ROOT / "configs" / "eval_genes.csv"
    gene_classes: Path = ROOT / "configs" / "gene_classes.csv"
    control_split: Path = ROOT / "configs" / "control_split.json"
    library_groups: Path = ROOT / "configs" / "library_groups.csv"
    seeds: Path = ROOT / "configs" / "seeds.json"
    gates: Path = ROOT / "results" / "gates"
    registry: Path = ROOT / "results" / "registry.csv"
    model_dir: Path = ROOT / "cache" / "models"
    frozen_hashes: Path = ROOT / "configs" / "frozen_hashes.json"
    criteria: Path = ROOT / "configs" / "transfer_criteria.json"
    harmonised_dir: Path = guard.HARMONISED_DIR
    comparators: Path = ROOT / "cache" / "transfer" / "comparators.parquet"
    paralog_pairs: Path = ROOT / "cache" / "priors" / "paralog_pairs.parquet"
    external_expr: Path = ROOT / "cache" / "external_expr" / "harmonised.parquet"
    published_sets: Path = ROOT / "configs" / "published_genesets.csv"
    predictable_genes: Path = ROOT / "configs" / "predictable_genes.csv"  # D22 tiers + confidence (DESIGN §4)
    cv_cache: Path = ROOT / "cache" / "runs"
    out_dir: Path = ROOT / "results" / "test"
    pred_dir: Path = ROOT / "cache" / "test" / "preds"
    extra: dict = field(default_factory=dict)

    def guard_kw(self):
        kw = dict(splits_path=self.splits, log_path=self.access_log, repo_root=self.root)
        if self.crispr is not None:
            kw["crispr_path"] = self.crispr
        return kw

    def rel(self, p) -> str:
        try:
            return Path(p).resolve().relative_to(Path(self.root).resolve()).as_posix()
        except ValueError:
            return Path(p).as_posix()

    def with_(self, **kw):
        return replace(self, **kw)


# ============================================================================ finalists
@dataclass
class Finalist:
    role: str          # headline | T | H3_noEC | B0 | B1 | B2
    name: str
    config_path: str
    config_hash: str
    cfg: dict
    cv_run_id: str = ""
    own_predictable_genes: str = ""   # optional secondary tier table of this finalist (D22f: m1_E_resB0)

    @property
    def tag(self):
        return f"{self.name}_{self.config_hash[:8]}"


def config_hash(cfg) -> str:
    return RUN.config_hash(cfg)


def load_finalists(path, root=ROOT) -> dict:
    """{role: Finalist}; raises on schema or hash problems (see configs/finalists.example.json)."""
    d = json.loads(Path(path).read_text(encoding="utf-8"))
    raw = {"headline": d.get("headline"), "T": d.get("T"), "H3_noEC": d.get("H3_noEC")}
    for b in ("B0", "B1", "B2"):
        raw[b] = (d.get("baselines") or {}).get(b)
    adds = d.get("additional") or []
    if not isinstance(adds, list):
        raise ValueError("'additional' must be a list of finalist entries")
    for e in adds:
        raw[ADDITIONAL_PREFIX + str(e.get("name"))] = e
    missing = [k for k, v in raw.items() if not v]
    if missing:
        raise ValueError(f"finalists file lacks {missing}")
    out = {}
    for role, e in raw.items():
        for k in ("name", "config", "config_hash"):
            if not e.get(k):
                raise ValueError(f"finalist {role}: missing '{k}'")
        cfg = json.loads((Path(root) / e["config"]).read_text(encoding="utf-8"))
        h = config_hash(cfg)
        if h != e["config_hash"]:
            raise ValueError(f"finalist {role} ({e['name']}): config hash {h} != recorded {e['config_hash']}")
        if cfg.get("name") != e["name"]:
            raise ValueError(f"finalist {role}: name {e['name']} != config name {cfg.get('name')}")
        out[role] = Finalist(role, e["name"], e["config"], h, cfg, e.get("cv_run_id", ""),
                             e.get("own_predictable_genes", ""))
    if len({f.config_hash for f in out.values()}) != len(out):
        raise ValueError("finalists must have distinct config hashes (one access each)")
    for b in ("B0", "B1", "B2"):
        if out[b].cfg["model"] != b:
            raise ValueError(f"baseline {b} has model {out[b].cfg['model']}")
    if out["T"].cfg["model"] != "T" or out["H3_noEC"].cfg["model"] != "T":
        raise ValueError("T and H3_noEC must be model T (DESIGN §5)")
    return out


def check_noec(f: Finalist, ec_lines=EC_DEV):
    """D20d erratum: the H3 no-EC variant is T trained without exactly the 4 EC dev lines."""
    excl = sorted(map(str, f.cfg.get("exclude_from_training", [])))
    if excl != sorted(ec_lines):
        raise ValueError(f"H3 no-EC variant must exclude exactly the EC dev lines {sorted(ec_lines)} (D20d); "
                         f"config excludes {excl}")
    return {k: v for k, v in f.cfg.items() if k not in NOEC_FREE_KEYS}


NOEC_FREE_KEYS = ("name", "exclude_from_training", "expected", "threads", "requires_gate")


def check_noec_matches_T(noec: Finalist, t: Finalist):
    """The no-EC variant must be T in every fitting-relevant key (only the exclusion differs)."""
    a = check_noec(noec)
    b = {k: v for k, v in t.cfg.items() if k not in NOEC_FREE_KEYS}
    if a != b:
        diff = sorted(k for k in set(a) | set(b) if a.get(k) != b.get(k))
        raise ValueError(f"H3 no-EC config differs from T beyond the exclusion: {diff}")
    return True


# ============================================================================ git / tree
def _git(args, root):
    return subprocess.run(["git"] + list(args), cwd=root, capture_output=True, text=True)


def is_tracked(root, rel) -> bool:
    return (Path(root) / rel).exists() and _git(["ls-files", "--error-unmatch", rel], root).returncode == 0


def strict_dirty(root, allow_untracked=()) -> list:
    """Porcelain lines for modified tracked files and untracked non-ignored files, except untracked
    files under ``allow_untracked`` prefixes (the evaluation's own outputs, on --resume)."""
    r = _git(["status", "--porcelain", "--untracked-files=all"], root)
    if r.returncode != 0:
        return [f"git status failed: {r.stderr.strip()}"]
    bad = []
    for line in r.stdout.splitlines():
        if not line.strip():
            continue
        path = line[3:].strip().strip('"')
        if line.startswith("??") and any(path.startswith(a) for a in allow_untracked):
            continue
        bad.append(line)
    return bad


def git_head(root) -> str:
    return _git(["rev-parse", "HEAD"], root).stdout.strip() or "unknown"


def commit_access_log(paths: Paths, message: str) -> str | None:
    """Commit results/test_access_log.csv alone (see module doc). Returns the commit hash, or None if
    the log is untracked (nothing to keep clean). Raises if another tracked file is modified."""
    rel = paths.rel(paths.access_log)
    if not is_tracked(paths.root, rel):
        return None
    st = _git(["status", "--porcelain", "--untracked-files=no"], paths.root).stdout.splitlines()
    if not any(l.strip() and l[3:].strip() == rel for l in st):
        return git_head(paths.root)  # the log is unchanged: nothing to commit
    other = [l for l in st if l.strip() and l[3:].strip() != rel]
    if other:
        raise RuntimeError(f"tracked files other than the access log are modified: {other}")
    if not st:
        return git_head(paths.root)
    for args in (["add", "--", rel], ["commit", "-q", "-m", message, "--", rel]):
        r = _git(args, paths.root)
        if r.returncode != 0:
            raise RuntimeError(f"git {' '.join(args)} failed: {r.stderr.strip()}")
    return git_head(paths.root)


def commit_new_file(paths: Paths, path, message: str) -> str:
    """Commit one evaluation output (e.g. the predictions manifest, before any access: tamper-evident)."""
    rel = paths.rel(path)
    if not _git(["status", "--porcelain", "--", rel], paths.root).stdout.strip():
        return git_head(paths.root)  # already committed and unchanged
    for args in (["add", "--", rel], ["commit", "-q", "-m", message, "--", rel]):
        r = _git(args, paths.root)
        if r.returncode != 0:
            raise RuntimeError(f"git {' '.join(args)} failed: {r.stderr.strip()}")
    return git_head(paths.root)


def ensure_log_committed(paths: Paths, message: str, retries=5, wait=2.0):
    """Commit the access log if it is modified (before the next access / at the end); retried, then fatal."""
    last = None
    for i in range(retries):
        try:
            return commit_access_log(paths, message)
        except RuntimeError as e:  # e.g. index.lock held by an editor, or another tracked file modified
            last = e
            time.sleep(wait * (i + 1))
    raise RuntimeError(f"access log could not be committed ({last}); stopping BEFORE the next access. Commit "
                       f"{paths.rel(paths.access_log)} alone, then re-run with --resume.")


def sha256_lf(path) -> str:
    return hashlib.sha256(Path(path).read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def sha256_file(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def frozen_ok(paths: Paths) -> list:
    """Frozen evaluation modules whose sha256 differs from configs/frozen_hashes.json."""
    ref = json.loads(Path(paths.frozen_hashes).read_text())["files"]
    return [f for f, h in ref.items() if sha256_lf(Path(paths.root) / f) != h]


def gate_passed(paths: Paths, gate) -> bool:
    p = Path(paths.gates) / f"{gate}.json"
    return p.exists() and bool(json.loads(p.read_text()).get("pass"))


# ============================================================================ access ledger
class AccessLedger:
    """Every target access of this process; a repeated (config_hash, scope) is refused here."""

    def __init__(self, path=None):
        self.rows = []
        self.path = path

    def check(self, config_hash, scope):
        if any(r["config_hash"] == config_hash and r["scope"] == scope and r["granted"] for r in self.rows):
            raise guard.BlindingError(f"second access for {config_hash[:8]} scope {scope} refused by the ledger")

    def record(self, finalist, scope, n_ids, granted, note="", commit=None):
        self.rows.append({"timestamp": dt.datetime.now().isoformat(timespec="seconds"), "finalist": finalist.name,
                          "role": finalist.role, "config_hash": finalist.config_hash, "scope": scope,
                          "n_ids": int(n_ids), "granted": bool(granted), "note": note, "log_commit": commit})
        if self.path:
            Path(self.path).parent.mkdir(parents=True, exist_ok=True)
            Path(self.path).write_text(json.dumps(self.rows, indent=1))

    def n_granted(self, config_hash=None, scope=None):
        return sum(1 for r in self.rows if r["granted"] and (config_hash is None or r["config_hash"] == config_hash)
                   and (scope is None or r["scope"] == scope))


def guarded_targets(finalist: Finalist, ids, paths: Paths, ledger: AccessLedger, dry: bool, screen=None,
                    scope_label=None, commit_log=True) -> pd.DataFrame:
    """The one target request of ``finalist`` for ``ids`` (scope test / setA:<screen>).

    dry=True: every ID must be a dev line (checked BEFORE the guard is called), so the guard never
    sees a non-dev ID. Otherwise: after a granted access the access log is committed (module doc).
    """
    roles = guard.roles_of(ids, paths.splits)
    if dry and (roles != "dev").any():
        raise guard.BlindingError(f"dry run requested non-dev IDs: {roles[roles != 'dev'].value_counts().to_dict()}")
    scope = scope_label or ("dev-standin" if dry else (guard._scope(roles, screen)))
    ledger.check(finalist.config_hash, scope)
    if not dry and commit_log:  # earlier rows of the tracked log are committed before this access
        ensure_log_committed(paths, f"G4/G6 access log (before {finalist.name} {scope})")
    try:
        Y = guard.load_targets(list(ids), config_hash=finalist.config_hash, screen=screen, **paths.guard_kw())
    except guard.BlindingError as e:
        ledger.record(finalist, scope, len(ids), False, note=str(e)[:300])
        raise
    ledger.record(finalist, scope, len(ids), True)
    return Y


# ============================================================================ populations
@dataclass
class Population:
    train_ids: pd.Index         # lines the refit / fold models may train on (before exclusions)
    heldout_ids: pd.Index       # lines evaluated (test set, or dev fold 0 in the dry run)
    folds: pd.Series            # train_ids -> fold (fold models, D7)
    lineage: pd.Series          # every ModelID in splits.csv
    patient: pd.Series
    label: str


def population(paths: Paths, dry: bool, dry_fold=0) -> Population:
    sp = pd.read_csv(paths.splits, index_col=0)
    sp.index = sp.index.astype(str)
    dev = sp[sp["role"] == "dev"]
    if dry:
        ho = dev.index[dev["cv_fold"].astype(int) == dry_fold]
        tr = dev.index[dev["cv_fold"].astype(int) != dry_fold]
        label = f"dry-run: dev fold {dry_fold} stands in for the test set; training = other dev folds"
    else:
        ho = sp.index[sp["role"] == "test"]
        tr = dev.index
        label = "test set (refit on all dev lines)"
    return Population(pd.Index(tr), pd.Index(ho), dev.loc[tr, "cv_fold"].astype(int),
                      sp["OncotreeLineage"].astype(str), sp["PatientID"].astype(str), label)


def default_feature_loader(names, ids):
    """scripts/run.py convention: X = concat of cache matrices prefixed '<name>:' (None if no features)."""
    if not names:
        return None
    mats = [io.load_matrix(fn).reindex(ids).add_prefix(f"{fn}:") for fn in names]
    X = pd.concat(mats, axis=1)
    if X.isna().any().any():
        raise ValueError("features have NaN for some requested lines")
    return X


# ============================================================================ fitting (scripts/run.py protocol)
def _common(paths: Paths):
    eval_genes = pd.Index(pd.read_csv(paths.eval_genes, dtype={"entrez": str})["entrez"])
    neg, _ = genesets.load_hart_neg()
    seed = int(json.loads(Path(paths.seeds).read_text())["split_seed"])
    return eval_genes, neg, seed


def _check_fit_log(start, paths: Paths):
    ids = set().union(*[set(e["ids"]) for e in RUN.FIT_LOG[start:]]) if len(RUN.FIT_LOG) > start else set()
    roles = guard.roles_of(sorted(ids), paths.splits)
    if (roles != "dev").any():
        raise RUN.LeakageError("non-dev IDs reached a fit")
    return len(ids)


def refit_predict(f: Finalist, pop: Population, paths: Paths, feature_loader=None, predict_ids=None,
                  params_extra=None):
    """Tune + fit ``f`` on pop.train_ids (minus exclusions) and predict ``predict_ids`` (default held-out).

    Same protocol as one outer fold of scripts/run.py (inner grouped tuning, sentinel). Seed =
    split_seed (+0). Returns dict(model, pred, chosen, fit_s, predict_s, n_train_ids).
    M5 (m5_avg): ``fit`` only stores the training data (average mode); its inputs (M1, M2 with nested local
    lambda, M3 with nested early stopping, exactly as in CV) are refitted on the training lines inside
    ``predict``, so the refit time is fit_s + predict_s. ``params_extra`` (not part of the config hash) only
    switches execution options such as M5's prediction cache (``{"cache": False}``: a genuine refit).
    """
    feature_loader = feature_loader or default_feature_loader
    cfg = f.cfg
    if params_extra:
        cfg = dict(cfg, params={**cfg.get("params", {}), **params_extra})
    eval_genes, neg, seed = _common(paths)
    excl = pd.Index(map(str, cfg.get("exclude_from_training", [])))
    tr = pop.train_ids[~pop.train_ids.isin(excl)]
    ho = pd.Index(pop.heldout_ids if predict_ids is None else predict_ids)
    Y_tr = guard.load_targets(list(tr), config_hash=f.config_hash, **paths.guard_kw()).astype(np.float64)
    X = feature_loader(cfg.get("features", []), tr.append(ho))
    sent = RUN.Sentinel(ho.union(excl), "refit")
    n0 = len(RUN.FIT_LOG)
    model, pred, best, table, t_fit, t_pred = RUN.fit_predict_fold(
        cfg["model"], cfg.get("params", {}), cfg.get("grid", {}), Y_tr, ho, pop.lineage, X, pop.patient,
        sent, eval_genes, neg, cfg.get("inner_folds", 5), seed, cfg.get("inner_min_obs", 30))
    n_train = _check_fit_log(n0, paths)
    for e in getattr(model, "input_fit_log_", None) or []:  # M5: its input refits happen inside predict
        RUN.FIT_LOG.append({"tag": f"m5_input:{e['input']}", "outer_fold": "refit", "ids": e["ids"]})
        bad = set(e["ids"]) & set(map(str, ho.union(excl)))
        if bad or (guard.roles_of(e["ids"], paths.splits) != "dev").any():
            raise RUN.LeakageError(f"M5 input {e['input']} was fitted on held-out / non-dev lines")
    return {"model": model, "pred": pred, "chosen": best, "fit_s": t_fit, "predict_s": t_pred,
            "n_train_ids": n_train, "train_ids": tr}


def _fold_meta(f: Finalist, folds: pd.Series, dry: bool):
    sig = hashlib.sha256(json.dumps(sorted(map(str, folds.index)) + [int(x) for x in folds.sort_index()]).encode())
    return {"config_hash": f.config_hash, "folds_sha256": sig.hexdigest(), "dry": bool(dry)}


def fold_model_tag(f: Finalist, dry: bool):
    return f"{'dry_' if dry else ''}{f.tag}"


def load_or_fit_fold_models(f: Finalist, pop: Population, paths: Paths, dry: bool, feature_loader=None,
                            target_loader=None, log=print):
    """The D7 fold models of ``f``: one per fold of pop.folds, trained on the other folds (exclusions
    applied). Loaded from cache/models/<tag>/ if persisted with the same config and folds, else fitted
    with scripts/run.py ``run_cv`` (same seeds as the CV run: seed + fold) and saved (modelstore).
    Returns (list of (fold, model, meta), info dict)."""
    tag = fold_model_tag(f, dry)
    d = Path(paths.model_dir) / tag
    meta = _fold_meta(f, pop.folds, dry)
    tl = target_loader or (lambda ids: guard.load_targets(list(ids), **paths.guard_kw()))
    mp = d / "meta.json"
    if mp.exists() and json.loads(mp.read_text()) == meta:
        log(f"  fold models of {f.name}: loaded from {paths.rel(d)}")
        return modelstore.load_fold_models(tag, paths.model_dir, tl), {"source": "persisted", "tag": tag}
    feature_loader = feature_loader or default_feature_loader
    eval_genes, neg, seed = _common(paths)
    cfg = {k: v for k, v in f.cfg.items() if k != "reference_baselines"}
    tr = pop.folds.index
    Y = guard.load_targets(list(tr), config_hash=f.config_hash, **paths.guard_kw()).astype(np.float64)
    X = feature_loader(cfg.get("features", []), tr)
    n0 = len(RUN.FIT_LOG)
    res = RUN.run_cv(Y, pop.folds, pop.patient.reindex(tr), pop.lineage.reindex(tr), X, cfg, eval_genes, neg,
                     seed, log=log)
    _check_fit_log(n0, paths)
    feats = [] if X is None else list(X.columns)
    for k, m in res["models"].items():
        if hasattr(m, "compact"):
            m.compact()
        modelstore.save_fold_model(tag, k, m, res["train_ids"][k], feats, paths.model_dir)
    mp.write_text(json.dumps(meta))
    info = {"source": "fitted", "tag": tag, "fit_s": res["fit_s"], "oof_pred": res["pred"]}
    return modelstore.load_fold_models(tag, paths.model_dir, tl), info


CV_REPRO_TOL = 1e-3


def fold_models_oof(models, f: Finalist, pop: Population, feature_loader=None) -> pd.DataFrame:
    """Each fold model predicts its own held-out fold (pop.folds == k): the run's out-of-fold predictions."""
    feature_loader = feature_loader or default_feature_loader
    parts = []
    for k, m, _ in models:
        ids = pd.Index(pop.folds.index[pop.folds == k])
        X = feature_loader(f.cfg.get("features", []), ids)
        lin = pop.lineage.reindex(ids)
        parts.append(m.predict(ids, lineage=lin, X=X.loc[ids] if isinstance(X, pd.DataFrame) else X))
    return pd.concat(parts)


def cv_reproduction(f: Finalist, oof: pd.DataFrame, paths: Paths) -> dict:
    """Fold models' OOF predictions vs the registered CV run's cache/runs/<cv_run_id>/pred.parquet."""
    if not f.cv_run_id:
        return {"status": "no cv_run_id in finalists file", "ok": False}
    p = Path(paths.cv_cache) / f.cv_run_id / "pred.parquet"
    if not p.exists():
        return {"cv_run_id": f.cv_run_id, "status": "cached CV predictions not found", "ok": False}
    P = pd.read_parquet(p)
    P.columns = P.columns.astype(str)
    o = oof.copy()
    o.columns = o.columns.astype(str)
    ids, g = o.index.intersection(P.index), o.columns.intersection(P.columns)
    d = float(np.abs(o.loc[ids, g].to_numpy(np.float64) - P.loc[ids, g].to_numpy(np.float64)).max())
    return {"cv_run_id": f.cv_run_id, "n_lines": int(len(ids)), "n_genes": int(len(g)), "max_abs_diff": d,
            "ok": bool(len(ids) == len(oof) and d <= CV_REPRO_TOL)}


def predict_folds(models, ids, lineage, X) -> dict:
    """{fold: DataFrame [ids x gene]} from loaded fold models."""
    ids = pd.Index(ids)
    lin = None if lineage is None else pd.Series(lineage).reindex(ids)
    out = {}
    for k, m, _ in models:
        Xp = X.loc[ids] if isinstance(X, pd.DataFrame) else X
        out[k] = m.predict(ids, lineage=lin, X=Xp)
    return out


def peak_rss():
    return RUN.peak_rss_mb()


def write_json(path, obj):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(obj, indent=1, default=_json_default))


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, (pd.Index, np.ndarray)):
        return list(o)
    if isinstance(o, Path):
        return o.as_posix()
    return str(o)


def clean_floats(o):
    """Recursively replace non-finite floats with None (strict JSON)."""
    if isinstance(o, dict):
        return {str(k): clean_floats(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean_floats(v) for v in o]
    if isinstance(o, (float, np.floating)):
        return float(o) if np.isfinite(o) else None
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o
