"""Test-set guard for CRISPR targets (FROZEN: hashed in configs/frozen_hashes.json; DESIGN.md §9).

All CRISPR target access goes through this module (``io.load_matrix('crispr')`` refuses without
``unguarded=True``, reserved for the legacy pre-split audit; tests/test_no_direct_crispr_reads.py
fails if any other module reads the CRISPR file).

Roles come from ``configs/splits.csv``: 'dev', 'test', 'noncancerous' (set A); IDs absent from it
are 'unassigned'. External published screens (set B) are files, read via ``load_external_screen``.

Full target release for non-dev data (``load_targets`` / ``load_external_screen``) requires ALL of
  (a) env ESSML_UNBLIND == "1";
  (b) ``configs/finalists.json`` exists and is tracked in git (``git ls-files``);
  (c) a clean git tree (no modified tracked files, so finalists.json equals its committed version);
  (d) ``config_hash`` is listed in finalists.json (any "config_hash" value in the JSON);
  (e) no earlier granted access for the same (config_hash, scope): scope 'test' for the cancer test
      set, 'setA:<screen>' / 'setB:<screen>' / 'other:<screen>' otherwise (one access per finalist,
      and per finalist x screen; staged per screen, so a request may not mix scopes);
  (f) for set A / set B / other scopes only: ``configs/harmonisation_frozen.json`` and
      ``configs/transfer_criteria.json`` are tracked in git.
Logged exceptions (no unblinding needed):
  - ``observed_counts(ids)``: per-gene non-NaN COUNTS only (access_type 'mask_only');
  - control-gene reads (CEGv2 / NEGv1 / scaling-set = dev common essentials) for harmonisation and QC:
    ``load_control_targets(ids)`` for DepMap rows and ``load_external_screen(..., controls_only=True)``
    (access_type 'control_only').
Every non-dev access (granted or refused) is appended to ``results/test_access_log.csv``.
``screened_model_ids()`` returns the CRISPR row index only (no values).

Only requested rows are read from parquet (row filter), so refused rows never enter memory.
Features (omics) of test lines are inputs and are not guarded here.
"""
from __future__ import annotations

import csv
import datetime as _dt
import inspect
import json
import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SPLITS = ROOT / "configs" / "splits.csv"
FINALISTS = "configs/finalists.json"  # relative to repo root
TRANSFER_FREEZE = ["configs/harmonisation_frozen.json", "configs/transfer_criteria.json"]
CRISPR_PARQUET = ROOT / "cache" / "crispr.parquet"
ACCESS_LOG = ROOT / "results" / "test_access_log.csv"
LOG_COLS = ["timestamp", "caller", "access_type", "scope", "config_hash", "n_ids", "n_genes",
            "roles", "granted", "reason"]


class BlindingError(PermissionError):
    """Raised when non-dev CRISPR targets are requested without a valid unblinding."""


def roles_of(model_ids, splits_path=SPLITS) -> pd.Series:
    """Role per ModelID ('dev'/'test'/'noncancerous'; 'unassigned' if not in splits.csv)."""
    s = pd.read_csv(splits_path, index_col=0)["role"]
    ids = pd.Index([str(i) for i in model_ids])
    return pd.Series(s.reindex(ids).fillna("unassigned").to_numpy(), index=ids, name="role")


def _git(args, repo_root):
    return subprocess.run(["git"] + args, cwd=repo_root, capture_output=True, text=True)


def _tracked(rel, repo_root) -> bool:
    return (Path(repo_root) / rel).exists() and \
        _git(["ls-files", "--error-unmatch", rel], repo_root).returncode == 0


def finalist_hashes(repo_root=ROOT, rel=FINALISTS) -> set:
    """All string values stored under a "config_hash" key anywhere in finalists.json."""
    out = set()

    def walk(o):
        if isinstance(o, dict):
            for k, v in o.items():
                if k == "config_hash" and isinstance(v, str):
                    out.add(v)
                else:
                    walk(v)
        elif isinstance(o, list):
            for v in o:
                walk(v)

    walk(json.loads((Path(repo_root) / rel).read_text()))
    return out


def _granted_before(config_hash, scope, log_path) -> bool:
    if not Path(log_path).exists():
        return False
    log = pd.read_csv(log_path, dtype=str)
    if "config_hash" not in log or "scope" not in log:
        return False
    m = ((log["access_type"] == "targets") & (log["granted"] == "True")
         & (log["config_hash"] == config_hash) & (log["scope"] == scope))
    return bool(m.any())


def unblinding_check(config_hash=None, scope="test", repo_root=ROOT, log_path=ACCESS_LOG,
                     rel=FINALISTS):
    """(allowed, reason) for releasing non-dev targets of ``scope`` (conditions a-f above)."""
    repo_root = Path(repo_root)
    if scope.endswith(":"):
        return False, "screen name required for per-screen accounting"
    if os.environ.get("ESSML_UNBLIND") != "1":
        return False, "ESSML_UNBLIND!=1"
    if not (repo_root / rel).exists():
        return False, "finalists.json missing"
    try:
        if not _tracked(rel, repo_root):
            return False, "finalists.json not tracked"
        st = _git(["status", "--porcelain", "--untracked-files=no"], repo_root)
        if st.returncode != 0 or st.stdout.strip():
            return False, "dirty tree"
        if scope != "test" and not all(_tracked(f, repo_root) for f in TRANSFER_FREEZE):
            return False, "transfer freeze files not committed"
    except OSError:  # git not available -> cannot verify -> refuse
        return False, "git unavailable"
    if not config_hash or config_hash not in finalist_hashes(repo_root, rel):
        return False, "config hash not a finalist"
    if _granted_before(config_hash, scope, log_path):
        return False, "finalist already accessed this scope"
    return True, "ok"


def _caller(depth=2) -> str:
    st = inspect.stack()
    if len(st) <= depth:
        return "?"
    f = st[depth]
    return f"{Path(f.filename).name}:{f.function}:{f.lineno}"


def _log(access_type, roles: pd.Series, granted: bool, log_path=ACCESS_LOG, caller=None,
         config_hash="", reason="", scope="", n_genes=""):
    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    vc = roles.value_counts()
    row = {"timestamp": _dt.datetime.now().isoformat(timespec="seconds"),
           "caller": caller or _caller(3), "access_type": access_type, "scope": scope,
           "config_hash": config_hash or "", "n_ids": int(len(roles)), "n_genes": n_genes,
           "roles": ";".join(f"{k}:{int(v)}" for k, v in sorted(vc.items())), "granted": granted,
           "reason": reason}
    new = not log_path.exists()
    with open(log_path, "a", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=LOG_COLS)
        if new:
            w.writeheader()
        w.writerow(row)


def _check_paths(splits_path, crispr_path):
    """Refuse a substitute splits file combined with the real CRISPR cache (red-team G0 #1).

    Toy paths are allowed for tests only when *both* are substitutes.
    """
    real_crispr = Path(crispr_path).resolve() == CRISPR_PARQUET.resolve()
    real_splits = Path(splits_path).resolve() == Path(SPLITS).resolve()
    if real_crispr and not real_splits:
        raise BlindingError("the real CRISPR cache may only be read with the committed configs/splits.csv")


def _read_rows(model_ids, crispr_path, columns=None):
    ids = [str(i) for i in model_ids]
    df = pd.read_parquet(crispr_path, filters=[("ModelID", "in", ids)],
                         columns=None if columns is None else list(columns))
    missing = sorted(set(ids) - set(df.index))
    if missing:
        raise KeyError(f"{len(missing)} ModelIDs not in CRISPR matrix, e.g. {missing[:3]}")
    return df.loc[ids]


def _scope(roles: pd.Series, screen):
    """Access scope of a request; a request may not mix the test set with other non-dev data."""
    nd = set(roles[roles != "dev"])
    if not nd:
        return "dev"
    if "test" in nd:
        if len(nd) > 1:
            raise ValueError("stage access per scope: do not mix test with other non-dev IDs")
        return "test"
    if len(nd) > 1:
        raise ValueError("stage access per screen: do not mix set A and unassigned IDs")
    # no screen name -> scope ends with ':' and is refused (logged) by unblinding_check
    return {"noncancerous": "setA", "unassigned": "other"}[nd.pop()] + f":{screen or ''}"


def load_targets(model_ids, config_hash=None, screen=None, splits_path=SPLITS,
                 crispr_path=CRISPR_PARQUET, log_path=ACCESS_LOG, repo_root=ROOT) -> pd.DataFrame:
    """CRISPR gene effects (float32, index=ModelID in request order, columns=Entrez strings).

    Raises BlindingError for any non-dev ID unless unblinding is allowed (see module doc).
    """
    _check_paths(splits_path, crispr_path)
    roles = roles_of(model_ids, splits_path)
    if len(roles) and roles.index.duplicated().any():
        raise ValueError("duplicated ModelIDs requested")
    scope = _scope(roles, screen)
    if scope != "dev":
        ok, reason = unblinding_check(config_hash, scope, repo_root, log_path)
        _log("targets", roles, ok, log_path, caller=_caller(2), config_hash=config_hash,
             reason=reason, scope=scope)
        if not ok:
            bad = roles[roles != "dev"]
            raise BlindingError(
                f"refused CRISPR targets for {len(bad)} non-dev IDs "
                f"({bad.value_counts().to_dict()}), scope {scope}: {reason}. Requires "
                f"ESSML_UNBLIND=1, a clean tree and a finalist config hash in a committed {FINALISTS}")
    return _read_rows(roles.index, crispr_path)


def observed_counts(model_ids, splits_path=SPLITS, crispr_path=CRISPR_PARQUET,
                    log_path=ACCESS_LOG) -> pd.Series:
    """Per-gene number of non-NaN CRISPR entries over ``model_ids`` (int Series, index=Entrez).

    Mask-only: values are never returned. Allowed without unblinding; non-dev access is logged
    as access_type='mask_only'.
    """
    _check_paths(splits_path, crispr_path)
    roles = roles_of(model_ids, splits_path)
    if (roles != "dev").any():
        _log("mask_only", roles, True, log_path, caller=_caller(2), reason="counts only")
    df = _read_rows(roles.index, crispr_path)
    counts = pd.Series(np.isfinite(df.to_numpy()).sum(0), index=df.columns, name="n_obs")
    del df
    return counts


def control_genes(gene_classes_csv=ROOT / "configs" / "gene_classes.csv") -> set:
    """Control set = Hart CEGv2 ∪ Hart NEGv1 ∪ dev common essentials (scaling set; DESIGN §7.2)."""
    from . import genesets
    ceg, _ = genesets.load_hart_ceg()
    neg, _ = genesets.load_hart_neg()
    out = set(ceg) | set(neg)
    if Path(gene_classes_csv).exists():
        gc = genesets.load_gene_classes(gene_classes_csv)
        out |= set(gc.index[gc["common_essential"]])
    return out


def _check_controls(genes, controls):
    genes = [str(g) for g in genes]
    bad = sorted(set(genes) - set(controls))
    if bad:
        raise BlindingError(f"{len(bad)} requested genes are not control genes, e.g. {bad[:3]}")
    return genes


def load_control_targets(model_ids, genes=None, controls=None, splits_path=SPLITS,
                         crispr_path=CRISPR_PARQUET, log_path=ACCESS_LOG) -> pd.DataFrame:
    """CRISPR values of CONTROL genes only (default: all control genes present), DepMap dev lines only.

    DepMap test/set-A lines are refused (red-team G0 #2): the control-gene exception of DESIGN §9 covers
    external screens (load_external_screen), not held-out DepMap lines, whose control genes carry most
    P2a positives.
    """
    import pyarrow.parquet as pq
    controls = control_genes() if controls is None else set(map(str, controls))
    if genes is None:  # gene names from the parquet schema only (no values read)
        genes = [g for g in pq.ParquetFile(crispr_path).schema_arrow.names if g in controls]
    genes = _check_controls(genes, controls)
    _check_paths(splits_path, crispr_path)
    roles = roles_of(model_ids, splits_path)
    if (roles != "dev").any():
        _log("control_only", roles, False, log_path, caller=_caller(2), reason="refused: non-dev DepMap IDs",
             n_genes=len(genes))
        raise BlindingError("control-gene reads are allowed for DepMap dev lines only")
    return _read_rows(roles.index, crispr_path, columns=genes)


def load_external_screen(path, screen, reader, config_hash=None, controls_only=False,
                         controls=None, log_path=ACCESS_LOG, repo_root=ROOT) -> pd.DataFrame:
    """Read a set-B external screen through the guard.

    reader(path) -> DataFrame with Entrez-string gene COLUMNS (rows = screen arms / samples).
    controls_only=True returns control-gene columns only (logged 'control_only', no unblinding);
    otherwise the full table is released only per conditions a-f with scope 'setB:<screen>'.
    """
    scope = f"setB:{screen}"
    roles = pd.Series(["external"], index=[screen])
    if controls_only:
        controls = control_genes() if controls is None else set(map(str, controls))
        df = reader(path)
        keep = [c for c in df.columns if str(c) in controls]
        _log("control_only", roles, True, log_path, caller=_caller(2), reason="control genes",
             scope=scope, n_genes=len(keep))
        return df[keep]
    ok, reason = unblinding_check(config_hash, scope, repo_root, log_path)
    _log("targets", roles, ok, log_path, caller=_caller(2), config_hash=config_hash,
         reason=reason, scope=scope)
    if not ok:
        raise BlindingError(f"refused external screen {screen}: {reason}")
    return reader(path)


def screened_model_ids(crispr_path=CRISPR_PARQUET) -> pd.Index:
    """ModelIDs with a CRISPR profile (row index only; no target values are read)."""
    import pyarrow.parquet as pq
    ids = pq.read_table(crispr_path, columns=["ModelID"]).column("ModelID").to_pylist()
    return pd.Index(ids, name="ModelID")
