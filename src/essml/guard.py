"""Test-set guard for CRISPR targets (FROZEN: hashed in configs/frozen_hashes.json; DESIGN.md §9).

All CRISPR target access goes through this module (``io.load_matrix('crispr')`` refuses without
``unguarded=True``, reserved for the pre-split audit scripts).

- ``load_targets(ids)`` returns CRISPR rows. Any ModelID whose role in ``configs/splits.csv`` is
  not 'dev' (test, noncancerous, or absent from splits.csv) is REFUSED unless
  (a) env ESSML_UNBLIND == "1", AND (b) ``configs/finalists.json`` exists, is tracked in git
  (``git ls-files``) and has no uncommitted changes vs HEAD.
- ``observed_counts(ids)`` returns per-gene non-NaN COUNTS only (never values). It is allowed
  without unblinding (needed to fix the evaluable gene set) and logged as access_type='mask_only'.
- Every access touching a non-dev ID (granted or refused) is appended to
  ``results/test_access_log.csv``: timestamp, caller, access_type, n_ids, roles, granted.

Only the requested rows are read from parquet (row filter), so refused rows never enter memory.
Features (omics) of test lines are inputs and are not guarded here.
"""
from __future__ import annotations

import csv
import datetime as _dt
import inspect
import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SPLITS = ROOT / "configs" / "splits.csv"
FINALISTS = "configs/finalists.json"  # relative to repo root
CRISPR_PARQUET = ROOT / "cache" / "crispr.parquet"
ACCESS_LOG = ROOT / "results" / "test_access_log.csv"
LOG_COLS = ["timestamp", "caller", "access_type", "n_ids", "roles", "granted"]


class BlindingError(PermissionError):
    """Raised when non-dev CRISPR targets are requested without a valid unblinding."""


def roles_of(model_ids, splits_path=SPLITS) -> pd.Series:
    """Role per ModelID ('dev'/'test'/'noncancerous'; 'unassigned' if not in splits.csv)."""
    s = pd.read_csv(splits_path, index_col=0)["role"]
    ids = pd.Index([str(i) for i in model_ids])
    return pd.Series(s.reindex(ids).fillna("unassigned").to_numpy(), index=ids, name="role")


def finalists_committed(repo_root=ROOT, rel=FINALISTS) -> bool:
    """True iff ``rel`` exists, is tracked (git ls-files) and equals its HEAD version."""
    repo_root = Path(repo_root)
    if not (repo_root / rel).exists():
        return False
    try:
        ls = subprocess.run(["git", "ls-files", "--error-unmatch", rel], cwd=repo_root,
                            capture_output=True, text=True)
        if ls.returncode != 0:
            return False
        diff = subprocess.run(["git", "diff", "--quiet", "HEAD", "--", rel], cwd=repo_root,
                              capture_output=True, text=True)
        return diff.returncode == 0
    except OSError:  # git not available -> cannot verify -> refuse
        return False


def unblinding_allowed(repo_root=ROOT) -> bool:
    return os.environ.get("ESSML_UNBLIND") == "1" and finalists_committed(repo_root)


def _caller(depth=2) -> str:
    st = inspect.stack()
    if len(st) <= depth:
        return "?"
    f = st[depth]
    return f"{Path(f.filename).name}:{f.function}:{f.lineno}"


def _log(access_type, roles: pd.Series, granted: bool, log_path=ACCESS_LOG, caller=None):
    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    vc = roles.value_counts()
    row = {"timestamp": _dt.datetime.now().isoformat(timespec="seconds"),
           "caller": caller or _caller(3), "access_type": access_type, "n_ids": int(len(roles)),
           "roles": ";".join(f"{k}:{int(v)}" for k, v in sorted(vc.items())), "granted": granted}
    new = not log_path.exists()
    with open(log_path, "a", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=LOG_COLS)
        if new:
            w.writeheader()
        w.writerow(row)


def _read_rows(model_ids, crispr_path):
    ids = [str(i) for i in model_ids]
    df = pd.read_parquet(crispr_path, filters=[("ModelID", "in", ids)])
    missing = sorted(set(ids) - set(df.index))
    if missing:
        raise KeyError(f"{len(missing)} ModelIDs not in CRISPR matrix, e.g. {missing[:3]}")
    return df.loc[ids]


def load_targets(model_ids, splits_path=SPLITS, crispr_path=CRISPR_PARQUET,
                 log_path=ACCESS_LOG, repo_root=ROOT) -> pd.DataFrame:
    """CRISPR gene effects (float32, index=ModelID in request order, columns=Entrez strings).

    Raises BlindingError for any non-dev ID unless unblinding is allowed (see module doc).
    """
    roles = roles_of(model_ids, splits_path)
    if len(roles) and roles.index.duplicated().any():
        raise ValueError("duplicated ModelIDs requested")
    if (roles != "dev").any():
        ok = unblinding_allowed(repo_root)
        _log("targets", roles, ok, log_path, caller=_caller(2))
        if not ok:
            bad = roles[roles != "dev"]
            raise BlindingError(
                f"refused CRISPR targets for {len(bad)} non-dev IDs "
                f"({bad.value_counts().to_dict()}); requires ESSML_UNBLIND=1 and a committed "
                f"{FINALISTS}")
    return _read_rows(roles.index, crispr_path)


def observed_counts(model_ids, splits_path=SPLITS, crispr_path=CRISPR_PARQUET,
                    log_path=ACCESS_LOG) -> pd.Series:
    """Per-gene number of non-NaN CRISPR entries over ``model_ids`` (int Series, index=Entrez).

    Mask-only: values are never returned. Allowed without unblinding; non-dev access is logged
    as access_type='mask_only'.
    """
    roles = roles_of(model_ids, splits_path)
    if (roles != "dev").any():
        _log("mask_only", roles, True, log_path, caller=_caller(2))
    df = _read_rows(roles.index, crispr_path)
    counts = pd.Series(np.isfinite(df.to_numpy()).sum(0), index=df.columns, name="n_obs")
    del df
    return counts
