"""Test guard: refuses non-dev targets, logs every non-dev access, mask-only exception (DESIGN §9)."""
import json
import subprocess

import numpy as np
import pandas as pd
import pytest

from essml import guard


@pytest.fixture
def toy(tmp_path):
    splits = pd.DataFrame({"ModelID": ["D1", "D2", "T1", "N1"], "PatientID": ["p1", "p2", "p3", "p4"],
                           "OncotreeLineage": ["a"] * 4, "role": ["dev", "dev", "test", "noncancerous"],
                           "cv_fold": [0, 1, -1, -1]})
    sp = tmp_path / "splits.csv"
    splits.to_csv(sp, index=False)
    Y = pd.DataFrame({"10": [1.0, 2.0, 3.0, 4.0, 5.0], "20": [np.nan, 1.0, np.nan, 1.0, 1.0]},
                     index=pd.Index(["D1", "D2", "T1", "N1", "X1"], name="ModelID"))
    cp = tmp_path / "crispr.parquet"
    Y.to_parquet(cp)
    repo = tmp_path / "repo"
    repo.mkdir()
    kw = dict(splits_path=sp, crispr_path=cp, log_path=tmp_path / "log.csv")
    return kw, repo


def _git(repo, *a):
    subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t", *a], cwd=repo, check=True,
                   capture_output=True)


def test_dev_allowed_not_logged(toy):
    kw, repo = toy
    Y = guard.load_targets(["D2", "D1"], repo_root=repo, **kw)
    assert list(Y.index) == ["D2", "D1"] and Y.loc["D1", "10"] == 1.0 and np.isnan(Y.loc["D1", "20"])
    assert not kw["log_path"].exists()


@pytest.mark.parametrize("ids", [["D1", "T1"], ["N1"], ["X1"]])  # test, set A, not in splits (set B)
def test_refuses_non_dev_and_logs(toy, ids, monkeypatch):
    kw, repo = toy
    monkeypatch.delenv("ESSML_UNBLIND", raising=False)
    with pytest.raises(guard.BlindingError):
        guard.load_targets(ids, repo_root=repo, **kw)
    log = pd.read_csv(kw["log_path"], dtype=str)
    assert len(log) == 1 and log.loc[0, "granted"] == "False"
    assert log.loc[0, "access_type"] == "targets" and int(log.loc[0, "n_ids"]) == len(ids)
    assert "test_guard.py" in log.loc[0, "caller"]


def test_unblind_env_alone_is_not_enough(toy, monkeypatch):
    kw, repo = toy
    monkeypatch.setenv("ESSML_UNBLIND", "1")
    with pytest.raises(guard.BlindingError, match="finalists.json missing"):
        guard.load_targets(["T1"], config_hash="abc", repo_root=repo, **kw)


def test_unblinding_conditions(toy, monkeypatch):
    kw, repo = toy
    _git(repo, "init", "-q")
    (repo / "configs").mkdir()
    fj = repo / "configs" / "finalists.json"
    fj.write_text(json.dumps({"headline": {"name": "m", "config_hash": "abc"}}))
    monkeypatch.setenv("ESSML_UNBLIND", "1")
    with pytest.raises(guard.BlindingError, match="not tracked"):      # untracked
        guard.load_targets(["T1"], config_hash="abc", repo_root=repo, **kw)
    _git(repo, "add", "configs/finalists.json")
    _git(repo, "commit", "-q", "-m", "finalists")
    with pytest.raises(guard.BlindingError, match="not a finalist"):   # wrong config hash
        guard.load_targets(["T1"], config_hash="zzz", repo_root=repo, **kw)
    fj.write_text(json.dumps({"headline": {"name": "m", "config_hash": "zzz"}}))
    with pytest.raises(guard.BlindingError, match="dirty tree"):       # modified tracked file
        guard.load_targets(["T1"], config_hash="zzz", repo_root=repo, **kw)
    _git(repo, "checkout", "--", "configs/finalists.json")
    monkeypatch.setenv("ESSML_UNBLIND", "0")
    with pytest.raises(guard.BlindingError):                            # env not set
        guard.load_targets(["T1"], config_hash="abc", repo_root=repo, **kw)
    monkeypatch.setenv("ESSML_UNBLIND", "1")
    Y = guard.load_targets(["T1", "D1"], config_hash="abc", repo_root=repo, **kw)
    assert Y.loc["T1", "10"] == 3.0
    with pytest.raises(guard.BlindingError, match="already accessed"):  # one access per finalist
        guard.load_targets(["T1"], config_hash="abc", repo_root=repo, **kw)
    log = pd.read_csv(kw["log_path"], dtype=str)
    assert list(log["granted"]) == ["False"] * 4 + ["True", "False"]
    assert log.iloc[4]["roles"] == "dev:1;test:1" and log.iloc[4]["config_hash"] == "abc"
    assert log.iloc[4]["scope"] == "test"
    # set A needs a screen name, and the transfer-freeze files committed, before release
    with pytest.raises(guard.BlindingError, match="screen name"):
        guard.load_targets(["N1"], config_hash="abc", repo_root=repo, **kw)
    with pytest.raises(guard.BlindingError, match="transfer freeze"):
        guard.load_targets(["N1"], config_hash="abc", screen="rpe1", repo_root=repo, **kw)
    for f in guard.TRANSFER_FREEZE:
        (repo / f).write_text("{}")
    _git(repo, "add", *guard.TRANSFER_FREEZE)
    _git(repo, "commit", "-q", "-m", "transfer freeze")
    Y = guard.load_targets(["N1"], config_hash="abc", screen="rpe1", repo_root=repo, **kw)
    assert Y.loc["N1", "10"] == 4.0
    with pytest.raises(guard.BlindingError, match="already accessed"):  # per finalist x screen
        guard.load_targets(["N1"], config_hash="abc", screen="rpe1", repo_root=repo, **kw)
    with pytest.raises(ValueError, match="mix"):
        guard.load_targets(["N1", "T1"], config_hash="abc", screen="x", repo_root=repo, **kw)


def test_external_screen(toy, tmp_path, monkeypatch):
    kw, repo = toy
    monkeypatch.delenv("ESSML_UNBLIND", raising=False)
    reader = lambda p: pd.DataFrame({"10": [-1.0], "20": [0.1], "30": [0.0]})  # noqa: E731
    c = guard.load_external_screen("x.csv", "s1", reader, controls_only=True, controls={"10", "30"},
                                   log_path=kw["log_path"], repo_root=repo)
    assert list(c.columns) == ["10", "30"]
    with pytest.raises(guard.BlindingError):
        guard.load_external_screen("x.csv", "s1", reader, config_hash="abc",
                                   log_path=kw["log_path"], repo_root=repo)
    log = pd.read_csv(kw["log_path"], dtype=str)
    assert list(log["access_type"]) == ["control_only", "targets"]
    assert list(log["scope"]) == ["setB:s1", "setB:s1"] and list(log["granted"]) == ["True", "False"]


def test_control_targets_only_control_genes(toy, monkeypatch):
    kw, repo = toy
    monkeypatch.delenv("ESSML_UNBLIND", raising=False)
    Y = guard.load_control_targets(["D1"], controls={"20"}, **kw)
    assert list(Y.columns) == ["20"] and list(Y.index) == ["D1"]
    with pytest.raises(guard.BlindingError):  # non-control gene
        guard.load_control_targets(["D1"], genes=["10"], controls={"20"}, **kw)
    with pytest.raises(guard.BlindingError):  # held-out DepMap line: refused and logged
        guard.load_control_targets(["T1", "D1"], controls={"20"}, **kw)
    log = pd.read_csv(kw["log_path"], dtype=str)
    assert log.loc[0, "access_type"] == "control_only" and log.loc[0, "granted"] == "False"


def test_fake_splits_with_real_crispr_refused(tmp_path):
    fake = tmp_path / "splits.csv"
    pd.DataFrame({"ModelID": ["ACH-000001"], "role": ["dev"]}).to_csv(fake, index=False)
    with pytest.raises(guard.BlindingError):
        guard.load_targets(["ACH-000001"], splits_path=fake)


def test_observed_counts_mask_only(toy, monkeypatch):
    kw, repo = toy
    monkeypatch.delenv("ESSML_UNBLIND", raising=False)
    c = guard.observed_counts(["T1", "N1"], **kw)
    assert c.to_dict() == {"10": 2, "20": 1}
    log = pd.read_csv(kw["log_path"], dtype=str)
    assert log.loc[0, "access_type"] == "mask_only" and log.loc[0, "granted"] == "True"


def test_screened_ids_and_loader_hook(toy):
    kw, _ = toy
    assert list(guard.screened_model_ids(kw["crispr_path"])) == ["D1", "D2", "T1", "N1", "X1"]
    from essml import io
    with pytest.raises(PermissionError):
        io.load_matrix("crispr")
