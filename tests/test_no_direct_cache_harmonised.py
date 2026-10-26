"""Static + behavioural check: cache/harmonised (set-B test profiles) is read only via essml.guard.

Any module under src/ or scripts/ other than guard.py that refers to the cache/harmonised directory
may write there, but a line that both calls a reader (read_parquet, read_csv, read_table, ParquetFile,
read_pickle, np.load) and mentions the directory (the literal path, or a variable assigned from it)
fails. Reads go through ``guard.load_harmonised``.
"""
import re
from pathlib import Path

import pandas as pd
import pytest

from essml import guard

ROOT = Path(__file__).resolve().parents[1]
DIR_DEF = re.compile(r"""^\s*(\w+)\s*=\s*.*["']cache["']\s*/\s*["']harmonised["']|cache/harmonised""", re.M)
VAR_DEF = re.compile(r"""^\s*(\w+)\s*=\s*.*["']cache["']\s*/\s*["']harmonised["']""", re.M)
READ = re.compile(r"read_parquet|read_csv|read_table|ParquetFile|read_pickle|np\.load|open\(")


def test_no_direct_cache_harmonised_reads():
    offenders, seen = [], []
    for p in sorted([*ROOT.glob("src/**/*.py"), *ROOT.glob("scripts/**/*.py")]):
        rel = p.relative_to(ROOT).as_posix()
        text = p.read_text(encoding="utf-8")
        if rel == "src/essml/guard.py" or not DIR_DEF.search(text):
            continue
        seen.append(rel)
        names = {m.group(1) for m in VAR_DEF.finditer(text)}
        tokens = re.compile(r"cache/harmonised|" + "|".join(rf"\b{re.escape(n)}\b" for n in names)
                            if names else r"cache/harmonised")
        for i, line in enumerate(text.splitlines(), 1):
            if "load_harmonised" in line or line.lstrip().startswith("#"):
                continue
            if READ.search(line) and tokens.search(line):
                offenders.append(f"{rel}:{i}")
    assert seen, "expected the harmonisation scripts to reference cache/harmonised"
    assert not offenders, f"direct cache/harmonised reads outside guard: {offenders}"


@pytest.fixture
def prof(tmp_path):
    d = tmp_path / "harm"
    d.mkdir()
    pd.DataFrame({"final": [1.0, 2.0, 3.0]}, index=pd.Index(["1", "2", "3"], name="entrez")
                 ).to_parquet(d / "scr.parquet")
    return d, tmp_path


def test_controls_only_filters_and_logs(prof):
    d, tmp = prof
    log = tmp / "log.csv"
    df = guard.load_harmonised("scr", controls_only=True, controls={"1", "3"}, harmonised_dir=d,
                               log_path=log)
    assert list(df.index) == ["1", "3"]
    row = pd.read_csv(log, dtype=str).iloc[0]
    assert row["access_type"] == "control_only" and row["scope"] == "setB:scr"


def test_target_read_refused_without_unblinding(prof, monkeypatch):
    d, tmp = prof
    monkeypatch.delenv("ESSML_UNBLIND", raising=False)
    repo = tmp / "repo"
    repo.mkdir()
    with pytest.raises(guard.BlindingError):
        guard.load_harmonised("scr", harmonised_dir=d, log_path=tmp / "log.csv", repo_root=repo)
    row = pd.read_csv(tmp / "log.csv", dtype=str).iloc[0]
    assert row["granted"] == "False" and row["access_type"] == "targets"


def test_aggregate_purpose_only_before_transfer_freeze(prof, monkeypatch):
    import subprocess
    d, tmp = prof
    monkeypatch.delenv("ESSML_UNBLIND", raising=False)
    repo = tmp / "repo"
    (repo / "configs").mkdir(parents=True)
    kw = dict(purpose="replicate_ceiling_aggregate", harmonised_dir=d, log_path=tmp / "log.csv",
              repo_root=repo)
    subprocess.run(["git", "init", "-q"], cwd=repo, check=True)
    assert len(guard.load_harmonised("scr", **kw)) == 3  # criteria not yet committed
    assert pd.read_csv(tmp / "log.csv", dtype=str).iloc[0]["access_type"] == "aggregate_only"
    (repo / "configs" / "transfer_criteria.json").write_text("{}")
    subprocess.run(["git", "add", "."], cwd=repo, check=True)
    subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t", "commit", "-qm", "x"], cwd=repo,
                   check=True)
    with pytest.raises(guard.BlindingError):  # frozen: full unblinding needed
        guard.load_harmonised("scr", **kw)
    with pytest.raises(guard.BlindingError):  # an unknown purpose never qualifies
        guard.load_harmonised("scr", **{**kw, "purpose": "other"})


def test_anchor_must_be_dev(prof):
    d, tmp = prof
    sp = tmp / "splits.csv"
    pd.DataFrame({"ModelID": ["D1", "T1"], "PatientID": ["p", "q"], "OncotreeLineage": ["a", "a"],
                  "role": ["dev", "test"], "cv_fold": [0, -1]}).to_csv(sp, index=False)
    kw = dict(harmonised_dir=d, log_path=tmp / "log.csv", splits_path=sp)
    assert len(guard.load_harmonised("scr", anchor_model_id="D1", **kw)) == 3
    with pytest.raises(guard.BlindingError):
        guard.load_harmonised("scr", anchor_model_id="T1", **kw)
