"""Git must never track raw DepMap data, caches or large files (PROMPT.md: never commit restricted data)."""
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAX_BYTES = 5_000_000
RAW_NAMES = re.compile(r"CRISPRGeneEffect|OmicsExpression|OmicsCN|OmicsSomatic|OmicsFusion|(^|/)Model(Condition)?\.csv$")
BINARY = re.compile(r"\.(parquet|feather|npy|npz|pkl|h5|h5ad)$")


def _tracked():
    out = subprocess.run(["git", "ls-files", "-z"], cwd=ROOT, capture_output=True, check=True).stdout
    return [p for p in out.decode().split("\0") if p]


def test_no_data_or_cache_paths():
    bad = [p for p in _tracked() if p.startswith(("data/", "cache/")) or RAW_NAMES.search(p) or BINARY.search(p)]
    assert not bad, bad


def test_no_large_files():
    big = [(p, (ROOT / p).stat().st_size) for p in _tracked() if (ROOT / p).is_file()
           and (ROOT / p).stat().st_size > MAX_BYTES]
    assert not big, big
