"""Static check (DESIGN v2 §9): only io.py and guard.py may read the CRISPR file directly.

Fails if any other module under src/ or scripts/ references the raw CRISPR file name, the cached
crispr parquet, or calls load_matrix with 'crispr'. Legacy allow-list (flagged, pre-guard code):
scripts/audit_data.py (pre-split data audit) and scripts/probe/ (flawed calibration probe).
"""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ALLOWED = {"src/essml/io.py", "src/essml/guard.py"}
LEGACY = {"scripts/audit_data.py"}
LEGACY_DIRS = ("scripts/probe/",)
PATTERN = re.compile(r"CRISPRGeneEffect|crispr\.parquet|CRISPR_PARQUET|_read_rows|load_matrix\(\s*['\"]crispr['\"]",
                     re.I)


def test_no_direct_crispr_reads():
    offenders, legacy_hits = [], []
    for p in sorted([*ROOT.glob("src/**/*.py"), *ROOT.glob("scripts/**/*.py")]):
        rel = p.relative_to(ROOT).as_posix()
        hits = [i + 1 for i, line in enumerate(p.read_text(encoding="utf-8").splitlines())
                if PATTERN.search(line)]
        if not hits or rel in ALLOWED:
            continue
        if rel in LEGACY or rel.startswith(LEGACY_DIRS):
            legacy_hits.append(rel)
            continue
        # a docstring/comment mention of the file name is still flagged: keep such modules clean
        offenders.append(f"{rel}:{hits}")
    assert not offenders, f"direct CRISPR reads outside io/guard: {offenders}"
    if legacy_hits:
        print(f"[flag] legacy modules reading CRISPR directly: {legacy_hits}")
