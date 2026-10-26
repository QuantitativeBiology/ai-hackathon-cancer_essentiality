"""Frozen evaluation code must match configs/frozen_hashes.json (DESIGN §9)."""
import json
from pathlib import Path

from freeze_eval import FROZEN, sha256_lf

ROOT = Path(__file__).resolve().parents[1]


def test_frozen_hashes():
    ref = json.loads((ROOT / "configs" / "frozen_hashes.json").read_text())["files"]
    assert sorted(ref) == sorted(FROZEN)
    changed = [f for f in FROZEN if sha256_lf(ROOT / f) != ref[f]]
    assert not changed, f"frozen evaluation code changed: {changed} (log in RESEARCH_TREE, re-run baselines)"
