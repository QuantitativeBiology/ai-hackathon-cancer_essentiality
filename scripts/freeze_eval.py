"""Freeze the evaluation code: write sha256 of the frozen modules to configs/frozen_hashes.json.

Hashes are computed on content with CRLF normalised to LF, so a Windows checkout with
core.autocrlf does not break the check. tests/test_frozen.py compares against this file.
Run only at a freeze (DESIGN §9); the PI tags the commit (freeze-eval).
"""
import datetime as dt
import hashlib
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FROZEN = ["src/essml/metrics.py", "src/essml/evaluate.py", "src/essml/genesets.py",
          "src/essml/guard.py"]
OUT = ROOT / "configs" / "frozen_hashes.json"


def sha256_lf(path):
    return hashlib.sha256(Path(path).read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def main():
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True,
                          text=True).stdout.strip()
    out = {"_comment": "sha256 (CRLF->LF normalised) of frozen evaluation modules; see DESIGN §9. "
                       "Any change requires a RESEARCH_TREE note and re-running baselines.",
           "frozen_at": dt.datetime.now().isoformat(timespec="seconds"),
           "parent_commit": head,
           "files": {f: sha256_lf(ROOT / f) for f in FROZEN}}
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out["files"], indent=1))


if __name__ == "__main__":
    main()
