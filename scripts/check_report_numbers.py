"""Check every tagged number in REPORT.md against the results files (DESIGN §9, gate G7).

    python scripts/check_report_numbers.py [REPORT.md] [--strict]

Tag convention (an HTML comment immediately before the number; invisible when rendered):

    P1 on the test set was <!--n:results/test/summary.json#table/0/P1_variable-->0.308 ...
    <!--n-src:test=results/test/summary.json-->            (alias, defined once anywhere in the file)
    ... <!--n:test#table/0/dP1_vs_B1_variable-->+0.108 ...
    ... <!--n:registry:20261026-222117_m1_E_resB0_91e52f9d#P1_variable-->0.3081 ...
    ... <!--n:test#table/0/frac_beat_B1_variable|pct-->94.4% ...

* KEY = SOURCE#PATH[|FLAG...]. SOURCE is a repo-relative JSON file, an alias from ``<!--n-src:...-->``,
  or ``registry:<run_id>`` (a row of results/registry.csv; PATH = column). PATH is '/'-separated
  (dict keys or list indices), so keys containing dots (``n_r_gt_0.3_variable``) need no quoting.
* FLAGS: ``pct`` (value x 100), ``abs`` (absolute value), ``neg`` (negated).
* The number may carry a sign (+, -, U+2212), thousands commas, an exponent and a trailing %.
  It must equal the source value rounded to the printed number of decimals:
  |printed - value| <= 0.5 x 10^-decimals (plus 1e-12 for float noise).
* Every run id cited anywhere in the file (YYYYMMDD-HHMMSS_<name>_<hash8>) must exist in
  results/registry.csv.
* ``--strict``: also fail on untagged numbers with >= 2 decimals in prose (outside code blocks).

Exit status 1 on any mismatch, unresolved key, malformed tag or unknown run id.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NUM = r"[+\-−]?(?:\d{1,3}(?:,\d{3})+|\d+)(?:\.\d+)?(?:[eE][+\-]?\d+)?%?"
TAG = re.compile(r"<!--\s*n:(?P<key>[^>]*?)\s*-->\s*(?P<num>" + NUM + r")?")
SRC = re.compile(r"<!--\s*n-src:(?P<alias>[\w.\-]+)\s*=\s*(?P<path>[^>]*?)\s*-->")
RUN_ID = re.compile(r"\b\d{8}-\d{6}_[A-Za-z0-9_\-]+?_[0-9a-f]{8}\b")
UNTAGGED = re.compile(r"(?<![\w.\-])" + r"[+\-−]?\d+\.\d{2,}" + r"(?![\w.])")


def parse_number(s: str):
    """(value, decimals, is_percent) of a printed number."""
    pct = s.endswith("%")
    t = s.rstrip("%").replace("−", "-").replace(",", "")
    mant = re.split(r"[eE]", t)[0]
    dec = len(mant.split(".")[1]) if "." in mant else 0
    exp = int(re.split(r"[eE]", t)[1]) if re.search(r"[eE]", t) else 0
    return float(t), max(dec - exp, 0), pct


def resolve(obj, path: str):
    cur = obj
    for part in [p for p in path.split("/") if p != ""]:
        if isinstance(cur, list):
            cur = cur[int(part)]
        elif isinstance(cur, dict):
            if part not in cur:
                raise KeyError(part)
            cur = cur[part]
        else:
            raise KeyError(part)
    return cur


class Sources:
    def __init__(self, root: Path, aliases: dict, registry: Path):
        self.root, self.aliases, self.registry = root, aliases, registry
        self._json, self._reg = {}, None

    def reg(self):
        if self._reg is None:
            self._reg = {}
            if self.registry.exists():
                with open(self.registry, newline="", encoding="utf-8") as fh:
                    for row in csv.DictReader(fh):
                        self._reg[row["run_id"]] = row
        return self._reg

    def value(self, source: str, path: str):
        if source.startswith("registry:"):
            row = self.reg().get(source.split(":", 1)[1])
            if row is None:
                raise KeyError(f"run id {source.split(':', 1)[1]} not in registry")
            return float(row[path])
        p = self.root / self.aliases.get(source, source)
        if p not in self._json:
            self._json[p] = json.loads(p.read_text(encoding="utf-8"))
        return resolve(self._json[p], path)


def strip_code(text: str) -> str:
    """Blank out fenced code blocks and inline code (keeps line numbers)."""
    text = re.sub(r"```.*?```", lambda m: re.sub(r"[^\n]", " ", m.group(0)), text, flags=re.S)
    return re.sub(r"`[^`\n]*`", lambda m: " " * len(m.group(0)), text)


def check(report: Path, root: Path = ROOT, registry: Path | None = None, strict=False) -> list[dict]:
    text = report.read_text(encoding="utf-8")
    registry = registry or root / "results" / "registry.csv"
    aliases = {m.group("alias"): m.group("path") for m in SRC.finditer(text)}
    src = Sources(root, aliases, registry)
    rows = []

    def line_of(pos):
        return text.count("\n", 0, pos) + 1

    for m in TAG.finditer(text):
        key, num = m.group("key").strip(), m.group("num")
        row = {"line": line_of(m.start()), "key": key, "printed": num, "ok": False, "msg": ""}
        rows.append(row)
        if num is None:
            row["msg"] = "tag not followed by a number"
            continue
        spec, *flags = key.split("|")
        if "#" not in spec:
            row["msg"] = "key must be SOURCE#PATH"
            continue
        source, path = spec.split("#", 1)
        try:
            v = src.value(source.strip(), path.strip())
        except (KeyError, IndexError, ValueError, FileNotFoundError, json.JSONDecodeError) as e:
            row["msg"] = f"unresolved: {type(e).__name__}: {e}"
            continue
        if v is None or isinstance(v, (dict, list, str, bool)):
            row["msg"] = f"source value is not a number: {v!r}"[:200]
            continue
        v = float(v)
        for fl in flags:
            fl = fl.strip()
            if fl == "pct":
                v *= 100
            elif fl == "abs":
                v = abs(v)
            elif fl == "neg":
                v = -v
            else:
                row["msg"] = f"unknown flag {fl}"
        if row["msg"]:
            continue
        pv, dec, pct = parse_number(num)
        if pct and "pct" not in [f.strip() for f in flags]:
            row["msg"] = "printed % but no |pct flag"
            continue
        tol = 0.5 * 10 ** (-dec) * (1 + 1e-9) + 1e-12
        row["value"] = v
        if not math.isfinite(v):
            row["msg"] = "source value not finite"
        elif abs(pv - v) <= tol:
            row["ok"] = True
        else:
            row["msg"] = f"mismatch: source {v:.{dec + 3}f}"
    reg = src.reg()
    for m in RUN_ID.finditer(text):
        rid = m.group(0)
        if rid not in reg:
            rows.append({"line": line_of(m.start()), "key": f"run id {rid}", "printed": rid, "ok": False,
                         "msg": "run id not in results/registry.csv"})
    if strict:
        clean = TAG.sub(lambda m: " " * len(m.group(0)), strip_code(text))
        clean = SRC.sub(lambda m: " " * len(m.group(0)), clean)
        clean = RUN_ID.sub(lambda m: " " * len(m.group(0)), clean)
        for m in UNTAGGED.finditer(clean):
            rows.append({"line": line_of(m.start()), "key": "(untagged)", "printed": m.group(0), "ok": False,
                         "msg": "untagged number with >= 2 decimals (--strict)"})
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("report", nargs="?", default=str(ROOT / "REPORT.md"))
    ap.add_argument("--strict", action="store_true")
    ap.add_argument("--root", default=str(ROOT))
    a = ap.parse_args(argv)
    report = Path(a.report)
    if not report.exists():
        print(f"{report} not found")
        return 1
    rows = check(report, Path(a.root), strict=a.strict)
    bad = [r for r in rows if not r["ok"]]
    n_tags = sum(1 for r in rows if not r["key"].startswith(("run id", "(untagged)")))
    for r in bad:
        print(f"FAIL line {r['line']}: {r['key']} printed {r['printed']}: {r['msg']}")
    print(f"{n_tags} tagged numbers, {len(bad)} failures")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
