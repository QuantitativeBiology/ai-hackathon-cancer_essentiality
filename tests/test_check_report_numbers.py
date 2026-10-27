"""scripts/check_report_numbers.py on toy report + results files."""
import json

import pytest

import check_report_numbers as C

REG = "run_id,status,P1_variable\n20261026-222117_m1_E_resB0_91e52f9d,ok,0.30812\n"


@pytest.fixture
def world(tmp_path):
    (tmp_path / "results" / "test").mkdir(parents=True)
    (tmp_path / "results" / "test" / "summary.json").write_text(json.dumps(
        {"table": [{"P1_variable": 0.30812, "dP1_vs_B1_variable": 0.10849, "frac": 0.9437, "n": 4362,
                    "n_r_gt_0.3_variable": 1121, "neg": -0.0196, "missing": None}]}))
    (tmp_path / "results" / "registry.csv").write_text(REG)
    return tmp_path


def _run(world, body, strict=False):
    rep = world / "REPORT.md"
    rep.write_text("<!--n-src:test=results/test/summary.json-->\n" + body, encoding="utf-8")
    return C.check(rep, world, strict=strict)


def test_all_tags_match(world):
    body = ("P1 <!--n:test#table/0/P1_variable-->0.308 (Δ <!--n:test#table/0/dP1_vs_B1_variable-->+0.108), "
            "<!--n:test#table/0/frac|pct-->94.4% of genes, <!--n:test#table/0/n-->4,362 genes, "
            "<!--n:test#table/0/n_r_gt_0.3_variable-->1121 with r > 0.3, R4 <!--n:test#table/0/neg-->−0.020, "
            "registry <!--n:registry:20261026-222117_m1_E_resB0_91e52f9d#P1_variable-->0.3081 "
            "(run 20261026-222117_m1_E_resB0_91e52f9d), full path <!--n:results/test/summary.json#table/0/P1_variable-->0.31")
    rows = _run(world, body)
    assert rows and all(r["ok"] for r in rows), [r for r in rows if not r["ok"]]


@pytest.mark.parametrize("body, msg", [
    ("<!--n:test#table/0/P1_variable-->0.309", "mismatch"),
    ("<!--n:test#table/0/P1_variable-->0.3082", "mismatch"),           # more decimals than the value supports
    ("<!--n:test#table/0/nope-->0.3", "unresolved"),
    ("<!--n:test#table/5/P1_variable-->0.3", "unresolved"),
    ("<!--n:test#table/0/missing-->0.3", "not a number"),
    ("<!--n:test#table/0/frac-->94.4%", "pct"),
    ("<!--n:test#table/0/P1_variable--> about", "not followed"),
    ("<!--n:testtable/0/P1_variable-->0.3", "SOURCE#PATH"),
    ("<!--n:test#table/0/dP1_vs_B1_variable-->-0.108", "mismatch"),    # sign checked
    ("see run 20261027-000000_ghost_deadbeef", "not in results/registry.csv"),
])
def test_failures_detected(world, body, msg):
    rows = _run(world, body)
    bad = [r for r in rows if not r["ok"]]
    assert bad and msg in bad[0]["msg"]


def test_rounding_boundary(world):
    assert all(r["ok"] for r in _run(world, "<!--n:test#table/0/P1_variable-->0.3081"))   # 0.30812 -> 0.3081
    assert all(r["ok"] for r in _run(world, "<!--n:test#table/0/frac-->0.94"))
    assert not all(r["ok"] for r in _run(world, "<!--n:test#table/0/frac-->0.95"))


def test_strict_flags_untagged_numbers_outside_code(world):
    body = "tagged <!--n:test#table/0/P1_variable-->0.308, untagged 0.276, code `0.123` and\n```\n0.555\n```\n"
    rows = _run(world, body, strict=True)
    bad = [r for r in rows if not r["ok"]]
    assert [r["printed"] for r in bad] == ["0.276"]
    assert all(r["ok"] for r in _run(world, body.replace(", untagged 0.276", ""), strict=True))


def test_main_exit_codes(world, capsys):
    rep = world / "REPORT.md"
    rep.write_text("<!--n:results/test/summary.json#table/0/P1_variable-->0.308", encoding="utf-8")
    assert C.main([str(rep), "--root", str(world)]) == 0
    rep.write_text("<!--n:results/test/summary.json#table/0/P1_variable-->0.400", encoding="utf-8")
    assert C.main([str(rep), "--root", str(world)]) == 1
    assert "FAIL" in capsys.readouterr().out
