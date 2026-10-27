"""G6 transfer evaluation (scripts/evaluate_transfer.py) on a toy repository with synthetic screens.

Known verdicts: one synthetic background per outcome (transfers / fails / inconclusive), the set-A-first
order, one guarded access per finalist x screen (and none for predictor-only screens, whose replicate
ceiling is frozen), D7 fold-model handling, the D20d no-EC variant, and refusal on a dirty tree.
"""
import json

import numpy as np
import pandas as pd
import pytest

import evaluate_transfer as ETR
import heldout_common as HC
import toy_heldout as TH


def _paths(w):
    return w["paths"].with_(out_dir=w["repo"] / "results" / "transfer",
                            extra={"test_dir": w["repo"] / "results" / "test"})


def _run(w, tw, stage, **kw):
    return ETR.run_g6(w["finalists"], _paths(w), dry=False, stage=stage, threads=2, x_builder=tw["x_builder"],
                      expression=tw["expression"], anchors=tw["anchors"], spec_override=tw["spec"],
                      comparators=tw["comparators"], feature_loader=w["loader"], log=lambda *_: None, **kw)


def _log(w):
    return pd.read_csv(w["paths"].access_log, dtype=str)


@pytest.fixture(scope="module")
def g6(tmp_path_factory):
    mp = pytest.MonkeyPatch()
    mp.setenv("ESSML_UNBLIND", "1")
    w = TH.make_world(tmp_path_factory.mktemp("g6"), n_genes=1200)
    tw = TH.make_transfer_world(w)
    sa = _run(w, tw, "A")
    TH.git(w["repo"], "add", "-A")                 # the PI commits the set-A results before set B
    TH.git(w["repo"], "commit", "-q", "-m", "G6 set A results")
    sb = _run(w, tw, "B")
    out = json.loads((w["repo"] / "results" / "transfer" / "setB" / "summary.json").read_text())
    yield w, tw, sa, sb, out
    mp.undo()


def test_known_verdicts_per_background(g6):
    _, _, _, _, out = g6
    v = out["verdicts"]
    assert v["RPE1"]["Q3a"]["verdict"] == "transfers"
    assert v["hPSC"]["Q3a"]["verdict"] == "fails"
    assert v["bgI"]["Q3a"]["verdict"] == "inconclusive"
    assert "anchor gate" in " ".join(v["bgI"]["Q3a"]["reasons"])
    assert out["overall"]["Q3a"]["verdict"] == "mixed"
    assert out["declared_backgrounds"]["primary_T"]["status"] == "inconclusive by declaration"


def test_ceiling_comparison_uses_frozen_values(g6):
    w, _, _, _, out = g6
    bg = json.loads((w["repo"] / "results" / "transfer" / "setB" / "RPE1.json").read_text())
    s = bg["statistics"]["Q3a_primary"]
    assert s["ceiling"]["source"] == "frozen criteria" and s["ceiling"]["MDE"] == 0.05
    assert s["T"]["point"] >= 0.5 * 0.02 and s["T"]["ci_lo"] > 0
    assert out["ceiling_comparison"]["RPE1"]["Q3a_primary"]["skill_over_ceiling"] > 0.5
    assert "T_ensemble" in s                                    # D7: ensemble is secondary
    assert set(bg["anchor_gate"]["Q3a_primary"]["anchors"]["A_good"]["oof_folds"]) == {0} or \
        len(bg["anchor_gate"]["Q3a_primary"]["anchors"]["A_good"]["oof_folds"]) == 1   # single OOF model


def test_one_access_per_finalist_and_screen(g6):
    w, tw, sa, sb, _ = g6
    log = _log(w)
    g = log[log["granted"] == "True"]
    assert not g.duplicated(["config_hash", "scope"]).any()
    assert (log["granted"] == "True").all()
    fin = HC.load_finalists(w["finalists"], w["repo"])
    planned = ETR.planned_accesses(tw["spec"]["crit"], "all", True, tw["spec"]["h2_arms"])
    assert len(g) == len(planned)
    assert {(fin[r].config_hash, s) for r, s in planned} == set(zip(g["config_hash"], g["scope"]))
    assert not g["scope"].str.endswith(":po").any()             # predictor-only screen never accessed
    assert set(g.loc[g["scope"] == "setA:A1", "config_hash"]) == {fin["T"].config_hash, fin["headline"].config_hash}
    msgs = TH.git(w["repo"], "log", "--format=%s").stdout.splitlines()
    assert sum(m.startswith(("G6 access log", "G4/G6 access log")) for m in msgs) == len(g)


def test_set_a_results_and_descriptive_labels(g6):
    w, _, sa, _, _ = g6
    a = json.loads((w["repo"] / "results" / "transfer" / "setA" / "summary.json").read_text())
    assert a["status"] == "complete" and sorted(a["lines"]) == sorted(w["nonc"])
    for key in ("Q3a_primary", "Q3b_primary"):
        assert set(a["statistics"][key]["per_line"]) == set(w["nonc"])
        assert a["statistics"][key]["descriptive"]["descriptive_verdict"] == "n/a"  # no G4 reference in this toy
    assert "secondary_headline" in a


def test_stage_b_requires_set_a_first(tmp_path, monkeypatch):
    monkeypatch.setenv("ESSML_UNBLIND", "1")
    w = TH.make_world(tmp_path, n_genes=1100, n_dev=60, n_test=20)
    tw = TH.make_transfer_world(w, n_boot=20)
    with pytest.raises(SystemExit, match="setA_done_first"):
        _run(w, tw, "B")
    assert len(_log(w)) == 0


def test_refuses_dirty_tree(tmp_path, monkeypatch):
    monkeypatch.setenv("ESSML_UNBLIND", "1")
    w = TH.make_world(tmp_path, n_genes=1100, n_dev=60, n_test=20)
    tw = TH.make_transfer_world(w, n_boot=20)
    (w["repo"] / "configs" / "seeds.json").write_text('{"split_seed": 3}')
    with pytest.raises(SystemExit, match="tree_clean"):
        _run(w, tw, "A")
    assert len(_log(w)) == 0


def test_no_ec_variant_and_h3(g6):
    w, _, _, _, out = g6
    h3 = out["hypotheses"]["H3"]
    assert h3["noEC_excluded_lines"] == HC.EC_DEV and "result" in h3
    fin = HC.load_finalists(w["finalists"], w["repo"])
    pop = HC.population(w["paths"], dry=False)
    models, _ = HC.load_or_fit_fold_models(fin["H3_noEC"], pop, w["paths"], dry=False, feature_loader=w["loader"],
                                           log=lambda *_: None)
    assert len(models) == 5
    for _, _, meta in models:
        assert not set(HC.EC_DEV) & set(meta["train_ids"])
    tm, _ = HC.load_or_fit_fold_models(fin["T"], pop, w["paths"], dry=False, feature_loader=w["loader"],
                                       log=lambda *_: None)
    assert all(set(HC.EC_DEV) & set(meta["train_ids"]) for _, _, meta in tm)   # T itself trains on them


def test_h1_h2_h4_reported(g6):
    _, _, _, _, out = g6
    hy = out["hypotheses"]
    assert hy["H4"]["result"] == "inconclusive"
    assert hy["H1"]["result"] in ("pass", "fail", "inconclusive") and "H2" in hy
    assert not hy["errors"]


def test_clopper_pearson():
    lo, hi = ETR.clopper_pearson(5, 10)
    assert 0.18 < lo < 0.19 and 0.81 < hi < 0.82
    assert ETR.clopper_pearson(0, 10)[0] == 0.0 and ETR.clopper_pearson(10, 10)[1] == 1.0
