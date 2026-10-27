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
    # set B order: every primary T / no-EC access precedes every secondary (headline, B2) access
    b = g[g["scope"].str.startswith("setB:")].reset_index(drop=True)
    sec = b["config_hash"].isin([fin["headline"].config_hash, fin["B2"].config_hash])
    assert sec.any() and b.index[~sec].max() < b.index[sec].min()
    msgs = TH.git(w["repo"], "log", "--format=%s").stdout.splitlines()
    assert sum(m.startswith(("G6 access log", "G4/G6 access log")) for m in msgs) == len(g)


def test_set_a_results_and_descriptive_labels(g6):
    w, _, sa, _, _ = g6
    a = json.loads((w["repo"] / "results" / "transfer" / "setA" / "summary.json").read_text())
    assert a["status"] == "complete" and sorted(a["lines"]) == sorted(w["nonc"])
    for key in ("Q3a_primary", "Q3b_primary"):
        assert set(a["statistics"][key]["per_line"]) == set(w["nonc"])
        d = a["statistics"][key]["descriptive"]
        assert a["statistics"][key]["test_reference_n"] == len(w["test"])
        assert d["descriptive_verdict"] != "n/a" and np.isfinite(d["median_test_percentile"])
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


def test_planned_accesses_skip_completed_parts(g6):
    _, tw, _, _, _ = g6
    crit, h2 = tw["spec"]["crit"], tw["spec"]["h2_arms"]
    full = ETR.planned_accesses(crit, "B", True, h2)
    part = ETR.planned_accesses(crit, "B", True, h2, done={("B", "RPE1"), ("Bsec", "RPE1")})
    rpe1 = {f"setB:{a}" for a in ETR.verdict_arms({"backgrounds": {"RPE1": crit["backgrounds"]["RPE1"]}})}
    assert rpe1 <= {s for _, s in full} and not rpe1 & {s for _, s in part}
    assert ETR.planned_accesses(crit, "A", True, h2, done={("A",)}) == []


def test_b2_external_profile_fills_unmeasured_genes_with_training_mean():
    from essml.baselines import KNNExpression
    rng = np.random.default_rng(0)
    cols = [f"expression:{g}" for g in range(30)]
    X = pd.DataFrame(rng.normal(size=(40, 30)), index=[f"L{i}" for i in range(40)], columns=cols)
    Y = pd.DataFrame(rng.normal(size=(40, 5)), index=X.index, columns=list("abcde"))
    m = KNNExpression(k=5, n_top=10).fit(Y, X=X)
    f = HC.Finalist("B2", "b2", "x", "0" * 64, {"model": "B2"})
    row = pd.Series(rng.normal(size=30), index=[str(g) for g in range(30)])
    row[[m.feat_[0].split(":")[1], m.feat_[1].split(":")[1]]] = np.nan      # unmeasured in the external profile
    pid, Xe = ETR.real_x_builder(f, m, "state:x", row, "WT")
    assert np.isfinite(Xe[list(m.feat_)].to_numpy()).all()
    assert Xe.loc[pid, m.feat_[0]] == pytest.approx(m.xm_[0])
    assert np.isfinite(m.predict([pid], X=Xe).to_numpy()).all()
    with pytest.raises(NotImplementedError):
        ETR.real_x_builder(HC.Finalist("M5", "m5", "x", "1" * 64, {"model": "M5"}), m, "state:x", row, "WT")


def test_resume_after_crash_inside_a_part_marks_it_spent(g6):
    """Simulate a crash after hPSC's screens were accessed but before its results were written: on resume the
    part is marked access_spent_incomplete (no second access); completed parts are not re-accessed."""
    w, tw, _, _, _ = g6
    d = w["repo"] / "results" / "transfer" / "setB"
    for f in ("hPSC.json", "hPSC_secondary.json", "bgI_secondary.json"):
        (d / f).unlink()
    TH.git(w["repo"], "add", "-A")                   # the PI commits what the crashed run left (as instructed)
    TH.git(w["repo"], "commit", "-q", "-m", "partial G6 set B")
    n0 = len(_log(w))
    s = _run(w, tw, "B", resume=True)
    assert len(_log(w)) == n0                                   # nothing accessed again
    assert json.loads((d / "hPSC.json").read_text())["status"] == "access_spent_incomplete"
    assert s["setB"]["verdicts"]["hPSC"]["Q3a"]["verdict"] == "inconclusive"
    assert s["setB"]["verdicts"]["RPE1"]["Q3a"]["verdict"] == "transfers"   # reloaded, not recomputed
    assert sorted(map(tuple, s["preflight"]["spent_incomplete_parts"]["parts"])) == [
        ("B", "hPSC"), ("Bsec", "bgI"), ("Bsec", "hPSC")]
