"""G4 test-set evaluation (scripts/evaluate_test.py, essml.heldout) on toy data only.

Known answers for the resampling layer, positive/negative controls, and the blinding protocol on a toy
git repository: one guard access per finalist, refusal on a dirty tree, refusal of a second run, and
recovery after a crash that happens after an access (no second access).
"""
import json

import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import average_precision_score

import evaluate_test as ET
import heldout_common as HC
import toy_heldout as TH
from essml import guard, heldout
from essml import metrics as M
from essml.evaluate import evaluate


# ============================================================================ resampling layer == frozen metrics
def _dup(A, w):
    return np.repeat(A, w.astype(int), axis=0)


def test_weighted_pearson_equals_frozen_on_duplicated_rows():
    rng = np.random.default_rng(1)
    n, G = 40, 30
    T = rng.normal(size=(n, G))
    T[rng.random((n, G)) < 0.1] = np.nan
    P = T + rng.normal(0, 1, (n, G))
    P[:, 0] = 3.0                                  # constant prediction -> r = 0
    P = np.where(np.isnan(P), 0.0, P)
    W = rng.integers(0, 3, (5, n)).astype(float)
    W[:, :3] = 1
    R = heldout.weighted_pearson_cols(P, T, W, min_obs=2)
    for b in range(W.shape[0]):
        ref = M.per_gene_pearson(_dup(P, W[b]), _dup(T, W[b]), min_obs=2)
        assert np.allclose(R[b], ref, atol=1e-10, equal_nan=True)
    assert np.all(R[:, 0] == 0.0)
    cm = heldout.weighted_pearson_class_means(P, T, W, {"a": np.ones(G, bool)}, min_obs=2)
    assert np.allclose(cm["a"], np.nanmean(R, 1))


def test_weighted_ap_equals_sklearn_on_duplicated_rows():
    rng = np.random.default_rng(2)
    n, G = 50, 8
    T = rng.normal(-0.4, 0.6, (n, G))
    T[rng.random((n, G)) < 0.1] = np.nan
    P = np.round(T + rng.normal(0, 0.5, (n, G)), 1)   # ties
    P = np.where(np.isnan(P), 0.0, P)
    W = rng.integers(0, 3, (4, n)).astype(float)
    ap, prev = heldout.weighted_ap_cols(P, T, W, np.ones(G, bool))
    for b in range(4):
        for g in range(G):
            o = ~np.isnan(T[:, g])
            y = _dup((T[o, g] < -0.5)[:, None], W[b, o])[:, 0]
            s = _dup((-P[o, g])[:, None], W[b, o])[:, 0]
            if y.any():
                assert ap[b, g] == pytest.approx(average_precision_score(y, s), abs=1e-12)
                assert prev[b, g] == pytest.approx(y.mean())


def test_patient_weights_stratified_and_reproducible():
    lin = pd.Series(["a"] * 6 + ["b"] * 4, index=[f"l{i}" for i in range(10)])
    pat = pd.Series(["p0", "p0", "p1", "p2", "p3", "p4", "q0", "q1", "q1", "q2"], index=lin.index)
    W1 = heldout.patient_bootstrap_weights(lin, pat, 200, 7)
    W2 = heldout.patient_bootstrap_weights(lin, pat, 200, 7)
    assert np.array_equal(W1, W2)
    assert np.array_equal(W1[:, 0], W1[:, 1]) and np.array_equal(W1[:, 7], W1[:, 8])  # same patient
    # stratum a has 5 patients, b has 3: every replicate draws exactly that many in each stratum
    assert np.all(W1[:, [0, 2, 3, 4, 5]].sum(1) == 5) and np.all(W1[:, [6, 7, 9]].sum(1) == 3)
    assert not np.array_equal(W1, heldout.patient_bootstrap_weights(lin, pat, 200, 8))


# ============================================================================ positive / negative controls
@pytest.fixture(scope="module")
def synth():
    """Lines x genes with shared line factors (genes correlated, as in DepMap), lineages and patient pairs."""
    rng = np.random.default_rng(5)
    n, G = 120, 200
    ids = [f"L{i}" for i in range(n)]
    genes = [str(g) for g in range(G)]
    F, L = rng.normal(size=(n, 5)), rng.normal(0, 0.5, (5, G))
    T = pd.DataFrame(-0.3 + F @ L + rng.normal(0, 0.6, (n, G)), index=ids, columns=genes)
    T.iloc[:, :30] -= 1.5
    gc = pd.DataFrame({"all": True, "variable": True, "selective": np.arange(G) >= 30,
                       "rare_selective": False, "common_essential": np.arange(G) < 30, "nonessential": False},
                      index=genes)
    lin = pd.Series(rng.choice(["a", "b", "c"], n), index=ids)
    pat = pd.Series([f"p{i // 2}" for i in range(n)], index=ids)  # pairs of lines share a patient
    W = heldout.patient_bootstrap_weights(lin, pat, 200, 20261026)
    b0 = pd.DataFrame(np.broadcast_to(T.mean().to_numpy(), T.shape), index=ids, columns=genes)
    return T, gc, W, b0


def _boot(pred, T, gc, W, b0):
    ev = evaluate(pred, T, gc, None, baseline_preds={"B0": b0}, min_obs=50, spearman=False)
    return ev, heldout.bootstrap_test_metrics(pred, T, gc, W, ev, baselines={"B0": b0}, min_obs=50)


def test_perfect_predictor_p1_one_and_ci_contains_truth(synth):
    T, gc, W, b0 = synth
    ev, b = _boot(T.copy(), T, gc, W, b0)
    assert ev["summary"]["variable"]["P1"] == pytest.approx(1.0)
    assert b["P1_variable"]["ci_lo"] == pytest.approx(1.0) and b["P1_variable"]["ci_hi"] == pytest.approx(1.0)
    assert b["P2b_nap_prev05"]["ci_lo"] <= 1.0 + 1e-12 <= b["P2b_nap_prev05"]["ci_hi"] + 2e-12
    assert b["P2b_nap_prev05"]["point"] == pytest.approx(1.0)


def test_constant_predictor_p1_zero(synth):
    T, gc, W, b0 = synth
    ev, b = _boot(b0.copy(), T, gc, W, b0)
    assert ev["summary"]["variable"]["P1"] == 0.0
    assert b["P1_variable"]["ci_lo"] == 0.0 and b["P1_variable"]["ci_hi"] == 0.0
    assert b["dP1_vs_B0_variable"]["point"] == 0.0


def test_shuffled_predictor_ci_includes_zero(synth):
    """Row-shuffled predictions have no skill: the CI covers 0 (in >= 75% of 20 shuffles; nominal 95%,
    the bootstrap conditions on one shuffle), while the unshuffled predictor is clearly positive."""
    T, gc, W, b0 = synth
    Tm = T.to_numpy()
    allg = {"v": np.ones(T.shape[1], bool)}
    cover = []
    for s in range(20):
        rng = np.random.default_rng(100 + s)
        P = (T + rng.normal(0, 0.5, T.shape)).to_numpy()[rng.permutation(len(T))]
        boot = heldout.weighted_pearson_class_means(P, Tm, W, allg, min_obs=50)["v"]
        b = heldout.summarise(np.nanmean(M.per_gene_pearson(P, Tm, 50)), boot)
        cover.append(b["ci_lo"] < 0 < b["ci_hi"])
    assert cover[0] and np.mean(cover) >= 0.75
    _, bt = _boot(T + np.random.default_rng(1).normal(0, 0.5, T.shape), T, gc, W, b0)
    assert bt["P1_variable"]["ci_lo"] > 0.3


def test_unit_weights_reproduce_frozen_point_estimates(synth):
    T, gc, _, b0 = synth
    P = T + np.random.default_rng(3).normal(0, 0.8, T.shape)
    ev = evaluate(P, T, gc, None, baseline_preds={"B0": b0}, min_obs=50, spearman=False)
    b = heldout.bootstrap_test_metrics(P, T, gc, np.ones((2, len(T))), ev, baselines={"B0": b0}, min_obs=50)
    for k, c in (("P1_variable", "variable"), ("P1_selective", "selective")):
        assert b[k]["ci_lo"] == pytest.approx(ev["summary"][c]["P1"], abs=1e-12)
    assert b["P2a"]["ci_lo"] == pytest.approx(ev["summary"]["all"]["P2a"], abs=1e-12)
    assert b["P2b_nap_prev05"]["ci_lo"] == pytest.approx(ev["summary"]["selective"]["P2b_nap_prev05"], abs=1e-12)
    assert b["dP2a_vs_B0"]["ci_lo"] == pytest.approx(ev["summary"]["all"]["dP2a_vs_B0"], abs=1e-12)


# ============================================================================ protocol on a toy repository
def _run(w, **kw):
    return ET.run_g4(w["finalists"], w["paths"], dry=kw.pop("dry", False), threads=2, n_boot=kw.pop("n_boot", 30),
                     isolate=False, feature_loader=w["loader"], log=lambda *_: None, min_obs=40, **kw)


def _log(w):
    return pd.read_csv(w["paths"].access_log, dtype=str)


@pytest.fixture(scope="module")
def g4(tmp_path_factory):
    mp = pytest.MonkeyPatch()
    mp.setenv("ESSML_UNBLIND", "1")
    w = TH.make_world(tmp_path_factory.mktemp("g4"), n_genes=400)
    s = _run(w)
    yield w, s
    mp.undo()


def test_one_access_per_finalist(g4):
    w, s = g4
    log = _log(w)
    fin = HC.load_finalists(w["finalists"], w["repo"])
    assert len(log) == 6 and (log["granted"] == "True").all() and (log["scope"] == "test").all()
    assert sorted(log["config_hash"]) == sorted(f.config_hash for f in fin.values())
    assert s["n_accesses_granted"] == 6
    assert all(r["status"] == "complete" for r in s["table"])


def test_access_log_committed_after_each_access_and_tree_clean(g4):
    w, _ = g4
    msgs = TH.git(w["repo"], "log", "--format=%s").stdout.splitlines()
    assert sum(m.startswith("G4/G6 access log") for m in msgs) == 6
    st = TH.git(w["repo"], "status", "--porcelain", "--untracked-files=no").stdout.strip()
    assert st == ""                                      # outputs are untracked; nothing tracked was modified


def test_predictions_saved_before_access_and_targets_not_on_disk(g4):
    w, _ = g4
    man = json.loads((w["paths"].out_dir / "predictions_manifest.json").read_text())
    first_access = _log(w)["timestamp"].min()
    assert all(e["timestamp"] <= first_access for e in man.values())
    truth = w["Y"].loc[w["test"]]
    for p in list(w["repo"].rglob("*.parquet")):
        if p.name == "crispr.parquet":
            continue
        d = pd.read_parquet(p)
        if d.shape == truth.shape:
            assert not np.allclose(d.to_numpy(np.float64), truth.to_numpy(np.float64), equal_nan=True), p


def test_toy_numbers_sensible(g4):
    _, s = g4
    t = {r["role"]: r for r in s["table"]}
    assert t["B0"]["P1_variable"] == 0.0
    assert t["headline"]["P1_variable"] > t["B1"]["P1_variable"] > 0
    lo, hi = t["headline"]["dP1_vs_B1_variable_ci"]
    assert 0 < lo <= t["headline"]["dP1_vs_B1_variable"] <= hi
    assert s["library_groups_reproduced"] is True


def test_per_tier_reporting_and_tier_table_hash(g4):
    """D22: P1 / P2b per tier, Spearman(dev mean_r, held-out r) with a CI, sha256 of the tier table everywhere."""
    w, s = g4
    sha = HC.sha256_file(w["paths"].predictable_genes)
    assert s["predictable_genes"]["sha256"] == sha and s["preflight"]["predictable_genes_committed"]["pass"]
    head = json.loads((w["paths"].out_dir / "toy_m1" / "summary.json").read_text())
    assert head["predictable_genes"]["sha256"] == sha
    pt = head["per_tier"]
    assert pt["predictable_genes"]["sha256"] == sha
    assert set(pt["tiers"]) == {"high", "medium", "low"}
    hi, lo = pt["tiers"]["high"], pt["tiers"]["low"]
    assert hi["n_genes"] > 0 and lo["n_genes"] > 0
    assert hi["P1"] > lo["P1"]                       # planted-signal genes are the toy 'high' tier
    assert hi["P1_ci"]["ci_lo"] <= hi["P1"] <= hi["P1_ci"]["ci_hi"]
    # per-tier point estimates equal the frozen per-gene output restricted to the tier
    pg = pd.read_csv(w["paths"].out_dir / "toy_m1" / "per_gene.csv.gz", dtype={"gene": str})
    m = (pg["tier"] == "high") & (pg["n_obs"] >= 40)
    assert np.isclose(pg.loc[m, "r"].mean(), hi["P1"], atol=1e-4)
    rho = pt["spearman_dev_mean_r_vs_heldout_r"]["evaluable"]
    assert rho["point"] > 0.3 and rho["ci_lo"] <= rho["point"] <= rho["ci_hi"]
    assert rho["ci_method"] == "normal" and "ci_percentile" in rho                # duplicate-line bias (as P2b)
    assert np.isclose(rho["point"] - rho["ci_lo"], rho["ci_hi"] - rho["point"])
    row = {r["role"]: r for r in s["table"]}["headline"]
    assert row["P1_high"] == hi["P1"] and row["spearman_dev_mean_r_vs_r"] == rho["point"]


def test_outputs_carry_confidence_and_tier(g4):
    w, _ = g4
    pg = pd.read_csv(w["paths"].out_dir / "toy_m1" / "per_gene.csv.gz", dtype={"gene": str})
    assert {"confidence", "tier", "mean_r_dev_cv"} <= set(pg.columns)
    ref = pd.read_csv(w["paths"].predictable_genes, dtype={"entrez": str}).set_index("entrez")
    assert (pg["tier"].to_numpy() == ref["tier"].reindex(pg["gene"]).fillna("not_evaluable").to_numpy()).all()
    man = json.loads((w["paths"].out_dir / "predictions_manifest.json").read_text())
    for e in man.values():
        tidy = pd.read_parquet(w["repo"] / e["predictions_with_confidence"]["file"])
        assert {"ModelID", "entrez", "predicted_effect", "confidence", "tier"} <= set(tidy.columns)
        assert e["predictions_with_confidence"]["predictable_genes_sha256"] == HC.sha256_file(w["paths"].predictable_genes)


def _failed_checks(exc):
    """The list of failed preflight checks on the first line of the refusal message."""
    first = str(exc.value).splitlines()[0]
    return first


def test_refuses_without_committed_tier_table(tmp_path, monkeypatch):
    monkeypatch.setenv("ESSML_UNBLIND", "1")
    w = TH.make_world(tmp_path, n_genes=200, n_dev=60, n_test=20)
    (w["repo"] / ".gitignore").write_text("\n".join(["cache/", "*.parquet", "*.pkl", "configs/predictable_genes.csv", ""]))
    TH.git(w["repo"], "rm", "-q", "--cached", "configs/predictable_genes.csv")
    TH.git(w["repo"], "commit", "-q", "-am", "untrack the tier table")          # tree clean, table untracked
    with pytest.raises(SystemExit) as e:
        _run(w)
    assert _failed_checks(e) == "G4 refused to start: ['predictable_genes_committed']"   # tree_clean passed
    assert len(_log(w)) == 0 and not (w["paths"].out_dir / "predictions_manifest.json").exists()


def test_refuses_modified_tier_table(tmp_path, monkeypatch):
    monkeypatch.setenv("ESSML_UNBLIND", "1")
    w = TH.make_world(tmp_path, n_genes=200, n_dev=60, n_test=20)
    p = w["paths"].predictable_genes
    p.write_text(p.read_text().replace(",high,", ",low,"))
    with pytest.raises(SystemExit) as e:
        _run(w)
    assert "predictable_genes_committed" in _failed_checks(e)
    assert len(_log(w)) == 0


def test_predictable_check_unit(tmp_path):
    w = TH.make_world(tmp_path, n_genes=200, n_dev=60, n_test=20)
    paths, p = w["paths"], w["paths"].predictable_genes
    ok = ET.predictable_check(paths)
    assert ok["pass"] and ok["sha256"] == HC.sha256_file(p) and ok["problems"] == []
    assert not ET.predictable_check(paths, headline_cv_run_id="other_run")["pass"]      # wrong source run
    assert ET.predictable_check(paths, headline_cv_run_id="toy_cv")["pass"]
    assert not ET.predictable_check(paths, dry=True)["pass"]                             # primary table in a dry run
    t = pd.read_csv(p)
    for col, val in (("folds_used", "1;2;3;4"), ("table_role", "secondary (m1)")):     # wrong table, committed
        t2 = t.copy()
        t2[col] = val
        t2.to_csv(p, index=False)
        TH.git(w["repo"], "commit", "-q", "-am", f"{col}")
        assert not ET.predictable_check(paths)["pass"]
    t2 = t.copy()
    t2["folds_used"], t2["table_role"] = "1;2;3;4", "dry_run"
    t2.to_csv(p, index=False)
    TH.git(w["repo"], "commit", "-q", "-am", "dry table")
    assert ET.predictable_check(paths, dry=True)["pass"] and not ET.predictable_check(paths)["pass"]
    p.write_text("\n".join(["not,a,table", "1,2,3", ""]))                                               # unreadable + modified
    r = ET.predictable_check(paths)
    assert not r["pass"] and any("unreadable" in x for x in r["problems"])
    TH.git(w["repo"], "commit", "-q", "-am", "garbage")
    assert not ET.predictable_check(paths)["pass"]                                       # committed but unreadable
    p.unlink()
    assert not ET.predictable_check(paths)["pass"]


def test_dry_run_reads_the_folds_1_to_4_tier_table():
    t = pd.read_csv(ET.DRY_TIER, dtype=str)
    assert set(t["folds_used"]) == {"1;2;3;4"} and set(t["table_role"]) == {"dry_run"}
    assert "r_fold0" not in t.columns


def test_second_run_refused_without_new_access(g4):
    w, _ = g4
    n0 = len(_log(w))
    with pytest.raises(SystemExit, match="guard_would_grant"):
        _run(w)
    assert len(_log(w)) == n0


def test_resume_after_complete_run_accesses_nothing(g4):
    w, _ = g4
    TH.git(w["repo"], "add", "-A")                       # the PI commits the G4 outputs
    TH.git(w["repo"], "commit", "-q", "-m", "G4 results")
    n0 = len(_log(w))
    s = _run(w, resume=True)
    assert len(_log(w)) == n0 and all(r["status"] == "complete" for r in s["table"])


def test_refuses_dirty_tree(tmp_path, monkeypatch):
    monkeypatch.setenv("ESSML_UNBLIND", "1")
    w = TH.make_world(tmp_path, n_genes=200, n_dev=60, n_test=20)
    (w["repo"] / "configs" / "seeds.json").write_text('{"split_seed": 1}')        # modified tracked file
    with pytest.raises(SystemExit, match="tree_clean"):
        _run(w)
    TH.git(w["repo"], "checkout", "--", "configs/seeds.json")
    (w["repo"] / "stray.txt").write_text("x")                                   # untracked file
    with pytest.raises(SystemExit, match="tree_clean"):
        _run(w)
    assert len(_log(w)) == 0 and not (w["paths"].out_dir / "predictions_manifest.json").exists()


def test_refuses_without_unblind(tmp_path, monkeypatch):
    monkeypatch.delenv("ESSML_UNBLIND", raising=False)
    w = TH.make_world(tmp_path, n_genes=200, n_dev=60, n_test=20)
    with pytest.raises(SystemExit, match="guard_would_grant"):
        _run(w)
    assert len(_log(w)) == 0


def test_crash_after_access_needs_no_second_access(tmp_path, monkeypatch):
    monkeypatch.setenv("ESSML_UNBLIND", "1")
    w = TH.make_world(tmp_path, n_genes=300)

    def boom_core(ctx, f, pred, truth, refs):
        if f.role == "T":
            raise RuntimeError("simulated crash after access")
        return ET.step_core(ctx, f, pred, truth, refs)

    def boom_late(ctx, f, pred, truth, refs):
        raise RuntimeError("late step fails")

    steps = [("core", boom_core)] + ET.STEPS[1:] + [("late", boom_late)]
    s = _run(w, steps=steps)
    st = {r["role"]: r["status"] for r in s["table"]}
    assert st["T"] == "core_failed" and st["headline"] == "complete_with_errors"
    hl = json.loads((w["paths"].out_dir / "toy_m1" / "summary.json").read_text())
    assert "core" in hl and "bootstrap" in hl and "late" in hl["errors"]   # earlier results kept
    assert (_log(w)["granted"] == "True").sum() == 6
    # resume (after the PI commits the outputs): predictions reused, complete finalists skipped, T's spent
    # access is not retried and does not block the others
    TH.git(w["repo"], "add", "-A")
    TH.git(w["repo"], "commit", "-q", "-m", "partial G4 results")
    s2 = _run(w, resume=True)
    st2 = {r["role"]: r["status"] for r in s2["table"]}
    assert st2["T"] == "access_spent_incomplete" and st2["headline"] == "complete_with_errors"
    assert s2["preflight"]["spent_incomplete"]["finalists"] == ["toy_t"]
    assert (_log(w)["granted"] == "True").sum() == 6 and len(_log(w)) == 6


def test_ledger_refuses_second_request_before_guard(tmp_path, monkeypatch):
    monkeypatch.setenv("ESSML_UNBLIND", "1")
    w = TH.make_world(tmp_path, n_genes=200, n_dev=60, n_test=20)
    f = HC.load_finalists(w["finalists"], w["repo"])["B0"]
    led = HC.AccessLedger()
    HC.guarded_targets(f, w["test"], w["paths"], led, dry=False)
    n = len(_log(w))
    with pytest.raises(guard.BlindingError, match="ledger"):
        HC.guarded_targets(f, w["test"], w["paths"], led, dry=False)
    assert len(_log(w)) == n                       # the guard never saw the second request


def test_dry_run_never_shows_non_dev_ids_to_guard(tmp_path):
    w = TH.make_world(tmp_path, n_genes=200, n_dev=60, n_test=20)
    f = HC.load_finalists(w["finalists"], w["repo"])["B0"]
    with pytest.raises(guard.BlindingError, match="dry run requested non-dev"):
        HC.guarded_targets(f, w["dev"][:3] + w["test"][:1], w["paths"], HC.AccessLedger(), dry=True)
    assert len(_log(w)) == 0
    Y = HC.guarded_targets(f, w["dev"][:3], w["paths"], HC.AccessLedger(), dry=True, scope_label="x")
    assert Y.shape[0] == 3 and len(_log(w)) == 0


def test_finalists_schema_checks(tmp_path):
    w = TH.make_world(tmp_path, n_genes=200, n_dev=60, n_test=20)
    d = json.loads(w["finalists"].read_text())
    d["T"]["config_hash"] = "0" * 64
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps(d))
    with pytest.raises(ValueError, match="config hash"):
        HC.load_finalists(bad, w["repo"])
    fin = HC.load_finalists(w["finalists"], w["repo"])
    HC.check_noec(fin["H3_noEC"])
    fin["H3_noEC"].cfg["exclude_from_training"] = HC.EC_DEV[:3]
    with pytest.raises(ValueError, match="D20d"):
        HC.check_noec(fin["H3_noEC"])


def test_example_finalists_file_matches_configs():
    from pathlib import Path
    root = Path(__file__).resolve().parents[1]
    fin = HC.load_finalists(root / "configs" / "finalists.example.json", root)
    assert fin["headline"].name == "m1_E_resB0" and fin["T"].name == "t_E_status"
    assert HC.check_noec_matches_T(fin["H3_noEC"], fin["T"])


def test_p2b_ci_is_centred_on_point_despite_duplicate_bias(synth):
    """Resampled AP is biased upward (duplicated positives credit their copies); the P2b CI is point +/- z SE."""
    T, gc, W, b0 = synth
    P = T + np.random.default_rng(4).normal(0, 1.0, T.shape)
    _, b = _boot(P, T, gc, W, b0)
    r = b["P2b_nap_prev05"]
    assert r["ci_method"] == "normal" and r["ci_lo"] < r["point"] < r["ci_hi"]
    assert r["ci_hi"] - r["point"] == pytest.approx(r["point"] - r["ci_lo"])
    assert r["bias"] > 0                                    # the artefact the normal interval avoids
    assert not b["P1_variable"]["point_outside_percentile_ci"]
