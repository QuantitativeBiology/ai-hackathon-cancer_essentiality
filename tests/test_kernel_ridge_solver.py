"""Solver checks moved from the independent G2 review (reviews/g2_review/check_solver.py, part 1;
review findings 3-4) plus G3 additions. Fast (seconds).

1. p > n (the production regime: ~18k features, ~800 rows) with merged NaN patterns (intersection
   rows), an unmerged near-duplicate, a fallback pattern and patient groups that straddle NaN-pattern
   boundaries: closed-form leave-group-out r equals brute force for the B0 residual AND for a B1
   leave-line-out residual offset (R != Y, fixed targets); predictions at new points equal primal
   ridge.
2. Real path (X=None: raw blocks from cache/features through FoldTransformer with role checks):
   perturbing the feature rows of every non-training line leaves the fitted model and its
   predictions on training rows unchanged (M1 B0 residual; M1b with tuned block weights).
3. M1b block-weight tuning: <= 25 criterion evaluations, starts at the anchor-only model, gives a
   pure-noise block weight 0 and a signal block a positive weight.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from essml.baselines import leave_line_out_lineage_residual
from essml.models import kernel_ridge as KR
from essml.models.kernel_ridge import fit_pattern_ridge, predict_pattern_ridge, tune_block_weights

ROOT = Path(__file__).resolve().parents[1]
FEAT = ROOT / "cache" / "features"


# ---------------------------------------------------------------------------- 1. p > n brute force
def _ridge_dual(Z, y, lam):
    zm, ym = Z.mean(0), y.mean()
    Zc = Z - zm
    a = np.linalg.solve(Zc @ Zc.T + lam * np.eye(len(y)), y - ym)
    beta = Zc.T @ a
    return lambda Zn: ym + (Zn - zm) @ beta


def _ridge_primal(Z, y, lam):
    zm, ym = Z.mean(0), y.mean()
    Zc = Z - zm
    beta = np.linalg.solve(Zc.T @ Zc + lam * np.eye(Z.shape[1]), Zc.T @ (y - ym))
    return lambda Zn: ym + (Zn - zm) @ beta


def _make(seed=3, n=90, p=150, G=10):
    """reviews/g2_review/check_solver.py ``make`` (verbatim design)."""
    rng = np.random.default_rng(seed)
    Z = rng.normal(size=(n, p)) * rng.uniform(0.3, 2, size=p)
    lin = rng.choice(list("abcd"), size=n, p=[0.4, 0.3, 0.2, 0.1])
    leff = {k: rng.normal(size=G) for k in "abcd"}
    Y = Z[:, :20] @ rng.normal(size=(20, G)) * 0.15 + np.array([leff[l] for l in lin]) \
        + rng.normal(size=(n, G))
    Y[:, 3] += 40.0
    Y[:30] += 25.0  # NaN-region rows hold very different values (centring bugs would show)
    obs = np.ones((n, G), bool)
    obs[:30, 4] = obs[:30, 5] = False   # pattern A
    obs[:29, 8] = False                 # A's rows + row 29 -> merges into A (tol = 1)
    obs[50, 6] = False                  # full minus one row: merge guard keeps it separate
    obs[60:63, 7] = False               # differs by 3 -> own cluster
    obs[30:, 9] = False                 # 30 rows -> fallback (min_rows 40)
    groups = np.arange(n).astype(object)
    for a, b in [(28, 31), (29, 30), (59, 61), (62, 70), (10, 11), (12, 13)]:
        groups[b] = groups[a]           # patient groups straddling the NaN boundaries
    groups[14] = groups[12]
    return Z, np.where(obs, Y, np.nan), lin, groups.astype(str)


@pytest.mark.parametrize("resid", ["B0", "B1"])
def test_p_gt_n_group_loo_and_prediction_match_brute_force(resid):
    Z, Y, lin, groups = _make()
    n, G = Y.shape
    grid = np.array([1.0, 30.0, 300.0, 3000.0])
    ids = pd.Index([f"L{i}" for i in range(n)])
    if resid == "B1":  # leave-line-out lineage residual: R != Y
        R = leave_line_out_lineage_residual(pd.DataFrame(Y, index=ids), pd.Series(lin, index=ids),
                                            5.0).to_numpy()
        assert not np.allclose(np.nan_to_num(R), np.nan_to_num(Y - np.nanmean(Y, 0)))
    else:
        R = Y - np.nanmean(Y, 0)
    sol = fit_pattern_ridge(Z @ Z.T, R, Y, np.zeros(G, int), groups=groups, lam_grid=grid,
                            merge_frac=0.012, min_rows=40)
    rep = sol["report"]
    assert rep["n_patterns"] == 6 and rep["n_patterns_merged"] == 1 and rep["n_fallback_genes"] == 1
    cl_of = {int(g): c for c in sol["clusters"] for g in c["genes"]}
    Znew = np.random.default_rng(9).normal(size=(5, Z.shape[1]))
    Pnew = predict_pattern_ridge(Znew @ Z.T, sol["clusters"], G)
    for g in range(G):
        if g not in cl_of:
            assert np.all(Pnew[:, g] == 0)  # fallback: offset only
            continue
        rows = cl_of[g]["rows"]
        f = _ridge_primal(Z[rows], R[rows, g], sol["lam"][g])
        assert np.allclose(f(Znew), Pnew[:, g], atol=1e-8)
        gr = groups[rows]
        for li, lam in enumerate(grid):
            p_fix = np.empty(len(rows))
            for gid in np.unique(gr):
                te = gr == gid
                trr = rows[~te]
                p_fix[te] = (Y[rows[te], g] - R[rows[te], g]) + _ridge_dual(Z[trr], R[trr, g], lam)(Z[rows[te]])
            assert abs(sol["r_grid"][li, g] - np.corrcoef(p_fix, Y[rows, g])[0, 1]) < 1e-8


# ---------------------------------------------------------------------------- 3. tuning (toy)
def test_block_weight_tuning_budget_start_and_direction():
    rng = np.random.default_rng(4)
    n, G = 120, 30
    Zs = rng.normal(size=(n, 40))      # signal block
    Zn = rng.normal(size=(n, 400))     # pure-noise block (large, would dominate at weight 1)
    Ze = rng.normal(size=(n, 200))     # anchor block with weak signal
    Y = Zs[:, :5] @ rng.normal(size=(5, G)) + 0.3 * Ze[:, :3] @ rng.normal(size=(3, G)) \
        + rng.normal(size=(n, G))
    Y[rng.random(Y.shape) < 0.02] = np.nan
    Kb = {"E": Ze @ Ze.T / 200, "S": Zs @ Zs.T / 40, "N": Zn @ Zn.T / 400}
    R = Y - np.nanmean(Y, 0)
    theta, rep = tune_block_weights(Kb, R, Y, np.ones(G, int), lam_grid=np.logspace(-2, 3, 11),
                                    merge_frac=0.0, min_rows=20, max_eval=25)
    assert rep["tune_n_eval"] <= 25 and len(rep["tune_path"]) == rep["tune_n_eval"]
    assert rep["tune_path"][0]["theta"] == {"E": 1.0, "S": 0.0, "N": 0.0}  # anchor-only start
    assert theta["E"] == 1.0 and theta["S"] > 0 and theta["N"] == 0.0
    assert rep["tune_score"] == max(p["score"] for p in rep["tune_path"])
    _, rep2 = tune_block_weights(Kb, R, Y, np.ones(G, int), lam_grid=np.logspace(-2, 3, 11),
                                 merge_frac=0.0, min_rows=20, max_eval=4)
    assert rep2["tune_n_eval"] == 4


# ---------------------------------------------------------------------------- 2. real path
def _real_blocks(names, n_cols=600):
    import pyarrow.parquet as pq
    out = {}
    for b in names:
        path = FEAT / f"{b}.parquet"
        f = pq.ParquetFile(path)
        idx = json.loads(f.schema_arrow.metadata[b"pandas"])["index_columns"]
        cols = [c for c in f.schema_arrow.names if c not in idx][:n_cols]
        out[b] = pd.read_parquet(path, columns=cols)
    return out


@pytest.mark.skipif(not (FEAT / "E.parquet").exists(), reason="feature cache not built")
@pytest.mark.parametrize("params", [{"modalities": ["E"], "residual": "B0"},
                                    {"modalities": ["E", "C"], "residual": "B0", "tune_weights": True,
                                     "max_tune_eval": 6}])
def test_real_path_held_out_feature_rows_cannot_change_fit(monkeypatch, params):
    raw = _real_blocks(KR.raw_names(KR.block_specs(params["modalities"])))
    sp = pd.read_csv(ROOT / "configs" / "splits.csv", index_col=0)
    dev = sp.index[(sp["role"] == "dev") & sp.index.isin(raw["E"].index)]
    rng = np.random.default_rng(8)
    tr = pd.Index(sorted(rng.choice(dev, size=110, replace=False)))
    E = raw["E"].loc[tr].to_numpy()
    E = (E - E.mean(0)) / (E.std(0) + 1e-9)
    G = 25
    Y = pd.DataFrame(E[:, :10] @ rng.normal(size=(10, G)) * 0.3 + rng.normal(size=(110, G)),
                     index=tr, columns=[f"g{j}" for j in range(G)])
    Y.iloc[:20, :4] = np.nan  # a second NaN pattern
    Y.iloc[:, :12] *= 3       # some 'variable' genes for the tuning criterion
    lin = sp.loc[tr, "OncotreeLineage"]
    groups = sp.loc[tr, "PatientID"].to_numpy()
    kw = dict(params, lam_grid=[1e2, 1e3, 1e4], min_rows=40)

    def fit_predict(blocks):
        monkeypatch.setattr(KR, "_RAW_CACHE", dict(blocks))
        m = KR.KernelRidgeM1(**kw).fit(Y, lineage=lin, groups=groups)
        return m, m.predict(tr).to_numpy()

    m1, P1 = fit_predict(raw)
    pert = {}
    for b, df in raw.items():  # every non-training row (held-out dev, test, others) is perturbed
        d = df.copy()
        other = ~d.index.isin(tr)
        vals = rng.normal(size=(int(other.sum()), d.shape[1])) * 50 + 7
        vals[rng.random(vals.shape) < 0.1] = np.nan
        d.loc[other] = vals.astype(d.to_numpy().dtype)
        pert[b] = d
    m2, P2 = fit_predict(pert)
    assert np.array_equal(P1, P2)
    assert m1.lam_by_group_ == m2.lam_by_group_
    assert getattr(m1, "theta_", None) == getattr(m2, "theta_", None)
    if params.get("tune_weights"):
        assert m1.report_["tune_n_eval"] <= 6 and m1.theta_["E"] == 1.0
