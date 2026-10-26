"""Brute-force check of the pattern-wise kernel ridge + closed-form leave-group-out (reviewer code).

Covers what tests/test_kernel_ridge.py does not combine: p > n (the real regime: 18k features,
~800 rows), merged NaN patterns (intersection rows), a fallback pattern, patient groups that
straddle NaN-pattern boundaries, and a B1 leave-line-out residual offset (R != Y).

For each gene we compare the code's LOO criterion r(lambda) with
  (a) "fixed-target" brute force: drop patient group, primal ridge (unpenalised intercept) on the
      remaining cluster rows of the precomputed residual R, add back the offset y - R;
  (b) "full-pipeline" brute force (B1 only): drop the group, REBUILD B1 and its leave-line-out
      residuals from the remaining lines, fit, and predict the group with B1-without-group.
(a) must match to machine precision; (b) shows how far the closed form is from a true LOGO of
the whole B1-residual pipeline. Prediction at new points is checked against primal ridge.

    python reviews/g2_review/check_solver.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

for v in ["OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"]:
    os.environ[v] = "2"
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from essml.baselines import LineageMean, leave_line_out_lineage_residual  # noqa: E402
from essml.models.kernel_ridge import fit_pattern_ridge, predict_pattern_ridge  # noqa: E402


def ridge(Z, y, lam):
    zm, ym = Z.mean(0), y.mean()
    Zc = Z - zm
    # dual form (p > n): beta = Zc' (Zc Zc' + lam I)^-1 (y - ym)
    a = np.linalg.solve(Zc @ Zc.T + lam * np.eye(len(y)), y - ym)
    beta = Zc.T @ a
    return lambda Zn: ym + (Zn - zm) @ beta


def ridge_primal(Z, y, lam):
    zm, ym = Z.mean(0), y.mean()
    Zc = Z - zm
    beta = np.linalg.solve(Zc.T @ Zc + lam * np.eye(Z.shape[1]), Zc.T @ (y - ym))
    return lambda Zn: ym + (Zn - zm) @ beta


def corr(a, b):
    return float(np.corrcoef(a, b)[0, 1])


def make(seed=3, n=90, p=150, G=10):
    rng = np.random.default_rng(seed)
    Z = rng.normal(size=(n, p)) * rng.uniform(0.3, 2, size=p)
    lin = rng.choice(list("abcd"), size=n, p=[0.4, 0.3, 0.2, 0.1])
    leff = {k: rng.normal(size=G) for k in "abcd"}
    Y = Z[:, :20] @ rng.normal(size=(20, G)) * 0.15 + np.array([leff[l] for l in lin]) \
        + rng.normal(size=(n, G))
    Y[:, 3] += 40.0
    Y[:30] += 25.0  # NaN-region rows have very different values (centring bugs would show)
    obs = np.ones((n, G), bool)
    obs[:30, 4] = obs[:30, 5] = False            # pattern A
    obs[:29, 8] = False                          # A's rows plus row 29 -> merges into A (tol=1);
                                                 # gene 8 is fitted on A's intersection rows
    obs[50, 6] = False                           # full minus one row: the merge guard keeps it
                                                 # separate (would shrink the 4-gene full cluster)
    obs[60:63, 7] = False                        # differs by 3 -> own cluster
    obs[30:, 9] = False                          # 30 rows -> fallback (min_rows 40)
    # patient groups: mostly singletons, some pairs/triples straddling row 30 / 60 boundaries
    groups = np.arange(n).astype(object)
    for a, b in [(28, 31), (29, 30), (59, 61), (62, 70), (10, 11), (12, 13)]:
        groups[b] = groups[a]
    groups[14] = groups[12]
    return Z, np.where(obs, Y, np.nan), obs, lin, groups.astype(str)


def main():
    Z, Y, obs, lin, groups = make()
    n, G = Y.shape
    K = Z @ Z.T
    grid = np.array([1.0, 30.0, 300.0, 3000.0])
    ids = pd.Index([f"L{i}" for i in range(n)])
    Ydf = pd.DataFrame(Y, index=ids)
    lser = pd.Series(lin, index=ids)
    kappa = 5.0
    worst = {"B0_fixed": 0.0, "B1_fixed": 0.0, "pred": 0.0}
    b1_gap, argmax_flip = [], 0
    for resid in ["B0", "B1"]:
        if resid == "B1":
            R = leave_line_out_lineage_residual(Ydf, lser, kappa).to_numpy()
        else:
            R = Y - np.nanmean(Y, 0)
        sol = fit_pattern_ridge(K, R, Y, np.zeros(G, int), groups=groups, lam_grid=grid,
                                merge_frac=0.012, min_rows=40)
        rep = sol["report"]
        assert rep["n_patterns"] == 6 and rep["n_patterns_merged"] == 1 and rep["n_fallback_genes"] == 1, rep
        cl_of = {int(g): c for c in sol["clusters"] for g in c["genes"]}
        Znew = np.random.default_rng(9).normal(size=(5, Z.shape[1]))
        Pnew = predict_pattern_ridge(Znew @ Z.T, sol["clusters"], G)
        for g in range(G):
            if g not in cl_of:
                assert np.all(Pnew[:, g] == 0)  # fallback: offset only
                continue
            rows = cl_of[g]["rows"]
            # prediction at new points with the selected lambda
            lam_g = sol["lam"][g]
            f = ridge_primal(Z[rows], R[rows, g], lam_g)
            worst["pred"] = max(worst["pred"], float(np.max(np.abs(f(Znew) - Pnew[:, g]))))
            gr = groups[rows]
            bf_fixed, bf_full = [], []
            for lam in grid:
                p_fix = np.empty(len(rows))
                p_full = np.empty(len(rows))
                for gid in np.unique(gr):
                    te = gr == gid
                    trr = rows[~te]
                    p_fix[te] = (Y[rows[te], g] - R[rows[te], g]) + ridge(Z[trr], R[trr, g], lam)(Z[rows[te]])
                    if resid == "B1":
                        # full pipeline: rebuild B1 + LLO residuals without the group (all observed rows)
                        keep = ids[(groups != gid)]
                        Yk = Ydf.loc[keep, [g]]
                        Rk = leave_line_out_lineage_residual(Yk, lser.loc[keep], kappa)
                        b1 = LineageMean(kappa).fit(Yk, lser.loc[keep])
                        off = b1.predict(ids[rows[te]], lser.loc[ids[rows[te]]]).to_numpy()[:, 0]
                        Rk_rows = Rk.loc[ids[trr]].to_numpy()[:, 0]
                        p_full[te] = off + ridge(Z[trr], Rk_rows, lam)(Z[rows[te]])
                bf_fixed.append(corr(p_fix, Y[rows, g]))
                if resid == "B1":
                    bf_full.append(corr(p_full, Y[rows, g]))
            code = sol["r_grid"][:, g]
            worst[f"{resid}_fixed"] = max(worst[f"{resid}_fixed"], float(np.max(np.abs(code - bf_fixed))))
            if resid == "B1":
                b1_gap.append(np.max(np.abs(code - np.array(bf_full))))
                argmax_flip += int(np.argmax(code) != np.argmax(bf_full))
    print("max |code - brute force| LOGO r (B0 residual, fixed targets):", f"{worst['B0_fixed']:.2e}")
    print("max |code - brute force| LOGO r (B1 residual, fixed targets):", f"{worst['B1_fixed']:.2e}")
    print("max |code - primal ridge| prediction at new points:          ", f"{worst['pred']:.2e}")
    print("B1 residual: max |closed form - full-pipeline LOGO| r per gene:",
          f"{max(b1_gap):.4f} (median {np.median(b1_gap):.4f}); argmax lambda differs for "
          f"{argmax_flip}/{len(b1_gap)} genes")
    ok = worst["B0_fixed"] < 1e-8 and worst["B1_fixed"] < 1e-8 and worst["pred"] < 1e-8
    print("PASS" if ok else "FAIL")
    realistic_gap()


def realistic_gap(seed=5, n=300, p=400, G=8, kappa=10.0):
    """B1 closed form vs full-pipeline LOGO on a less pathological toy: lineage-informative
    features, 6 lineages of 20-90 lines, ~5% patient pairs, no row shifts, all observed."""
    rng = np.random.default_rng(seed)
    sizes = [90, 70, 50, 40, 30, 20]
    lin = np.repeat(list("abcdef"), sizes)
    cent = {k: rng.normal(size=p) * 0.7 for k in "abcdef"}
    Z = np.array([cent[l] for l in lin]) + rng.normal(size=(n, p))
    leff = {k: rng.normal(size=G) * 0.5 for k in "abcdef"}
    Y = (Z[:, :30] @ rng.normal(size=(30, G))) * 0.05 + np.array([leff[l] for l in lin]) \
        + rng.normal(size=(n, G))
    groups = np.arange(n).astype(object)
    for a in rng.choice(n - 1, size=8, replace=False):
        groups[a + 1] = groups[a]
    groups = groups.astype(str)
    ids = pd.Index([f"L{i}" for i in range(n)])
    Ydf, lser = pd.DataFrame(Y, index=ids), pd.Series(lin, index=ids)
    grid = np.array([10.0, 100.0, 1000.0, 1e4, 1e5])
    R = leave_line_out_lineage_residual(Ydf, lser, kappa).to_numpy()
    sol = fit_pattern_ridge(Z @ Z.T, R, Y, np.zeros(G, int), groups=groups, lam_grid=grid,
                            merge_frac=0.0, min_rows=10)
    P = np.zeros((len(grid), n, G))
    for gid in np.unique(groups):
        te = groups == gid
        keep = ids[~te]
        Rk = leave_line_out_lineage_residual(Ydf.loc[keep], lser.loc[keep], kappa).to_numpy()
        off = LineageMean(kappa).fit(Ydf.loc[keep], lser.loc[keep]).predict(
            ids[te], lser.loc[ids[te]]).to_numpy()
        Zt = Z[~te]
        zm = Zt.mean(0)
        Zc = Zt - zm
        Kc = Zc @ Zc.T
        rm = Rk.mean(0)
        for li, lam in enumerate(grid):
            a = np.linalg.solve(Kc + lam * np.eye(len(Zt)), Rk - rm)
            P[li, te] = off + rm + ((Z[te] - zm) @ Zc.T) @ a
    full = np.array([[corr(P[li, :, g], Y[:, g]) for g in range(G)] for li in range(len(grid))])
    code = sol["r_grid"]
    gap = np.abs(code - full)
    print("realistic toy (n=300, kappa=10): B1 closed-form vs full-pipeline LOGO r: max gap "
          f"{gap.max():.4f}, median {np.median(gap):.4f}; mean-over-genes criterion by lambda: "
          f"code {code.mean(1).round(4).tolist()} vs full {full.mean(1).round(4).tolist()}; "
          f"selected lambda code {grid[np.argmax(code.mean(1))]:g} vs full {grid[np.argmax(full.mean(1))]:g}")


if __name__ == "__main__":
    main()
