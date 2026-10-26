"""M1 / M1b multi-output linear-kernel ridge, R1 and transfer model T (DESIGN.md §5, §10).

Pipeline inside one training fold (training IDs only reach any statistic):

1. Features. ``features.FoldTransformer`` (fitted on training IDs) standardises every raw block;
   standardised values are clipped to +-``clip`` (R1 algorithms review #11: rare binary features
   would otherwise reach z ~ 13). Columns are grouped into *modalities* (``MODALITY_BLOCKS``); each
   modality b is scaled by w_b = sqrt(S / p_b) (p_b = its number of fitted columns), so the linear
   kernel is the sum of per-modality kernels with default weights 1/sqrt(p_b) (M1b). The common
   constant S is p_E when E is used (single-block E is then exactly the plain linear kernel of
   M1) and 1e4 otherwise; it only shifts the lambda scale.
2. Target. M1/M1b: residual from B1 computed leave-line-out (``baselines.
   leave_line_out_lineage_residual``), kappa chosen by grouped inner CV of B1 on the training fold
   (same criterion as scripts/run.py: mean per-gene r on the inner-training variable genes).
   R1 / T: residual from B0 (= raw y with an unpenalised intercept; no lineage).
3. NaN-pattern-wise exact fit. Genes are grouped by their observed-row pattern on the training
   fold; patterns whose row sets differ by <= 1% of training lines are merged (the merged cluster
   uses the intersection of its members' rows, and no member loses more than 1% of its rows);
   clusters with < ``min_rows`` rows get the offset (B1/B0) only. Each cluster's kernel is the
   submatrix of ONE fold-level K, centred on the cluster's rows (= ridge with an unpenalised
   intercept fitted on those rows).
4. One float64 ``eigh`` per cluster; closed-form leave-one-out (leave-patient-out for training
   lines sharing a PatientID: block formula) over lambda in ``LAM_GRID``. lambda is shared within
   gene groups (variable = SD > 0.25 on the training fold / other; R1: one lambda for all genes)
   and chosen to maximise the mean per-gene Pearson r of the LOO predictions of y.
5. Prediction: offset(line) + rbar_c + (k(x, S_c) - kbar_c) @ alpha_c.

Saved fold models (``modelstore``): ``compact()`` drops the training design matrix (re-built from
cache/features on demand) and stores alpha as float32. ``predict(..., X=dict)`` takes raw blocks
(name -> DataFrame, rows = new lines) for external cells; without X the raw blocks come from
cache/features.
"""
from __future__ import annotations

import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import eigh

from essml import features as F
from essml import genesets, metrics
from essml.baselines import GeneMean, LineageMean, leave_line_out_lineage_residual

ROOT = Path(__file__).resolve().parents[3]
LAM_GRID = np.logspace(2, 7, 16)
KAPPA_GRID = [0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 10000]
STATUS_COLS = ["TP53_mut", "RB_altered"]  # model T: engineered-genotype flags only
MODALITY_BLOCKS = {
    "E": ["E"],
    "C": ["C", "C_instab"],
    "M": ["M", "M_burden"],
    "H": ["H"],
    "F": ["F_pair", "F_gene"],
    "Meta": ["meta"],
    "composites": ["composite", "p53score", "hallmark"],
    "P": ["P"],
    "Meta_prov": ["meta_prov"],
    "status": ["status"],
}
_RAW_CACHE: dict[str, pd.DataFrame] = {}


# ============================================================================ raw feature access
def load_raw_blocks(names, feature_dir=F.FEATURE_DIR) -> dict[str, pd.DataFrame]:
    """Raw (unfitted, row-wise) blocks from cache/features, memoised per process."""
    out = {}
    for n in names:
        if n not in _RAW_CACHE:
            if n == "status":
                _RAW_CACHE[n] = load_raw_blocks(["composite"], feature_dir)["composite"][STATUS_COLS]
            else:
                _RAW_CACHE[n] = F.load_raw([n], feature_dir)[n]
        out[n] = _RAW_CACHE[n]
    return out


def block_specs(modalities, gene_sets=None):
    """FoldTransformer specs for the modalities' raw blocks (DESIGN §2 defaults)."""
    names = [b for m in modalities for b in MODALITY_BLOCKS[m]]
    need_sets = any(b in ("p53score", "hallmark") for b in names)
    if need_sets and gene_sets is None:
        gene_sets = F.load_gene_sets()
    allspec = F.default_specs(gene_sets if need_sets else None,
                              include_provenance="meta_prov" in names)
    allspec["status"] = F.BlockSpec()
    return {b: allspec[b] for b in names}


def raw_names(specs):
    """Raw blocks needed by the specs (score blocks read their source block)."""
    out = []
    for b, s in specs.items():
        out.append(s.source if s.kind == "score" else b)
    return list(dict.fromkeys(out))


# ============================================================================ NaN-pattern clusters
def cluster_patterns(obs: np.ndarray, tol: int, min_rows: int):
    """Group genes (columns of the boolean observed mask ``obs``, n x G) by NaN pattern.

    Patterns are visited by decreasing gene count. A pattern joins the existing cluster with the
    smallest symmetric row difference if that difference is <= tol and, after intersecting the
    rows, no member (including the new one) has lost more than tol of its own observed rows.
    A join may remove rows from a cluster's existing genes only if the gene-rows lost there
    (removed rows x cluster genes) do not exceed tol x the joining pattern's genes, so a large
    cluster is never shrunk to absorb a handful of genes (they form their own cluster instead).
    Returns (clusters, fallback_genes, report); each cluster is {"rows": int idx, "genes": int idx}.
    Clusters with fewer than ``min_rows`` rows go to the fallback (offset only).
    """
    n, G = obs.shape
    pats, inv = np.unique(obs.T, axis=0, return_inverse=True)
    inv = np.asarray(inv).ravel()
    sizes = np.bincount(inv, minlength=len(pats))
    order = np.argsort(-sizes, kind="stable")
    cl = []  # {"rows": bool (n,), "members": [pattern idx]}
    n_merged = 0
    for p in order:
        rp = pats[p]
        best, best_d = None, None
        for c in cl:
            d = int((rp ^ c["rows"]).sum())
            if d > tol or (best_d is not None and d >= best_d):
                continue
            new = c["rows"] & rp
            removed = int((c["rows"] & ~rp).sum())  # rows the existing members would lose
            if removed * c["n_genes"] > tol * int(sizes[p]):
                continue  # never shrink a big cluster for a few genes (cost in gene-rows)
            if all(int((pats[q] & ~new).sum()) <= tol for q in c["members"] + [p]):
                best, best_d = c, d
        if best is None:
            cl.append({"rows": rp.copy(), "members": [p], "n_genes": int(sizes[p])})
        else:
            best["rows"] &= rp
            best["members"].append(p)
            best["n_genes"] += int(sizes[p])
            n_merged += 1
    clusters, fallback = [], []
    lost_rows = 0
    for c in cl:
        genes = np.flatnonzero(np.isin(inv, c["members"]))
        rows = np.flatnonzero(c["rows"])
        if len(rows) < min_rows:
            fallback.append(genes)
            continue
        lost_rows += int(sum((obs[:, g].sum() - len(rows)) for g in genes))
        clusters.append({"rows": rows, "genes": genes})
    fb = np.concatenate(fallback) if fallback else np.array([], dtype=int)
    rep = {"n_patterns": int(len(pats)), "n_clusters": len(clusters), "n_patterns_merged": n_merged,
           "n_fallback_genes": int(len(fb)), "cluster_rows": [int(len(c["rows"])) for c in clusters],
           "cluster_genes": [int(len(c["genes"])) for c in clusters],
           "gene_rows_dropped_by_merging": lost_rows}
    return clusters, fb, rep


# ============================================================================ core solver
def _colwise_r(P, T):
    """Pearson r per column (no NaN). Constant columns -> 0."""
    Pc = P - P.mean(0)
    Tc = T - T.mean(0)
    num = (Pc * Tc).sum(0)
    den = np.sqrt((Pc ** 2).sum(0) * (Tc ** 2).sum(0))
    with np.errstate(invalid="ignore", divide="ignore"):
        r = num / den
    return np.where(den > 0, r, 0.0)


def _group_blocks(groups):
    """Index arrays of groups with >= 2 members (positions into the given array)."""
    if groups is None:
        return []
    g = pd.Series(np.asarray(groups))
    return [np.asarray(ix) for ix in g.groupby(g).indices.values() if len(ix) > 1]


def loo_cluster(Kc, Rc, Yc, groups, lam_grid):
    """eigh of the centred cluster kernel + closed-form leave-one(-group)-out over lam_grid.

    Kc: centred kernel (n x n); Rc: centred residual targets (n x G); Yc: original targets on the
    same rows (criterion); groups: labels (n,) or None. Returns (d, V, U, r_grid (L x G), E_last).
    The LOO hat matrix of ridge with an unpenalised intercept is 11'/n + V diag(d/(d+lam)) V'.
    """
    n = Kc.shape[0]
    d, V = eigh(Kc)
    d = np.clip(d, 0.0, None)
    U = V.T @ Rc
    V2 = V ** 2
    blocks = _group_blocks(groups)
    r_grid = np.empty((len(lam_grid), Rc.shape[1]))
    # r(y - e, y) from per-gene sums: cov(p, y) = var(y) - cov(e, y),
    # var(p) = var(y) - 2 cov(e, y) + var(e)   (p = y - e; no (n x G) temporaries beyond e)
    yc = Yc - Yc.mean(0)
    syy = np.einsum("ij,ij->j", yc, yc)
    for li, lam in enumerate(lam_grid):
        s = d / (d + lam)
        E = V @ (s[:, None] * U)
        np.subtract(Rc, E, out=E)
        E *= (1.0 / (1.0 - (1.0 / n + V2 @ s)))[:, None]
        for ix in blocks:
            Vg = V[ix]
            Hgg = 1.0 / n + (Vg * s) @ Vg.T
            Fg = Vg @ (s[:, None] * U)
            E[ix] = np.linalg.solve(np.eye(len(ix)) - Hgg, Rc[ix] - Fg)
        se = E.sum(0)
        sey = np.einsum("ij,ij->j", E, yc)
        see = np.einsum("ij,ij->j", E, E) - se ** 2 / n
        cov = syy - sey
        var_p = syy - 2 * sey + see
        with np.errstate(invalid="ignore", divide="ignore"):
            r = cov / np.sqrt(var_p * syy)
        r_grid[li] = np.where((syy > 0) & (var_p > 1e-300), r, 0.0)
    return d, V, U, r_grid


def fit_pattern_ridge(K, R, Y, gene_group, groups=None, lam_grid=LAM_GRID, merge_frac=0.01,
                      min_rows=50, log=None):
    """Pattern-wise exact kernel ridge for all genes with grouped-LOO lambda.

    K: (n, n) linear kernel on training rows (uncentred). R: (n, G) targets to regress (NaN =
    unobserved). Y: (n, G) original targets for the LOO criterion (same NaN pattern; may equal R).
    gene_group: (G,) int labels; lambda is shared within each label. groups: (n,) patient labels
    for leave-group-out (None = singletons).
    Returns dict(clusters=[{rows, genes, alpha, rbar, kbar}], fallback, lam (per gene),
    lam_by_group, score_by_group (L x n_groups), r_grid (L x G), report).
    """
    t0 = time.perf_counter()
    R = np.asarray(R, np.float64)
    Y = np.asarray(Y, np.float64)
    n, G = R.shape
    obs = ~np.isnan(R)
    if not np.array_equal(obs, ~np.isnan(Y)):
        raise ValueError("R and Y must share the NaN pattern")
    tol = int(np.floor(merge_frac * n))
    clusters, fb, rep = cluster_patterns(obs, tol, min_rows)
    r_grid = np.full((len(lam_grid), G), np.nan)
    keep = []
    grp_arr = None if groups is None else np.asarray(groups)
    for c in clusters:
        rows, genes = c["rows"], c["genes"]
        Ks = K[np.ix_(rows, rows)]
        kbar = Ks.mean(0)
        Kc = Ks - kbar[None, :] - kbar[:, None] + kbar.mean()
        Rs = R[np.ix_(rows, genes)]
        rbar = Rs.mean(0)
        Ys = Y[np.ix_(rows, genes)]
        d, V, U, rg = loo_cluster(Kc, Rs - rbar, Ys, None if grp_arr is None else grp_arr[rows],
                                  lam_grid)
        r_grid[:, genes] = rg
        keep.append((c, d, V, U, rbar, kbar))
    gene_group = np.asarray(gene_group)
    labels = sorted(set(gene_group.tolist()))
    fitted = np.zeros(G, bool)
    for c in clusters:
        fitted[c["genes"]] = True
    score = np.full((len(lam_grid), len(labels)), np.nan)
    lam_by_group = {}
    for j, lab in enumerate(labels):
        m = fitted & (gene_group == lab)
        if m.any():
            score[:, j] = r_grid[:, m].mean(1)
            lam_by_group[lab] = float(lam_grid[int(np.argmax(score[:, j]))])  # first maximum
    lam_gene = np.array([lam_by_group.get(l, np.nan) for l in gene_group])
    out_clusters = []
    for c, d, V, U, rbar, kbar in keep:
        lg = lam_gene[c["genes"]]
        alpha = np.empty_like(U)
        for lam in np.unique(lg):
            m = lg == lam
            alpha[:, m] = V @ (U[:, m] / (d + lam)[:, None])
        out_clusters.append({"rows": c["rows"], "genes": c["genes"], "alpha": alpha,
                             "rbar": rbar, "kbar": kbar})
    rep = dict(rep, n_train=n, merge_tol_rows=tol, solve_s=round(time.perf_counter() - t0, 2),
               lam_by_group={str(k): v for k, v in lam_by_group.items()},
               lam_at_grid_edge={str(k): bool(v in (lam_grid[0], lam_grid[-1]))
                                 for k, v in lam_by_group.items()})
    if log:
        log(rep)
    return {"clusters": out_clusters, "fallback": fb, "lam": lam_gene, "lam_by_group": lam_by_group,
            "score_by_group": score, "labels": labels, "r_grid": r_grid, "report": rep}


def predict_pattern_ridge(k_new, clusters, G):
    """(m, G) residual predictions; genes outside every cluster get 0 (offset only)."""
    out = np.zeros((k_new.shape[0], G))
    for c in clusters:
        a = c["alpha"].astype(np.float64, copy=False)
        out[:, c["genes"]] = c["rbar"] + (k_new[:, c["rows"]] - c["kbar"]) @ a
    return out


# ============================================================================ B1 kappa by grouped CV
def tune_kappa(Y: pd.DataFrame, lineage: pd.Series, groups, seed, n_splits=5, grid=KAPPA_GRID,
               inner_min_obs=30):
    """Grouped inner CV of B1 on the training fold (criterion = scripts/run.py's)."""
    from sklearn.model_selection import StratifiedGroupKFold
    ids = Y.index
    grp = (pd.Series(np.arange(len(ids)), index=ids) if groups is None
           else pd.Series(np.asarray(groups), index=ids))
    sgk = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    scores = np.zeros((len(grid), n_splits))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        splits = list(sgk.split(np.zeros(len(ids)), lineage.loc[ids].to_numpy(), grp.to_numpy()))
    for j, (a, b) in enumerate(splits):
        Yi = Y.iloc[a]
        var = genesets.gene_classes(Yi, set())["variable"]
        g = var.index[var.to_numpy()]
        Tva = Y.iloc[b][g].to_numpy(np.float64)
        base = LineageMean(0.0).fit(Yi[g], lineage.loc[Yi.index])
        for i, k in enumerate(grid):
            base.kappa = float(k)
            base.means_ = base._shrunk(base.sum_, base.cnt_, base.mu_)
            P = base.predict(Y.index[b], lineage.loc[Y.index[b]]).to_numpy()
            scores[i, j] = np.nanmean(metrics.per_gene_pearson(P, Tva, inner_min_obs))
    mean = scores.mean(1)
    return float(grid[int(np.argmax(mean))]), mean


# ============================================================================ model classes
class KernelRidgeM1:
    """M1 (modalities=['E']) / M1b (several modalities, 1/sqrt(p_b) block weights)."""
    name = "M1"
    uses_groups = True
    residual = "B1"
    lam_groups = "variable"

    def __init__(self, modalities=("E",), clip=8.0, lam_grid=None, merge_frac=0.01, min_rows=50,
                 kappa=None, seed=20261026, scale=None, permute_seed=None, residual=None):
        self.modalities = list(modalities)
        self.clip = float(clip)
        self.lam_grid = np.asarray(LAM_GRID if lam_grid is None else lam_grid, np.float64)
        self.merge_frac, self.min_rows = float(merge_frac), int(min_rows)
        self.kappa = kappa
        self.seed = int(seed)
        self.scale = scale
        self.permute_seed = permute_seed  # negative control: feature rows shuffled across lines
        if residual is not None:  # diagnostic override of the class default ("B1" / "B0")
            if residual not in ("B0", "B1"):
                raise ValueError("residual must be 'B0' or 'B1'")
            self.residual = residual

    # -- features
    def _raw(self, X, ids=None):
        if isinstance(X, dict):
            return X
        if isinstance(X, pd.DataFrame):  # pre-built single block (tests / generic use)
            return {"X": X}
        raw = load_raw_blocks(raw_names(self.specs_))
        if self.permute_seed is not None:  # same permutation of ModelIDs for every block
            allids = sorted(set().union(*[set(b.index) for b in raw.values()]))
            perm = np.random.default_rng(int(self.permute_seed)).permutation(len(allids))
            remap = dict(zip(allids, np.asarray(allids)[perm]))
            raw = {k: v.rename(index=remap) for k, v in raw.items()}
        return raw

    def _setup_specs(self, X):
        if isinstance(X, pd.DataFrame):
            self.specs_ = {"X": F.BlockSpec()}
            self.mod_blocks_ = {"X": ["X"]}
        else:
            self.specs_ = block_specs(self.modalities)
            self.mod_blocks_ = {m: MODALITY_BLOCKS[m] for m in self.modalities}
        self.frame_mode_ = isinstance(X, pd.DataFrame)

    def _design(self, raw, ids):
        Xdf = self.ft_.transform(raw, ids)
        if list(Xdf.columns) != self.columns_:
            raise ValueError("feature columns changed between fit and transform")
        Z = Xdf.to_numpy(np.float64)
        np.clip(Z, -self.clip, self.clip, out=Z)
        Z *= self.col_w_[None, :]
        return Z

    def _fit_features(self, X, train_ids):
        self._setup_specs(X)
        raw = self._raw(X)
        self.ft_ = F.FoldTransformer(self.specs_, check_roles=not self.frame_mode_)
        self.ft_.fit(raw, train_ids)
        Xdf = self.ft_.transform(raw, train_ids)
        self.columns_ = list(Xdf.columns)
        blk = np.array([c.split(":", 1)[0] for c in self.columns_])
        p_mod = {m: int(np.isin(blk, bs).sum()) for m, bs in self.mod_blocks_.items()}
        S = self.scale or (p_mod.get("E") or p_mod.get("X") or 1e4)
        w = np.zeros(len(blk))
        for m, bs in self.mod_blocks_.items():
            if p_mod[m]:
                w[np.isin(blk, bs)] = np.sqrt(S / p_mod[m])
        self.col_w_, self.p_mod_, self.S_ = w, p_mod, float(S)
        Z = Xdf.to_numpy(np.float64)
        np.clip(Z, -self.clip, self.clip, out=Z)
        Z *= w[None, :]
        return Z

    # -- offset
    def _offset_fit(self, Y, lineage, groups):
        if self.residual == "B1":
            k = self.kappa
            if k is None:
                k, self.kappa_scores_ = tune_kappa(Y, lineage, groups, self.seed)
            self.kappa_ = float(k)
            self.offset_ = LineageMean(self.kappa_).fit(Y, lineage)
            return leave_line_out_lineage_residual(Y, lineage, self.kappa_).to_numpy(np.float64)
        self.kappa_ = None
        self.offset_ = GeneMean().fit(Y)
        return Y.to_numpy(np.float64) - self.offset_.mu_[None, :]

    def _gene_groups(self, Y):
        if self.lam_groups == "single":
            return np.zeros(Y.shape[1], int)
        var = genesets.gene_classes(Y, set())["variable"].to_numpy()
        return var.astype(int)  # 1 = variable, 0 = other (training-fold classes)

    def fit(self, Y: pd.DataFrame, lineage=None, X=None, groups=None):
        t0 = time.perf_counter()
        ids = pd.Index([str(i) for i in Y.index])
        if groups is None and not isinstance(X, pd.DataFrame):
            sp = pd.read_csv(ROOT / "configs" / "splits.csv", index_col=0)
            groups = sp["PatientID"].reindex(ids).fillna(pd.Series(ids, index=ids)).to_numpy()
        self.train_ids_ = ids
        self.genes_ = Y.columns
        Z = self._fit_features(X, ids)
        t_feat = time.perf_counter() - t0
        K = Z @ Z.T
        lin = None if lineage is None else pd.Series(lineage).reindex(Y.index)
        R = self._offset_fit(Y, lin, groups)
        Yv = Y.to_numpy(np.float64)
        gg = self._gene_groups(Y)
        sol = fit_pattern_ridge(K, R, Yv, gg, groups=groups, lam_grid=self.lam_grid,
                                merge_frac=self.merge_frac, min_rows=self.min_rows)
        self.clusters_ = sol["clusters"]
        self.lam_ = sol["lam"]
        self.lam_by_group_ = sol["lam_by_group"]
        self.score_by_group_ = sol["score_by_group"]
        self.report_ = dict(sol["report"], kappa=self.kappa_, n_features=len(self.columns_),
                            p_modality=self.p_mod_, feature_s=round(t_feat, 2),
                            fit_s=round(time.perf_counter() - t0, 2))
        self.Ztr_ = Z
        return self

    def _train_design(self):
        if getattr(self, "Ztr_", None) is None:
            self.Ztr_ = self._design(self._raw(None), self.train_ids_)
        return self.Ztr_

    def predict_residual(self, index, X=None):
        index = pd.Index([str(i) for i in index])
        raw = self._raw(X)
        Zn = self._design(raw, index)
        kn = Zn @ self._train_design().T
        return predict_pattern_ridge(kn, self.clusters_, len(self.genes_))

    def predict(self, index, lineage=None, X=None) -> pd.DataFrame:
        index = pd.Index(index)
        if self.residual == "B1":
            off = self.offset_.predict(index, lineage=pd.Series(lineage).reindex(index)).to_numpy()
        else:
            off = self.offset_.predict(index).to_numpy()
        P = off + self.predict_residual(index, X)
        return pd.DataFrame(P, index=index, columns=self.genes_)

    def compact(self):
        """Before pickling: drop the training design (rebuilt from cache) and store alpha float32."""
        if not getattr(self, "frame_mode_", False):
            self.Ztr_ = None
        for c in self.clusters_:
            c["alpha"] = c["alpha"].astype(np.float32)
        return self


class KernelRidgeR1(KernelRidgeM1):
    """R1: Chang & Zhang-style multi-output ridge on raw y (no residualisation), one lambda."""
    name = "R1"
    residual = "B0"
    lam_groups = "single"


class KernelRidgeT(KernelRidgeM1):
    """Transfer model T: E + TP53/RB1 status flags, no lineage (B0 residual), grouped lambda."""
    name = "T"
    residual = "B0"
    lam_groups = "variable"

    def __init__(self, modalities=("E", "status"), **kw):
        bad = set(modalities) - {"E", "status"}
        if bad:
            raise ValueError(f"model T uses only externally available features, not {bad}")
        super().__init__(modalities=modalities, **kw)
