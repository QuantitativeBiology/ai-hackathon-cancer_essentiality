"""M4 pre-check and M4 factorised neural model (DESIGN.md §5, G3 nonlinear track).

All three models share the same inputs, fitted inside the outer training fold only:
    E (log2 TPM+1) -> FoldTransformer (training mean/SD, SD >= 0.1) -> clip +-8
    -> PCA on the TRAINING rows (``g3_common.TrainPCA``, 256 PCs; held-out rows only projected).
Target: y standardised per gene with the training mean and SD over observed rows (B0 residual
divided by the gene SD, so every gene has equal weight in the loss). Predictions are mapped back:
y_hat = mean_g + sd_g * out. NaN targets are never filled: the losses are masked.

- ``RidgePC`` ("M4pre_RidgePC"): multi-output ridge on the PCs, exact per NaN pattern (rows
  observed for the gene), unpenalised intercept. lambda (relative to n) and PC whitening are chosen
  on ONE inner grouped (PatientID) 80/20 split of the training lines; criterion = mean per-gene r
  on the inner validation rows over variable genes of the inner training part. Refit on all
  training lines with the chosen setting.
- ``NeuralPC`` (arch="mlp", "M4pre_MLP"): Linear(k, hidden) -> ReLU -> Dropout -> Linear(hidden, G),
  masked MSE, AdamW. Early stopping on the same inner split and criterion (evaluated every epoch,
  patience); the small grid (whitening x dropout) is chosen on the same split. Refit on all
  training lines for the chosen number of epochs.
- ``NeuralPC`` (arch="factorised", "M4"): y[c, g] = b_g + <u(x_c), v_g>, u = MLP encoder ->
  emb_dim; v_g (G x emb_dim) initialised from the top right singular vectors of the TRAINING
  target matrix (NaN -> 0 after standardisation), b_g = 0 (targets are centred).

``seed`` (configs/seeds.json model_seeds) sets torch initialisation and minibatch order only; the
inner split uses ``split_seed`` so all seeds see the same split.
"""
from __future__ import annotations

import os
import time

import numpy as np
import pandas as pd

from .g3_common import (CellFeatures, TrainPCA, default_groups, inner_split, masked_r,
                        masked_stats, variable_mask)

RIDGE_REL_LAM = np.logspace(-4, 2, 13)  # lambda / n_rows


# ============================================================================ shared PC inputs
class _PCInputs:
    uses_groups = True

    def _fit_inputs(self, Y, lineage, X, groups):
        ids = pd.Index([str(i) for i in Y.index])
        if groups is None:
            groups = default_groups(ids)
        self.train_ids_ = ids
        self.genes_ = Y.columns
        lin = None if lineage is None else pd.Series(lineage).reindex(Y.index).astype(str)
        tr, va = inner_split(ids, np.asarray(groups), lin, seed=self.split_seed)
        # inner phase: feature transform fitted on the inner-training lines only
        cf_in, Ztr_in = CellFeatures(("E",), permute_seed=self.permute_seed).fit(ids[tr], X)
        Zva_in = cf_in.transform(ids[va], X)
        # final fit: transform fitted on all (outer) training lines
        self.cf_, Z = CellFeatures(("E",), permute_seed=self.permute_seed).fit(ids, X)
        # target standardisation: inner phase uses inner-training statistics
        Yv = Y.to_numpy(np.float64)
        mu_in, sd_in = masked_stats(Yv[tr])
        Ys_in = ((Yv[tr] - mu_in) / sd_in, (Yv[va] - mu_in) / sd_in)
        self.mu_, self.sd_ = masked_stats(Yv)
        Ys = (Yv - self.mu_) / self.sd_
        return ids, (Ztr_in, Zva_in), Ys_in, Z, Ys, tr, va

    def _pcs(self, Z, whiten):
        return TrainPCA(self.n_pc, whiten=whiten).fit(Z)

    def _pred_inputs(self, index, X):
        Z = self.cf_.transform(pd.Index([str(i) for i in index]), X)
        return self.pca_.transform(Z)

    def _out(self, index, S):
        return pd.DataFrame(self.mu_ + self.sd_ * S, index=pd.Index(index), columns=self.genes_)


# ============================================================================ ridge on PCs
def ridge_patterns(P, Ys, lams, min_rows=10):
    """Exact per-NaN-pattern ridge with unpenalised intercept.

    P: (n, k) inputs; Ys: (n, G) targets with NaN; lams: absolute lambdas.
    Returns coef (L, k, G) and intercept (L, G); genes with < min_rows rows get coef 0 and
    intercept = observed mean (0 if none)."""
    P = np.asarray(P, np.float64)
    n, k = P.shape
    G = Ys.shape[1]
    obs = ~np.isnan(Ys)
    pats, inv = np.unique(obs.T, axis=0, return_inverse=True)
    inv = np.asarray(inv).ravel()
    coef = np.zeros((len(lams), k, G))
    icpt = np.zeros((len(lams), G))
    for j, rp in enumerate(pats):
        genes = np.flatnonzero(inv == j)
        rows = np.flatnonzero(rp)
        if len(rows) < min_rows:
            if len(rows):
                icpt[:, genes] = np.nanmean(Ys[np.ix_(rows, genes)], 0)
            continue
        Pr = P[rows]
        pm = Pr.mean(0)
        Pc = Pr - pm
        yr = Ys[np.ix_(rows, genes)]
        ym = yr.mean(0)
        d, V = np.linalg.eigh(Pc.T @ Pc)
        d = np.clip(d, 0, None)
        B = V.T @ (Pc.T @ (yr - ym))
        for li, lam in enumerate(lams):
            c = V @ (B / (d + lam)[:, None])
            coef[li][:, genes] = c
            icpt[li, genes] = ym - pm @ c
    return coef, icpt


class RidgePC(_PCInputs):
    name = "M4pre_RidgePC"

    def __init__(self, n_pc=256, rel_lam_grid=None, whiten_grid=(False, True), split_seed=20261026,
                 permute_seed=None, seed=0):
        self.n_pc = int(n_pc)
        self.rel_lam_grid = np.asarray(RIDGE_REL_LAM if rel_lam_grid is None else rel_lam_grid)
        self.whiten_grid = [bool(w) for w in whiten_grid]
        self.split_seed = int(split_seed)
        self.permute_seed = permute_seed
        self.seed = seed  # unused (deterministic); kept for config symmetry

    def fit(self, Y: pd.DataFrame, lineage=None, X=None, groups=None):
        t0 = time.perf_counter()
        ids, (Ztr_in, Zva_in), (Ytr_in, Yva_in), Z, Ys, tr, va = self._fit_inputs(Y, lineage, X, groups)
        var_in = variable_mask(Y.to_numpy(np.float64)[tr])
        table = []
        best = (-np.inf, None, None)
        for wh in self.whiten_grid:
            pca = self._pcs(Ztr_in, wh)
            Ptr, Pva = pca.transform(Ztr_in), pca.transform(Zva_in)
            lams = self.rel_lam_grid * len(tr)
            coef, icpt = ridge_patterns(Ptr, Ytr_in, lams)
            for li, rl in enumerate(self.rel_lam_grid):
                S = Pva @ coef[li] + icpt[li]
                sc = float(np.nanmean(masked_r(S[:, var_in], Yva_in[:, var_in])))
                table.append({"whiten": wh, "rel_lam": float(rl), "inner_P1": round(sc, 5)})
                if sc > best[0]:
                    best = (sc, wh, float(rl))
        _, wh, rl = best
        self.whiten_, self.rel_lam_ = wh, rl
        self.pca_ = self._pcs(Z, wh)
        coef, icpt = ridge_patterns(self.pca_.transform(Z), Ys, [rl * len(ids)])
        self.coef_, self.icpt_ = coef[0], icpt[0]
        self.report_ = {"chosen": {"whiten": wh, "rel_lam": rl}, "inner_P1": round(best[0], 5),
                        "grid": table, "n_pc": int(self.pca_.components_.shape[1]),
                        "n_inner_train": int(len(tr)), "n_inner_val": int(len(va)),
                        "fit_s": round(time.perf_counter() - t0, 1)}
        return self

    def predict(self, index, lineage=None, X=None) -> pd.DataFrame:
        P = self._pred_inputs(index, X).astype(np.float64)
        return self._out(index, P @ self.coef_ + self.icpt_)


# ============================================================================ neural (MLP / factorised)
def _torch():
    import torch
    nthr = int(os.environ.get("OMP_NUM_THREADS", "6") or 6)
    torch.set_num_threads(max(1, nthr))
    return torch


def _build_net(torch, k, G, arch, hidden, dropout, emb_dim, V_init=None):
    nn = torch.nn
    if arch == "mlp":
        return nn.Sequential(nn.Linear(k, hidden), nn.ReLU(), nn.Dropout(dropout), nn.Linear(hidden, G))

    class Factorised(nn.Module):
        def __init__(self):
            super().__init__()
            self.enc = nn.Sequential(nn.Linear(k, hidden), nn.ReLU(), nn.Dropout(dropout),
                                     nn.Linear(hidden, emb_dim))
            self.V = nn.Parameter(torch.randn(G, emb_dim) / np.sqrt(emb_dim))
            self.b = nn.Parameter(torch.zeros(G))
            if V_init is not None:
                with torch.no_grad():
                    self.V.copy_(torch.as_tensor(V_init, dtype=torch.float32))

        def forward(self, x):
            return self.enc(x) @ self.V.T + self.b
    return Factorised()


def svd_gene_init(Ys, emb_dim):
    """Gene embeddings from the training target matrix: top right singular vectors x S/sqrt(n)."""
    A = np.nan_to_num(np.asarray(Ys, np.float64))
    n = A.shape[0]
    K = A @ A.T
    d, U = np.linalg.eigh(K)
    order = np.argsort(d)[::-1][:emb_dim]
    d, U = np.clip(d[order], 1e-12, None), U[:, order]
    s = np.sqrt(d)
    W = (A.T @ U) / s           # (G, emb) right singular vectors
    return (W * (s / np.sqrt(n))).astype(np.float32)


def masked_mse(out, y, m):
    """Mean squared error over observed cells only (m = 1 observed, 0 masked; y finite)."""
    return (((out - y) ** 2) * m).sum() / m.sum().clamp(min=1.0)


def train_net(torch, P_tr, Y_tr, P_va, Y_va, var_va, *, arch, hidden, dropout, emb_dim, lr,
              weight_decay, batch, max_epochs, patience, seed, V_init=None, fixed_epochs=None):
    """Masked-MSE training. With validation data: early stopping on mean per-gene r (var_va genes);
    returns (net, best_epoch, best_score, history). With fixed_epochs: train exactly that many."""
    torch.manual_seed(int(seed))
    rng = np.random.default_rng(int(seed))
    k, G = P_tr.shape[1], Y_tr.shape[1]
    net = _build_net(torch, k, G, arch, hidden, dropout, emb_dim, V_init)
    opt = torch.optim.AdamW(net.parameters(), lr=lr, weight_decay=weight_decay)
    Xt = torch.as_tensor(np.asarray(P_tr, np.float32))
    M = ~np.isnan(Y_tr)
    Yt = torch.as_tensor(np.where(M, Y_tr, 0.0).astype(np.float32))
    Mt = torch.as_tensor(M.astype(np.float32))
    n = len(Xt)
    n_ep = int(fixed_epochs) if fixed_epochs is not None else int(max_epochs)
    best, best_ep, best_state, hist, wait = -np.inf, 0, None, [], 0
    if P_va is not None:
        Xv = torch.as_tensor(np.asarray(P_va, np.float32))
        Tv = np.asarray(Y_va, np.float64)[:, var_va]
    for ep in range(1, n_ep + 1):
        net.train()
        perm = rng.permutation(n)
        for s in range(0, n, batch):
            b = perm[s:s + batch]
            out = net(Xt[b])
            m = Mt[b]
            loss = masked_mse(out, Yt[b], m)
            opt.zero_grad()
            loss.backward()
            opt.step()
        if P_va is None:
            continue
        net.eval()
        with torch.no_grad():
            Sv = net(Xv).numpy()[:, var_va]
        sc = float(np.nanmean(masked_r(Sv, Tv)))
        hist.append(round(sc, 5))
        if sc > best:
            best, best_ep, wait = sc, ep, 0
            best_state = {k_: v.detach().clone() for k_, v in net.state_dict().items()}
        else:
            wait += 1
            if wait >= patience:
                break
    if best_state is not None:
        net.load_state_dict(best_state)
    net.eval()
    return net, best_ep, best, hist


class NeuralPC(_PCInputs):
    name = "M4pre_MLP"

    def __init__(self, arch="mlp", n_pc=256, hidden=512, emb_dim=64, grid=None, lr=1e-3,
                 weight_decay=1e-4, batch=64, max_epochs=200, patience=20, seed=0,
                 split_seed=20261026, gene_init="svd", permute_seed=None):
        if arch not in ("mlp", "factorised"):
            raise ValueError("arch must be 'mlp' or 'factorised'")
        self.arch = arch
        self.n_pc, self.hidden, self.emb_dim = int(n_pc), int(hidden), int(emb_dim)
        self.grid = grid or {"whiten": [False, True], "dropout": [0.2, 0.5]}
        self.lr, self.weight_decay = float(lr), float(weight_decay)
        self.batch, self.max_epochs, self.patience = int(batch), int(max_epochs), int(patience)
        self.seed, self.split_seed = int(seed), int(split_seed)
        self.gene_init = gene_init
        self.permute_seed = permute_seed

    def _points(self):
        pts = [{}]
        for key in sorted(self.grid):
            pts = [dict(p, **{key: v}) for p in pts for v in self.grid[key]]
        return pts

    def _v_init(self, Ys):
        if self.arch != "factorised" or self.gene_init != "svd":
            return None
        return svd_gene_init(Ys, self.emb_dim)

    def fit(self, Y: pd.DataFrame, lineage=None, X=None, groups=None):
        torch = _torch()
        t0 = time.perf_counter()
        ids, (Ztr_in, Zva_in), (Ytr_in, Yva_in), Z, Ys, tr, va = self._fit_inputs(Y, lineage, X, groups)
        var_in = variable_mask(Y.to_numpy(np.float64)[tr])
        table, best = [], (-np.inf, None, None)
        pcas = {}
        for p in self._points():
            wh = bool(p.get("whiten", False))
            if wh not in pcas:
                pca = self._pcs(Ztr_in, wh)
                pcas[wh] = (pca.transform(Ztr_in), pca.transform(Zva_in))
            Ptr, Pva = pcas[wh]
            t1 = time.perf_counter()
            _, ep, sc, hist = train_net(
                torch, Ptr, Ytr_in, Pva, Yva_in, var_in, arch=self.arch, hidden=self.hidden,
                dropout=float(p.get("dropout", 0.2)), emb_dim=self.emb_dim, lr=self.lr,
                weight_decay=float(p.get("weight_decay", self.weight_decay)), batch=self.batch,
                max_epochs=self.max_epochs, patience=self.patience, seed=self.seed,
                V_init=self._v_init(Ytr_in))
            table.append({**p, "best_epoch": ep, "inner_P1": round(sc, 5), "epochs_run": len(hist),
                          "s": round(time.perf_counter() - t1, 1), "history": hist})
            if sc > best[0]:
                best = (sc, p, ep)
        sc, p, ep = best
        self.chosen_ = dict(p, epochs=int(ep))
        wh = bool(p.get("whiten", False))
        self.pca_ = self._pcs(Z, wh)
        P_all = self.pca_.transform(Z)
        self.net_, _, _, _ = train_net(
            torch, P_all, Ys, None, None, None, arch=self.arch, hidden=self.hidden,
            dropout=float(p.get("dropout", 0.2)), emb_dim=self.emb_dim, lr=self.lr,
            weight_decay=float(p.get("weight_decay", self.weight_decay)), batch=self.batch,
            max_epochs=self.max_epochs, patience=self.patience, seed=self.seed,
            V_init=self._v_init(Ys), fixed_epochs=max(1, int(ep)))
        self.report_ = {"chosen": self.chosen_, "inner_P1": round(sc, 5), "grid": table,
                        "n_pc": int(self.pca_.components_.shape[1]), "arch": self.arch,
                        "n_inner_train": int(len(tr)), "n_inner_val": int(len(va)),
                        "torch_threads": int(torch.get_num_threads()),
                        "fit_s": round(time.perf_counter() - t0, 1)}
        return self

    def predict(self, index, lineage=None, X=None) -> pd.DataFrame:
        torch = _torch()
        P = self._pred_inputs(index, X)
        with torch.no_grad():
            S = self.net_(torch.as_tensor(np.asarray(P, np.float32))).numpy().astype(np.float64)
        return self._out(index, S)


class FactorisedM4(NeuralPC):
    name = "M4"

    def __init__(self, arch="factorised", **kw):
        super().__init__(arch=arch, **kw)
