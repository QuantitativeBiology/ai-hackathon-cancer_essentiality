"""G6 independent reproduction (reviewer's own code; essml.transfer / evaluate_transfer NOT imported).

Reads observed data directly (logged in results/test_access_log.csv as independent_reproduction_read):
  cache/harmonised/<arm>.parquet ('final', 'affine' for Bowden) and data/CRISPRGeneEffect.csv (set A + dev).
Per-fold predictions: fold_preds.parquet from predict_folds.py (persisted fold models; fold mean == saved table).
Writes only aggregate numbers to repro_out.json in this scratch directory.
"""
import json
import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score

ROOT = Path(r"C:\Users\ritab\Desktop\project")
HERE = Path(__file__).resolve().parent
DO_BOOT = "--boot" in sys.argv
OUT = {}

crit = json.loads((ROOT / "configs/transfer_criteria.json").read_text(encoding="utf-8"))
eval_genes = set(pd.read_csv(ROOT / "configs/eval_genes.csv", dtype={"entrez": str})["entrez"])
cs = json.loads((ROOT / "configs/control_split.json").read_text())
EXCL = set(cs["ceg_scaling"]) | set(cs["neg_scaling"]) | set(cs["scaling_set"])
NEG_EVAL = set(cs["neg_evaluation"])
CEG_ALL = set(cs["ceg_scaling"]) | set(cs["ceg_evaluation"])
gc = pd.read_csv(ROOT / "configs/gene_classes.csv", dtype={"entrez": str}).set_index("entrez")
FRAC = gc["frac_below"].astype(float)
COMP = pd.read_parquet(ROOT / "cache/transfer/comparators.parquet")
COMP.index = COMP.index.astype(str)
B0 = COMP["B0"]

# ------------------------------------------------------------------ metrics (own implementations)
def ap_w(score, y, w=None):
    """Average precision (sklearn step-wise semantics, tie groups) with integer multiplicity weights."""
    score = np.asarray(score, float)
    y = np.asarray(y, bool)
    w = np.ones(len(score)) if w is None else np.asarray(w, float)
    o = np.argsort(-score, kind="stable")
    s, yy, ww = score[o], y[o], w[o]
    last = np.r_[s[1:] != s[:-1], True]               # end of each tie group
    ctp = np.cumsum(ww * yy)[last]
    ctot = np.cumsum(ww)[last]
    if ctp[-1] <= 0:
        return np.nan
    dR = np.diff(np.r_[0.0, ctp]) / ctp[-1]
    with np.errstate(invalid="ignore", divide="ignore"):
        P = np.where(ctot > 0, ctp / ctot, 0.0)
    return float(np.sum(dR * P))


def recall_at_fpr_w(score, pos, neg, w=None, fpr=0.05):
    score = np.asarray(score, float)
    w = np.ones(len(score)) if w is None else np.asarray(w, float)
    ns, nw = score[neg], w[neg]
    o = np.argsort(-ns, kind="stable")
    ns, nw = ns[o], nw[o]
    n_neg, n_pos = nw.sum(), w[pos].sum()
    if n_neg <= 0 or n_pos <= 0:
        return np.nan
    k = math.floor(fpr * n_neg)
    j = int(np.argmax(np.cumsum(nw) > k))            # (k+1)-th highest negative, counting multiplicity
    thr = ns[j]
    return float(w[pos][score[pos] > thr].sum() / n_pos)


# sanity: own AP == sklearn on random data with ties and duplicated genes
rng = np.random.default_rng(0)
sc = np.round(rng.normal(size=3000), 1)
yy = rng.random(3000) < 0.1
ww = rng.integers(0, 3, 3000)
dup = np.repeat(np.arange(3000), ww)
OUT["selftest_ap_vs_sklearn"] = [abs(ap_w(sc, yy) - average_precision_score(yy, sc)),
                                 abs(ap_w(sc, yy, ww) - average_precision_score(yy[dup], sc[dup]))]


# ------------------------------------------------------------------ data
def arm(name, col="final"):
    d = pd.read_parquet(ROOT / "cache/harmonised" / f"{name}.parquet")
    d.index = d.index.astype(str)
    return d[col].astype(float)


def consensus(series_list):
    return pd.concat(series_list, axis=1).mean(axis=1, skipna=True).dropna()


def universe(obs):
    o = obs[np.isfinite(obs)]
    return pd.Index(sorted((set(o.index) & eval_genes) - EXCL))


FP = pd.read_parquet(HERE / "fold_preds.parquet")
def folds(tag, pid, tp):
    g = FP[(FP.tag == tag) & (FP.profile_id == pid) & (FP.tp53 == tp)]
    return {k: s.set_index("entrez")["pred"] for k, s in g.groupby("fold")}
T, NOEC = "t_E_status_24b46e2a", "t_E_status_noEC_fa0181c1"


def mean_dicts(ds):
    return {k: sum(d[k] for d in ds) / len(ds) for k in ds[0]}


# entries straight from the criteria file
def entries(bg, roles=("verdict",)):
    out = []
    for st, s in sorted(crit["backgrounds"][bg]["studies"].items()):
        if s["role"] not in roles:
            continue
        for tp, e in sorted(s["entries"].items()):
            out.append(dict(name=f"{st}:{tp}", study=st, tp53=tp, arms=e["arms"],
                            obs=consensus([arm(a) for a in e["arms"]])))
    return out


def entry_pred(bg, e, override=None):
    ei = crit["backgrounds"][bg]["expression_input_T"]
    profs = [override or ei.get("per_arm", {}).get(a) or ei["profile_id"] for a in e["arms"]]
    return mean_dicts([folds(T, p, e["tp53"]) for p in profs])


# ------------------------------------------------------------------ statistic per entry
def make_eval(e, q, U=None, frac_max=0.5, label="rank10"):
    U = universe(e["obs"]) if U is None else U
    x = e["obs"].reindex(U).to_numpy(float)
    ev = dict(U=U)
    if q == "Q3a":
        ev["y"] = (x <= np.quantile(x, 0.10)) if label == "rank10" else (x < -0.5)
        ev["valid"] = len(U) >= 1000 and ev["y"].sum() > 0
    else:
        neg = U.isin(list(NEG_EVAL))
        fr = FRAC.reindex(U).to_numpy(float)
        pos = (x < -0.5) & np.isfinite(fr) & (fr <= frac_max) & ~neg
        ev.update(pos=pos, neg=neg, n_pos=int(pos.sum()), n_neg=int(neg.sum()))
        ev["valid"] = len(U) >= 1000 and pos.sum() >= 10 and neg.sum() >= 100
    return ev


def metric(ev, q, pred, w=None):
    s = -pred.reindex(ev["U"]).to_numpy(float)
    assert np.isfinite(s).all()
    return ap_w(s, ev["y"], w) if q == "Q3a" else recall_at_fpr_w(s, ev["pos"], ev["neg"], w)


def entry_skill(ev, q, pf, comp, w=None):
    """D7: per-fold skill averaged (= mean fold metric - comparator metric)."""
    return float(np.mean([metric(ev, q, p, w) for p in pf.values()]) - metric(ev, q, comp, w))


def background(vals, studies):
    by = {}
    for (st, v) in vals:
        if np.isfinite(v):
            by.setdefault(st, []).append(v)
    sv = {s: float(np.mean(v)) for s, v in by.items()}
    return float(np.mean([sv[s] for s in studies if s in sv])), sv


# ------------------------------------------------------------------ paralog blocks and bootstrap (declared scheme)
PAIRS = pd.read_parquet(ROOT / "cache/priors/paralog_pairs.parquet")
def blocks_for(genes):
    genes = list(genes)
    idx = {g: i for i, g in enumerate(genes)}
    par = list(range(len(genes)))
    def f(i):
        while par[i] != i:
            par[i] = par[par[i]]
            i = par[i]
        return i
    p = PAIRS[np.minimum(PAIRS["identity"], PAIRS["rev_identity"]) >= 30.0]
    for a, b in zip(p["a"].astype(str), p["b"].astype(str)):
        if a in idx and b in idx:
            ra, rb = f(idx[a]), f(idx[b])
            if ra != rb:
                par[max(ra, rb)] = min(ra, rb)
    roots = [f(i) for i in range(len(genes))]
    code, out = {}, np.empty(len(genes), int)
    for i, r in enumerate(roots):          # label components by first appearance in sorted gene order
        out[i] = code.setdefault(r, len(code))
    return out, len(code)


class _Prep:
    """Fixed score order for one (entry, predictor); only weights change across replicates."""
    def __init__(self, ev, q, pred):
        s = -pred.reindex(ev["U"]).to_numpy(float)
        self.q = q
        if q == "Q3a":
            o = np.argsort(-s, kind="stable")
            ss = s[o]
            self.o, self.y = o, ev["y"][o].astype(float)
            self.last = np.r_[ss[1:] != ss[:-1], True]
        else:
            ni = np.flatnonzero(ev["neg"]); on = np.argsort(-s[ni], kind="stable")
            self.ni, self.ns = ni[on], s[ni][on]
            self.pi = np.flatnonzero(ev["pos"]); self.ps = s[self.pi]
    def __call__(self, w):
        if self.q == "Q3a":
            ww = w[self.o]
            ctp = np.cumsum(ww * self.y)[self.last]; ctot = np.cumsum(ww)[self.last]
            if ctp[-1] <= 0: return np.nan
            dR = np.diff(np.r_[0.0, ctp]) / ctp[-1]
            with np.errstate(invalid="ignore", divide="ignore"):
                P = np.where(ctot > 0, ctp / ctot, 0.0)
            return float(np.sum(dR * P))
        wn, wp = w[self.ni], w[self.pi]
        n_neg, n_pos = wn.sum(), wp.sum()
        if n_neg <= 0 or n_pos <= 0: return np.nan
        k = math.floor(0.05 * n_neg)
        thr = self.ns[int(np.argmax(np.cumsum(wn) > k))]
        return float(wp[self.ps > thr].sum() / n_pos)


def boot_background(evs, q, preds, comps, plan_genes, studies, n_boot=2000, seed=20261026):
    """Declared hierarchical bootstrap (studies, then paralog blocks), own implementation."""
    block_of, nb = blocks_for(plan_genes)
    pos_in_plan = {n: pd.Index(plan_genes).get_indexer(ev["U"]) for n, (st, ev) in evs.items()}
    prep = {n: ([_Prep(ev, q, p) for p in preds[n].values()], _Prep(ev, q, comps[n]))
            for n, (st, ev) in evs.items() if ev["valid"]}
    # check: unweighted prep == point metric
    for n, (mp, cp) in prep.items():
        st, ev = evs[n]
        assert abs((np.mean([m(np.ones(len(ev["U"]))) for m in mp]) - cp(np.ones(len(ev["U"])))) - entry_skill(ev, q, preds[n], comps[n])) < 1e-12
    ss = np.random.SeedSequence(seed).spawn(3)
    rs, rg = np.random.default_rng(ss[0]), np.random.default_rng(ss[1])
    studies = sorted(studies)
    reps = np.empty(n_boot)
    for r in range(n_boot):
        si = rs.integers(0, len(studies), len(studies))
        bc = np.bincount(rg.integers(0, nb, nb), minlength=nb)
        w = bc[block_of].astype(float)
        by = {}
        for n, (mp, cp) in prep.items():
            ww = w[pos_in_plan[n]]
            v = np.mean([m(ww) for m in mp]) - cp(ww)
            if np.isfinite(v):
                by.setdefault(evs[n][0], []).append(v)
        sv = {s_: np.mean(v) for s_, v in by.items()}
        d = [sv[studies[i]] for i in si if studies[i] in sv]
        reps[r] = np.mean(d) if d else np.nan
    ok = reps[np.isfinite(reps)]
    lo, hi = np.percentile(ok, [2.5, 97.5])
    return dict(ci_lo=float(lo), ci_hi=float(hi), se=float(np.std(ok, ddof=1)), n_valid=int(len(ok)),
                boot_mean=float(ok.mean()), n_blocks=nb)


# ------------------------------------------------------------------ set B
def ssmd(obs):
    c = obs.reindex(list(cs["ceg_evaluation"])).dropna()
    n = obs.reindex(list(cs["neg_evaluation"])).dropna()
    return float((c.mean() - n.mean()) / np.sqrt(c.var(ddof=1) + n.var(ddof=1))), float(n.median()), float(c.median())


RES = {}
for bg in ["RPE1", "hPSC"]:
    spec = crit["backgrounds"][bg]
    ver = entries(bg)
    allE = entries(bg, roles=("verdict", "predictor_only"))
    studies = sorted({e["study"] for e in ver})
    plan = sorted(set().union(*[set(universe(e["obs"])) for e in ver]))
    R = {"plan_n_genes": len(plan), "plan_n_blocks": blocks_for(plan)[1], "studies": studies}
    # QC tie-in to the frozen file
    R["consensus_qc"] = {e["name"]: ssmd(e["obs"]) for e in ver}
    R["final_eq_affine"] = {a: bool(np.allclose(arm(a, "final").to_numpy(), arm(a, "affine").to_numpy(), equal_nan=True))
                            for e in ver for a in e["arms"]}
    preds = {e["name"]: entry_pred(bg, e) for e in ver}
    b1col = spec["B1_nearest_column"]
    for key, (q, comp, fm, lab) in {"Q3a_primary": ("Q3a", B0, 0.5, "rank10"),
                                   "Q3a_effect_label": ("Q3a", B0, 0.5, "effect"),
                                   "Q3b_primary": ("Q3b", COMP[b1col], 0.5, "rank10"),
                                   "Q3b_vs_B0": ("Q3b", B0, 0.5, "rank10"),
                                   "Q3b_rare_vs_B1n": ("Q3b", COMP[b1col], 0.05, "rank10")}.items():
        evs = {e["name"]: (e["study"], make_eval(e, q, frac_max=fm, label=lab)) for e in ver}
        per = {n: entry_skill(ev, q, preds[n], comp) for n, (st, ev) in evs.items() if ev["valid"]}
        pt, sv = background([(evs[n][0], v) for n, v in per.items()], studies)
        ensv = {n: metric(ev, q, sum(preds[n].values()) / 5) - metric(ev, q, comp) for n, (st, ev) in evs.items()}
        ens_pt, _ = background([(evs[n][0], v) for n, v in ensv.items()], studies)
        absT = {n: float(np.mean([metric(ev, q, p) for p in preds[n].values()])) for n, (st, ev) in evs.items()}
        absC = {n: metric(ev, q, comp) for n, (st, ev) in evs.items()}
        # replicate ceiling point (LOSO consensus of the other studies' study profiles; B0 fill)
        prof = {}
        for e in allE:
            prof.setdefault(e["study"], []).append(e["obs"])
        prof = {s: consensus(v) for s, v in prof.items()}
        ceil = {}
        for n, (st, ev) in evs.items():
            others = [s for s in sorted(prof) if s != st]
            M = pd.concat([prof[s].reindex(ev["U"]) for s in others], axis=1)
            p = M.mean(axis=1, skipna=True).fillna(B0.reindex(ev["U"]))
            ceil[n] = metric(ev, q, p) - metric(ev, q, comp)
        cpt, csv_ = background([(evs[n][0], v) for n, v in ceil.items()], studies)
        R[key] = dict(point=pt, per_study=sv, per_entry=per, ensemble_point=ens_pt, absT=absT, absComp=absC,
                      ceiling_point=cpt, ceiling_per_entry=ceil,
                      n_pos={n: int(ev["y"].sum()) if q == "Q3a" else ev["n_pos"] for n, (st, ev) in evs.items()},
                      n_U={n: len(ev["U"]) for n, (st, ev) in evs.items()})
        if DO_BOOT and key in ("Q3a_primary", "Q3b_primary"):
            R[key]["boot"] = boot_background(evs, q, preds, {n: comp for n in evs}, plan, studies)
        # expression-swap null (primary statistics)
        if key in ("Q3a_primary", "Q3b_primary"):
            sw = {}
            for lab_, profs in spec["expression_swap_profiles"].items():
                for pid in profs:
                    sp = {e["name"]: entry_pred(bg, e, override=pid) for e in ver}
                    v = {n: entry_skill(ev, q, sp[n], comp) for n, (st, ev) in evs.items()}
                    swp, _ = background([(evs[n][0], x) for n, x in v.items()], studies)
                    sw[pid] = dict(point=swp)
                    if DO_BOOT:
                        sw[pid].update(boot_background(evs, q, sp, {n: comp for n in evs}, plan, studies))
            R[key]["swap"] = sw
    # diagnostics on the trivial component (Q3a primary)
    diag = {}
    exprs = pd.read_parquet(ROOT / "cache/external_expr/harmonised.parquet")
    exprs.columns = exprs.columns.astype(str)
    ce = set(gc.index[gc["common_essential"].astype(bool)])
    for e in ver:
        ev = make_eval(e, "Q3a")
        pids = [spec["expression_input_T"].get("per_arm", {}).get(a) or spec["expression_input_T"]["profile_id"] for a in e["arms"]]
        ex = exprs.loc[pids].mean(axis=0).reindex(ev["U"])
        pf = preds[e["name"]]
        # (i) trivial expression rule: B0, but 0 for genes not expressed (log2(TPM+1) < 1)
        triv = B0.reindex(ev["U"]).where(ex.fillna(0) >= 1.0, 0.0)
        d = {"triv_skill": metric(ev, "Q3a", triv) - metric(ev, "Q3a", B0),
             "frac_U_unexpressed": float((ex.fillna(0) < 1).mean())}
        # (ii) T skill restricted to expressed genes
        for nm, keep in {"expressed_only": ex.fillna(0).to_numpy() >= 1.0,
                         "no_depmap_common_ess": ~ev["U"].isin(list(ce)),
                         "no_CEGv2": ~ev["U"].isin(list(CEG_ALL))}.items():
            U2 = ev["U"][keep]
            ev2 = make_eval(e, "Q3a", U=U2)
            d[f"T_skill_{nm}"] = entry_skill(ev2, "Q3a", pf, B0)
            d[f"n_{nm}"] = int(len(U2))
        diag[e["name"]] = d
    R["diag_Q3a"] = diag
    RES[bg] = R
    print(bg, "done", flush=True)

OUT["setB"] = RES

# ------------------------------------------------------------------ hypotheses
rpe = {e["name"]: e for e in entries("RPE1")}
md = (sum(folds(T, "state:DepMap_RPE1", "WT").values()) / 5) - (sum(folds(T, "state:DepMap_RPE1", "null").values()) / 5)
od = rpe["brown2019:WT"]["obs"] - rpe["brown2019:null"]["obs"]
h1 = {g: dict(model=float(md.get(e, np.nan)), obs=float(od.get(e, np.nan))) for g, e in crit["hypotheses"]["H1"]["genes"].items()}
OUT["H1"] = dict(per_gene=h1, result="pass" if all(v["model"] < 0 and v["obs"] < 0 for v in h1.values()) else "fail")
g2 = [e for e in crit["hypotheses"]["H2"]["genes"].values() if e]
bd = (arm("bowden2020_WT_D19", "affine") - arm("bowden2020_KO_D19", "affine")).reindex(g2).dropna()
bwn = od.reindex(g2).dropna()
def h2(deltas):
    o = pd.concat(deltas, axis=1).mean(axis=1, skipna=True).dropna()
    g = [x for x in o.index if np.isfinite(md.get(x, np.nan))]
    m, oo = md.reindex(g).to_numpy(), o.reindex(g).to_numpy()
    rho = spearmanr(m, oo).statistic
    sg = float(np.mean(np.sign(m) == np.sign(oo)))
    return dict(n=len(g), spearman=float(rho), sign=sg, result="pass" if (len(g) >= 10 and rho > 0.3 and sg >= 0.7) else "fail")
OUT["H2"] = dict(pooled=h2([bwn, bd]), brown=h2([bwn]), bowden=h2([bd]))
mair = entries("hPSC")[0]
pT = mean_dicts([folds(T, "state:H1_iCas9_MEF", "WT"), folds(T, "state:H1_iCas9_laminin", "WT")])
pT = sum(pT.values()) / 5
pN = mean_dicts([folds(NOEC, "state:H1_iCas9_MEF", "WT"), folds(NOEC, "state:H1_iCas9_laminin", "WT")])
pN = sum(pN.values()) / 5
g3 = [e for e in crit["hypotheses"]["H3"]["genes"].values() if e]
dep = [g for g in g3 if np.isfinite(mair["obs"].get(g, np.nan)) and mair["obs"][g] < -0.5 and np.isfinite(pT.get(g, np.nan))]
OUT["H3"] = dict(n_obs_dep=len(dep), dep=dep, n_pred_T=int(sum(pT[g] < -0.5 for g in dep)),
                 n_pred_noEC=int(sum(pN[g] < -0.5 for g in dep)),
                 POU5F1_in_mair_arms={a: bool("5460" in arm(a, "raw").dropna().index) for a in mair["arms"]},
                 T_pred={g: round(float(pT[g]), 4) for g in g3}, noEC_pred={g: round(float(pN[g]), 4) for g in g3})
OUT["H3"]["result"] = ("inconclusive" if len(dep) < 3 else
                       ("pass" if OUT["H3"]["n_pred_T"] >= 0.5 * len(dep) else "fail"))

# ------------------------------------------------------------------ DepMap: set A + dev (comparator recomputation)
hdr = pd.read_csv(ROOT / "data/CRISPRGeneEffect.csv", nrows=0)
ent = [re.search(r"\((\d+)\)\s*$", c).group(1) if re.search(r"\((\d+)\)\s*$", c) else None for c in hdr.columns[1:]]
sp = pd.read_csv(ROOT / "configs/splits.csv", index_col=0)
dev = set(sp.index[sp.role == "dev"])
SETA = crit["set_A"]["tier_A1"]["lines"]
keep = dev | set(SETA)
chunks = []
for ch in pd.read_csv(ROOT / "data/CRISPRGeneEffect.csv", index_col=0, chunksize=200):
    chunks.append(ch[ch.index.isin(keep)])
Y = pd.concat(chunks)
Y.columns = ent
Y = Y.loc[:, [c is not None for c in Y.columns]]
Yd = Y.loc[sorted(dev & set(Y.index))]
OUT["n_dev_rows"] = len(Yd)
b0_re = Yd.mean(axis=0, skipna=True)
lin = sp["OncotreeLineage"].reindex(Yd.index)
cmp_ = {"B0": float(np.nanmax(np.abs(b0_re.reindex(COMP.index) - COMP["B0"])))}
for col in ["B1_CNS/Brain", "B1_Eye", "B1_Testis"]:
    l = col[3:]
    sub = Yd[lin == l]
    S, N = sub.sum(axis=0, skipna=True), sub.notna().sum(axis=0)
    b1 = (S + 10 * b0_re) / (N + 10)
    cmp_[col] = float(np.nanmax(np.abs(b1.reindex(COMP.index) - COMP[col])))
    cmp_[col + "_n_lines"] = int(len(sub))
fb = (Yd < -0.5).sum(axis=0) / Yd.notna().sum(axis=0)
cmp_["frac_below_maxdiff"] = float(np.nanmax(np.abs(fb.reindex(FRAC.index) - FRAC)))
cmp_["frac_below_crossing_0.5"] = int(((fb.reindex(FRAC.index) <= 0.5) != (FRAC <= 0.5)).sum())
OUT["comparator_recompute_maxabs"] = cmp_

YA = Y.loc[SETA]
in3 = YA.reindex(columns=sorted(eval_genes)).notna().sum(axis=0) >= 3
ref = pd.read_csv(ROOT / "results/test/t_E_status/setA_reference.csv", dtype={"ModelID": str})
SA = {}
for key, (q, comp) in {"Q3a_primary": ("Q3a", B0), "Q3b_primary": ("Q3b", COMP["B1_Eye"])}.items():
    rows = {}
    for ln in SETA:
        o = YA.loc[ln].dropna()
        U = pd.Index(sorted(set(universe(o)) & set(in3.index[in3])))
        ev = make_eval(dict(obs=o), q, U=U)
        pf = folds(T, f"line:{ln}", "WT")
        v = entry_skill(ev, q, pf, comp)
        dist = ref[f"{q}_fold"].to_numpy(float)
        dist = dist[np.isfinite(dist)]
        rows[ln] = dict(point=v, n_U=len(U), pct=float(100 * np.mean(dist <= v)))
        if DO_BOOT:
            rows[ln].update(boot_background({ln: (ln, ev)}, q, {ln: pf}, {ln: comp}, list(U), [ln]))
    SA[key] = dict(per_line=rows, median_pct=float(np.median([r["pct"] for r in rows.values()])))
OUT["setA"] = SA

(HERE / ("repro_out_boot.json" if DO_BOOT else "repro_out.json")).write_text(json.dumps(OUT, indent=1, default=float))
print("written")
