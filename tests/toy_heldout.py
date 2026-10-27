"""Toy repository for the G4/G6 evaluation-script tests (pattern of tests/test_guard.py).

Builds, under a pytest tmp_path, a git repository with its own splits.csv, CRISPR parquet, feature
matrix, gene lists, configs, finalists.json and a TRACKED access log, so the real guard runs its
full unblinding logic against toy data. Nothing under the real repository is read except the
Hart NEGv1 reference list (genesets.load_hart_neg, used by the launcher's inner tuning).
"""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

import heldout_common as HC
from essml import genesets

EC = HC.EC_DEV
N_GENES = 1200
N_FEAT = 40


def git(repo, *a, check=True):
    return subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t", *a], cwd=repo, check=check,
                          capture_output=True, text=True)


def _ids(prefix, n):
    return [f"{prefix}{i:04d}" for i in range(n)]


def make_world(tmp_path, seed=0, n_dev=160, n_test=60, n_genes=N_GENES, signal=1.0):
    rng = np.random.default_rng(seed)
    repo = Path(tmp_path) / "repo"
    (repo / "configs" / "runs").mkdir(parents=True)
    (repo / "results" / "gates").mkdir(parents=True)
    (repo / "cache").mkdir()
    dev = EC + _ids("D", n_dev - len(EC))
    test = _ids("T", n_test)
    nonc = _ids("N", 5)
    ids = dev + test + nonc
    lineages = np.array(["L0", "L1", "L2", "L3"])
    lin = rng.choice(lineages, len(ids)).astype(object)
    lin[len(dev) + len(test):] = "Eye"
    lin[: len(EC)] = "L3"
    pat = np.array([f"p{i}" for i in range(len(ids))], dtype=object)
    pat[1] = pat[0]                                      # shared patient inside dev (same fold below)
    pat[len(dev) + 1] = pat[len(dev)]                    # shared patient inside test
    folds = np.full(len(ids), -1)
    folds[: len(dev)] = np.arange(len(dev)) % 5
    folds[1] = folds[0]
    role = np.array(["dev"] * len(dev) + ["test"] * len(test) + ["noncancerous"] * len(nonc))
    splits = pd.DataFrame({"ModelID": ids, "PatientID": pat, "OncotreeLineage": lin, "role": role, "cv_fold": folds})
    splits.to_csv(repo / "configs" / "splits.csv", index=False)

    genes = [str(9_000_000 + g) for g in range(n_genes)]
    X = rng.normal(size=(len(ids), N_FEAT))
    B = np.zeros((N_FEAT, n_genes))
    act = rng.random(n_genes) < 0.5
    B[: 8, act] = rng.normal(0, 0.25, (8, int(act.sum())))
    mu = rng.normal(-0.1, 0.2, n_genes)
    mu[:150] = rng.normal(-1.6, 0.3, 150)               # common essentials
    lin_eff = {l: rng.normal(0, 0.15, n_genes) for l in list(lineages) + ["Eye"]}
    Y = mu + np.stack([lin_eff[l] for l in lin]) + signal * (X @ B) + rng.normal(0, 0.25, (len(ids), n_genes))
    nan_lines = rng.random(len(ids)) < 0.2               # a second "library": 5% of genes unobserved
    nan_genes = rng.random(n_genes) < 0.05
    Y[np.ix_(nan_lines, nan_genes)] = np.nan
    Ydf = pd.DataFrame(Y.astype(np.float32), index=pd.Index(ids, name="ModelID"), columns=genes)
    Ydf.to_parquet(repo / "cache" / "crispr.parquet")
    feat = pd.DataFrame(X, index=ids, columns=[f"e{j}" for j in range(N_FEAT)])

    # gene lists from dev lines only (as configs/gene_classes.csv / eval_genes.csv)
    Ydev = Ydf.loc[dev].astype(np.float64)
    gc = genesets.gene_classes(Ydev, neg_entrez=genes[150:330])
    v = Ydev.to_numpy()
    obs = ~np.isnan(v)
    gc_out = gc.copy()
    gc_out.insert(0, "symbol", [f"SYM{g}" for g in genes])
    gc_out["n_obs"] = obs.sum(0)
    gc_out["sd"] = np.nanstd(v, 0, ddof=1)
    gc_out["n_below"] = ((v < -0.5) & obs).sum(0)
    gc_out["frac_below"] = gc_out["n_below"] / gc_out["n_obs"]
    gc_out.index.name = "entrez"
    gc_out.reset_index().to_csv(repo / "configs" / "gene_classes.csv", index=False)
    pd.DataFrame({"entrez": genes, "symbol": [f"SYM{g}" for g in genes]}).to_csv(repo / "configs" / "eval_genes.csv", index=False)
    cs = {"ceg_scaling": genes[0:40], "ceg_evaluation": genes[40:80], "neg_scaling": genes[150:200],
          "neg_evaluation": genes[200:330], "scaling_set": genes[80:120]}
    (repo / "configs" / "control_split.json").write_text(json.dumps(cs))
    from make_library_groups import library_groups
    lg = library_groups(Ydf.loc[dev])
    pd.DataFrame({"lib_group": lg, "n_nan_genes": Ydf.loc[dev].isna().sum(1)}).rename_axis("ModelID").to_csv(
        repo / "configs" / "library_groups.csv")
    (repo / "configs" / "seeds.json").write_text(json.dumps({"split_seed": 20261026, "model_seeds": [0]}))
    (repo / "configs" / "frozen_hashes.json").write_text(json.dumps({"files": {}}))
    (repo / "results" / "gates" / "G3.json").write_text(json.dumps({"gate": "G3", "pass": True}))
    (repo / "results" / "gates" / "G4.json").write_text(json.dumps({"gate": "G4", "pass": True}))
    pub = genes[n_genes // 2: n_genes // 2 + 2 * (n_genes // 8)]
    pd.DataFrame({"set": ["deepdep_1298"] * (len(pub) // 2) + ["leap_1539"] * (len(pub) - len(pub) // 2),
                  "symbol": "x", "entrez": pub, "mapped_via": "toy"}).to_csv(repo / "configs" / "published_genesets.csv", index=False)
    (repo / "results" / "test_access_log.csv").write_text(",".join(["timestamp", "caller", "access_type", "scope",
                                                                     "config_hash", "n_ids", "n_genes", "roles",
                                                                     "granted", "reason"]) + "\n")
    (repo / ".gitignore").write_text("cache/\n*.parquet\n*.pkl\n")
    # D22 per-gene predictability table (toy: genes with planted expression signal tend to be 'high')
    mr = np.where(act, rng.uniform(0.1, 0.7, n_genes), rng.uniform(-0.05, 0.25, n_genes))
    sd = rng.uniform(0.02, 0.08, n_genes)
    tier = np.where(mr >= 0.3, "high", np.where(mr >= 0.15, "medium", "low")).astype(object)
    tier[-10:] = "not_evaluable"
    pgt = pd.DataFrame({"entrez": genes, "symbol": [f"SYM{g}" for g in genes], "mean_r": mr, "min_r": mr - 2 * sd,
                        "sd_r": sd, "confidence": np.clip(mr - sd, -1, 1), "tier": tier, "source_run_id": "toy_cv",
                        "threshold_high_mean_r": 0.3, "folds_used": "0;1;2;3;4", "table_role": "primary"})
    pgt.loc[pgt["tier"] == "not_evaluable", ["mean_r", "min_r", "sd_r", "confidence"]] = np.nan
    pgt.to_csv(repo / "configs" / "predictable_genes.csv", index=False)
    # the additional finalist's own secondary tier table (D22f), shifted so its partition differs
    own = pgt.assign(tier=np.roll(pgt["tier"].to_numpy(), 7), source_run_id="toy_cv_add",
                     table_role="secondary (toy additional finalist)")
    own.to_csv(repo / "configs" / "predictable_genes_add.csv", index=False)

    kappa = {"grid": {"kappa": [0, 10, 100]}}
    base = {"reference_baselines": {"B0": {}, "B1": kappa}, "inner_folds": 3, "min_obs": 20, "inner_min_obs": 10,
            "threads": 2, "requires_gate": "G0"}
    cfgs = {
        "headline": {"name": "toy_m1", "model": "M1", "features": ["toyE"], "params": {"residual": "B0"}, "grid": {}},
        "additional": {"name": "toy_r1", "model": "R1", "features": ["toyE"], "params": {}, "grid": {}},
        "T": {"name": "toy_t", "model": "T", "features": ["toyE"], "params": {}, "grid": {}},
        "H3_noEC": {"name": "toy_t_noEC", "model": "T", "features": ["toyE"], "params": {}, "grid": {},
                    "exclude_from_training": list(EC)},
        "B0": {"name": "toy_b0", "model": "B0", "features": [], "params": {}, "grid": {}},
        "B1": {"name": "toy_b1", "model": "B1", "features": [], "params": {}, "grid": kappa["grid"]},
        "B2": {"name": "toy_b2", "model": "B2", "features": ["toyE"], "params": {"n_top": 30}, "grid": {"k": [5, 10]}},
    }
    fin = {}
    for role, c in cfgs.items():
        c = {**c, **base}
        p = repo / "configs" / "runs" / f"{c['name']}.json"
        p.write_text(json.dumps(c, indent=1))
        fin[role] = {"name": c["name"], "config": f"configs/runs/{c['name']}.json", "config_hash": HC.config_hash(c)}
    # a registered toy CV run of T (independent run_cv), so G4's fold-model reproduction check has a reference
    cfg_t = {k: v for k, v in json.loads((repo / "configs" / "runs" / "toy_t.json").read_text()).items()
             if k != "reference_baselines"}
    sp_dev = splits.set_index("ModelID").loc[dev]
    neg, _ = genesets.load_hart_neg()
    cv = HC.RUN.run_cv(Ydf.loc[dev].astype(np.float64), sp_dev["cv_fold"].astype(int), sp_dev["PatientID"].astype(str),
                       sp_dev["OncotreeLineage"].astype(str), feat.loc[dev].add_prefix("toyE:"), cfg_t,
                       pd.Index(genes), neg, 20261026, log=lambda *_: None)
    (repo / "cache" / "runs" / "toy_cv_t").mkdir(parents=True)
    cv["pred"].astype(np.float32).to_parquet(repo / "cache" / "runs" / "toy_cv_t" / "pred.parquet")
    fin["T"]["cv_run_id"] = "toy_cv_t"
    fin["additional"].update(cv_run_id="toy_cv_add", own_predictable_genes="configs/predictable_genes_add.csv")
    finalists = {"headline": fin["headline"], "additional": [fin["additional"]], "T": fin["T"], "H3_noEC": fin["H3_noEC"],
                 "baselines": {b: fin[b] for b in ("B0", "B1", "B2")}}
    (repo / "configs" / "finalists.json").write_text(json.dumps(finalists, indent=1))
    git(repo, "init", "-q")
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-m", "toy world")

    paths = HC.Paths(root=repo, splits=repo / "configs" / "splits.csv", crispr=repo / "cache" / "crispr.parquet",
                     access_log=repo / "results" / "test_access_log.csv", eval_genes=repo / "configs" / "eval_genes.csv",
                     gene_classes=repo / "configs" / "gene_classes.csv", control_split=repo / "configs" / "control_split.json",
                     library_groups=repo / "configs" / "library_groups.csv", seeds=repo / "configs" / "seeds.json",
                     gates=repo / "results" / "gates", registry=repo / "results" / "registry.csv",
                     model_dir=repo / "cache" / "models", frozen_hashes=repo / "configs" / "frozen_hashes.json",
                     criteria=repo / "configs" / "transfer_criteria.json", harmonised_dir=repo / "cache" / "harmonised",
                     comparators=repo / "cache" / "transfer" / "comparators.parquet",
                     paralog_pairs=repo / "cache" / "paralog_pairs.parquet",
                     external_expr=repo / "cache" / "external_expr.parquet",
                     published_sets=repo / "configs" / "published_genesets.csv", cv_cache=repo / "cache" / "runs",
                     predictable_genes=repo / "configs" / "predictable_genes.csv",
                     out_dir=repo / "results" / "test", pred_dir=repo / "cache" / "test" / "preds")

    def loader(names, rows):
        if not names:
            return None
        return pd.concat([feat.reindex(pd.Index(rows)).add_prefix(f"{n}:") for n in names], axis=1)

    return {"repo": repo, "paths": paths, "loader": loader, "Y": Ydf, "X": feat, "splits": splits,
            "dev": dev, "test": test, "nonc": nonc, "genes": genes, "finalists": repo / "configs" / "finalists.json",
            "cs": cs, "gc": gc_out, "mu": mu, "B": B, "lineage": pd.Series(lin, index=ids)}


# ============================================================================ G6 toy world
def ceil(cs_, mde):
    return {k: {"ceiling_skill": cs_, "MDE": mde} for k in
            ("Q3a_primary", "Q3a_effect_label", "Q3b_primary", "Q3b_vs_B0", "Q3b_rare_vs_B1n")}


def make_transfer_world(w, n_boot=100):
    """Three synthetic verdict backgrounds with known answers, written as harmonised 'screens':

    RPE1  ("transfers"): arms = mu + x* B + noise; T learns mu + x B, so it beats B0 well above 50% of the
          frozen ceiling (0.02); the swap profile carries no information.
    hPSC  ("fails"): arms = mu + noise (no expression signal: T's deviations from B0 only add noise);
          upper CI < MDE (0.3), harmonisation passed, library-matched anchor = a dev line's own toy profile
          (T's out-of-fold model predicts it well -> gate positive).
    bgI   ("inconclusive"): same data type as hPSC, but its anchor profile is pure noise (gate not positive).
    Set A: the 5 toy noncancerous lines (their toy CRISPR rows; lineage Eye -> B1_Eye = B0 here).
    """
    repo, genes, rng = w["repo"], w["genes"], np.random.default_rng(3)
    mu, B = w["mu"], w["B"]
    hdir = repo / "cache" / "harmonised"
    hdir.mkdir(parents=True)
    xs = {"prof_T": rng.normal(size=N_FEAT), "prof_swap": rng.normal(size=N_FEAT), "prof_null": rng.normal(size=N_FEAT)}

    def arm(name, values):
        v = pd.Series(values, index=genes).astype(np.float32)
        pd.DataFrame({"raw": v, "affine": v, "final": v, "affine_alt": v}).rename_axis("entrez").to_parquet(
            hdir / f"{name}.parquet")

    signal = mu + xs["prof_T"] @ B
    for a in ("t1_wt", "t1_null", "t2_a", "t2_b", "t3"):
        arm(a, signal + rng.normal(0, 0.25, len(genes)))
    for a in ("f_a", "f_b", "po", "i_a"):
        arm(a, mu + rng.normal(0, 0.25, len(genes)))
    good, bad = w["dev"][10], w["dev"][11]
    anchors = {"A_good": (w["Y"].loc[good].astype(np.float64).dropna(), good),
               "A_bad": (pd.Series(rng.normal(-0.2, 0.6, len(genes)), index=genes), bad)}
    hv = {"passed": True}

    def bg(studies, prof, b1, anchors_, swap, c):
        return {"verdict_eligible": True, "studies": studies, "expression_input_T": prof, "B1_nearest_column": b1,
                "library_matched_anchors": anchors_, "expression_swap_profiles": swap,
                "harmonisation_validation": hv, "replicate_ceiling": c}

    crit = {
        "bootstrap": {"n_resamples": n_boot, "seed": 20261026},
        "inputs_sha256": {},
        "backgrounds": {
            "RPE1": bg({"s_t1": {"role": "verdict", "entries": {"WT": {"arms": ["t1_wt"]}, "null": {"arms": ["t1_null"]}}},
                        "s_t2": {"role": "verdict", "entries": {"WT": {"arms": ["t2_a", "t2_b"]}}},
                        "s_t3": {"role": "verdict", "entries": {"null": {"arms": ["t3"]}}}},
                       {"profile_id": "prof_T"}, "B1_L0", ["A_good"], {"swap": ["prof_swap"]}, ceil(0.02, 0.05)),
            "hPSC": bg({"s_fail": {"role": "verdict", "entries": {"WT": {"arms": ["f_a", "f_b"]}}},
                        "s_po": {"role": "predictor_only", "entries": {"WT": {"arms": ["po"]}}}},
                       {"per_arm": {"f_a": "prof_null", "f_b": "prof_null", "po": "prof_null"}}, "B1_L1", ["A_good"],
                       {"swap": ["prof_T"]}, ceil(0.2, 0.3)),
            "bgI": bg({"s_inc": {"role": "verdict", "entries": {"WT": {"arms": ["i_a"]}}}},
                      {"profile_id": "prof_null"}, "B1_L2", ["A_bad"], {}, ceil(0.2, 0.3)),
            "primary_T": {"verdict_eligible": False, "status": "inconclusive by declaration", "reason": "toy"}},
        "set_A": {"tier_A1": {"lines": list(w["nonc"])}},
        "hypotheses": {"H1": {"genes": {"SYMa": genes[500], "SYMb": genes[501]}},
                       "H1_prime": {"genes": {"SYMc": genes[502]}},
                       "H2": {"genes": {f"S{i}": genes[500 + i] for i in range(12)}},
                       "H3": {"genes": {f"P{i}": genes[600 + i] for i in range(6)}}},
    }
    h2_arms = {"brown_wt_entry": "s_t1:WT", "brown_null_entry": "s_t1:null", "bowden_wt": None, "bowden_ko": None,
               "mair_entry": "s_fail:WT"}
    Ydev = w["Y"].loc[w["dev"]].astype(np.float64)
    from essml.baselines import GeneMean, LineageMean
    b0 = pd.Series(GeneMean().fit(Ydev).mu_, index=genes)
    lm = LineageMean(10).fit(Ydev, w["lineage"].reindex(Ydev.index)).lineage_means()
    comp = pd.DataFrame({"B0": b0, "B1_L0": lm.loc["L0"], "B1_L1": lm.loc["L1"], "B1_L2": lm.loc["L2"], "B1_Eye": b0})
    ref = repo / "results" / "test" / "toy_t"                # stand-in for G4's set-A reference distribution
    ref.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"ModelID": w["test"], "Q3a_fold": rng.normal(0.05, 0.03, len(w["test"])),
                  "Q3b_fold": rng.normal(0.05, 0.05, len(w["test"]))}).to_csv(ref / "setA_reference.csv", index=False)
    for f in ("harmonisation_frozen.json", "transfer_criteria.json"):
        (repo / "configs" / f).write_text(json.dumps({"toy": True}))
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-m", "toy transfer freeze")
    feats = w["X"]

    def expression(pid):
        if str(pid).startswith("line:"):
            return feats.loc[str(pid)[5:]]
        return pd.Series(xs[pid], index=feats.columns)

    def x_builder(f, model, pid, row, tp53):
        key = f"profile::{pid}::{tp53}"
        return key, pd.DataFrame([row.to_numpy()], index=[key], columns=[f"toyE:{c}" for c in row.index])

    return {"spec": {"crit": crit, "h2_arms": h2_arms, "b1_setA": "B1_Eye"}, "comparators": comp, "anchors": anchors,
            "expression": expression, "x_builder": x_builder}
