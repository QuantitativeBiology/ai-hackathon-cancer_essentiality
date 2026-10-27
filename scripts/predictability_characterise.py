"""Characterise the per-gene predictability tiers on dev data (DESIGN §4, D22; methods: configs/predictability.json).

    python scripts/predictability_characterise.py --threads 4

Dev lines only; CRISPR values only via essml.guard (load_targets with dev IDs; dev-anchor external profiles
via guard.load_harmonised(anchor_model_id=...), which logs each read to the tracked access log).
Compares tiers high / medium / low on: gene effect (SD, mean, fraction < -0.5), gene class, screen
reliability (cross-library lineage-mean agreement; dev-anchor agreement), own-expression association,
paralog status; reports the committed biomarker-dependency pairs; Hallmark / Reactome enrichment.

Outputs: results/predictability/characterisation.json, gene_properties.csv.gz, biomarker_pairs.csv,
enrichment.csv.
"""
from __future__ import annotations

import os
import sys

if __name__ == "__main__":
    _n = sys.argv[sys.argv.index("--threads") + 1] if "--threads" in sys.argv else "4"
    for _v in ["OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"]:
        os.environ[_v] = _n

import argparse  # noqa: E402
import datetime as dt  # noqa: E402
import io  # noqa: E402
import json  # noqa: E402
import zipfile  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
from essml import genesets, guard  # noqa: E402
from essml import metrics as M  # noqa: E402
from essml import predictability as PR  # noqa: E402
from predictability_tiers import clean, describe, git_head  # noqa: E402

OUT = ROOT / "results" / "predictability"
SEED = 20261026
N_BOOT = 2000
TIERS3 = list(PR.EVALUABLE_TIERS)


# ============================================================================ statistics
def masked_pearson_cols(A, B, min_n):
    """Column-wise Pearson r over rows where both A and B are finite; NaN if < min_n such rows."""
    ok = np.isfinite(A) & np.isfinite(B)
    n = ok.sum(0)
    A0, B0 = np.where(ok, A, 0.0), np.where(ok, B, 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        ma, mb = A0.sum(0) / n, B0.sum(0) / n
        da, db = np.where(ok, A - ma, 0.0), np.where(ok, B - mb, 0.0)
        r = (da * db).sum(0) / np.sqrt((da ** 2).sum(0) * (db ** 2).sum(0))
    r[n < min_n] = np.nan
    return r, n


def compare_continuous(x: pd.Series, tier: pd.Series, mean_r: pd.Series, rng) -> dict:
    """Per-tier summary; high vs low and medium vs low: median difference (gene-bootstrap CI), rank-biserial,
    Mann-Whitney p; Kruskal-Wallis; Spearman with mean_r (gene-bootstrap CI)."""
    x = x.astype(np.float64)
    out = {"per_tier": {t: describe(x[tier == t]) for t in TIERS3}}
    groups = {t: x[(tier == t) & np.isfinite(x)].to_numpy() for t in TIERS3}
    for a, b in (("high", "low"), ("medium", "low"), ("high", "medium")):
        ga, gb = groups[a], groups[b]
        if len(ga) < 3 or len(gb) < 3:
            continue
        u = stats.mannwhitneyu(ga, gb, alternative="two-sided")
        rb = 2 * u.statistic / (len(ga) * len(gb)) - 1
        d = float(np.median(ga) - np.median(gb))
        bs = np.array([np.median(rng.choice(ga, len(ga))) - np.median(rng.choice(gb, len(gb)))
                       for _ in range(N_BOOT)])
        out[f"{a}_vs_{b}"] = {"median_diff": d, "median_diff_ci95": [float(np.quantile(bs, 0.025)),
                                                                     float(np.quantile(bs, 0.975))],
                              "rank_biserial": float(rb), "mannwhitney_p": float(u.pvalue)}
    if all(len(g) >= 3 for g in groups.values()):
        out["kruskal_p"] = float(stats.kruskal(*groups.values()).pvalue)
    ok = np.isfinite(x) & np.isfinite(mean_r) & tier.isin(TIERS3)
    xv, mv = x[ok].to_numpy(), mean_r[ok].to_numpy()
    rho = stats.spearmanr(xv, mv).statistic
    idx = np.arange(len(xv))
    bs = []
    for _ in range(N_BOOT // 4):  # 500 resamples: Spearman on ~18k genes
        i = rng.choice(idx, len(idx))
        bs.append(stats.spearmanr(xv[i], mv[i]).statistic)
    out["spearman_with_mean_r"] = {"rho": float(rho), "ci95": [float(np.quantile(bs, 0.025)),
                                                               float(np.quantile(bs, 0.975))],
                                   "n": int(ok.sum()), "n_boot": len(bs)}
    return out


def compare_binary(x: pd.Series, tier: pd.Series) -> dict:
    x = x.astype(bool)
    out = {"per_tier": {t: {"n": int((tier == t).sum()), "n_true": int(x[tier == t].sum()),
                            "fraction": float(x[tier == t].mean()) if (tier == t).any() else None}
                        for t in TIERS3}}
    for a, b in (("high", "low"), ("medium", "low")):
        k1, n1 = int(x[tier == a].sum()), int((tier == a).sum())
        k0, n0 = int(x[tier == b].sum()), int((tier == b).sum())
        tab = np.array([[k1, n1 - k1], [k0, n0 - k0]])
        _, p = stats.fisher_exact(tab)
        h = tab + 0.5  # Haldane correction for the Woolf interval
        lor = np.log(h[0, 0] * h[1, 1] / (h[0, 1] * h[1, 0]))
        se = np.sqrt((1 / h).sum())
        out[f"{a}_vs_{b}"] = {"odds_ratio": float(np.exp(lor)),
                              "or_ci95": [float(np.exp(lor - 1.96 * se)), float(np.exp(lor + 1.96 * se))],
                              "fisher_p": float(p)}
    return out


def bh(p):
    p = np.asarray(p, np.float64)
    n = len(p)
    o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1.0)
    return out


# ============================================================================ properties
def lineage_mean_reliability(Y, lineage, lib, groups, min_lines, min_obs_mean, min_lineages):
    """configs/predictability.json characterisation.screen_reliability (primary / secondary)."""
    g1, g2 = groups
    ct = pd.crosstab(lineage, lib)
    elig = [l for l in ct.index if ct.loc[l].get(g1, 0) >= min_lines and ct.loc[l].get(g2, 0) >= min_lines]
    mats = []
    for g in (g1, g2):
        rows = []
        for l in elig:
            ids = lineage.index[(lineage == l) & (lib == g)]
            V = Y.loc[ids].to_numpy(np.float64)
            n = np.isfinite(V).sum(0)
            with np.errstate(invalid="ignore", divide="ignore"):
                m = np.where(n > 0, np.nansum(V, 0) / np.maximum(n, 1), np.nan)
            rows.append(np.where(n >= min_obs_mean, m, np.nan))
        mats.append(np.vstack(rows))
    r, n = masked_pearson_cols(mats[0], mats[1], min_lineages)
    return pd.Series(r, index=Y.columns), pd.Series(n, index=Y.columns), elig


def anchor_agreement(Y_dev_loader, min_screens, log=print):
    """Per-gene Pearson across dev anchor screens of external (LOAO-corrected) vs DepMap effect."""
    a = pd.read_csv(ROOT / "results" / "harmonisation" / "anchor_agreement.csv")
    a = a[a["role"] == "anchor"][["id", "depmap_id"]].drop_duplicates()
    roles = guard.roles_of(sorted(set(a["depmap_id"])))
    a = a[a["depmap_id"].map(roles.to_dict()).eq("dev").to_numpy()]  # dev anchors only
    ext = {}
    for aid, mid in zip(a["id"], a["depmap_id"]):
        name = f"anchor_LOAO_location__{aid}"
        if not (guard.HARMONISED_DIR / f"{name}.parquet").exists():
            log(f"  anchor profile {name} missing; skipped")
            continue
        df = guard.load_harmonised(name, anchor_model_id=mid)
        s = df["loao_corrected"].astype(np.float64)
        s.index = s.index.astype(str)
        ext[aid] = (s, mid)
    ids = sorted({m for _, m in ext.values()})
    Yd = Y_dev_loader(ids)
    genes = Yd.columns
    E = np.vstack([s.reindex(genes).to_numpy() for s, _ in ext.values()])
    D = np.vstack([Yd.loc[m].to_numpy(np.float64) for _, m in ext.values()])
    r, n = masked_pearson_cols(E, D, min_screens)
    return pd.Series(r, index=genes), pd.Series(n, index=genes), {"n_screens": len(ext), "n_lines": len(ids),
                                                                  "screens": list(ext)}


def hgnc_symbol_map():
    h = pd.read_csv(genesets.HGNC_PATH, sep="\t", dtype=str,
                    usecols=["symbol", "entrez_id", "prev_symbol", "alias_symbol"])
    h = h[h["entrez_id"].notna()]
    appr = dict(zip(h["symbol"], h["entrez_id"]))
    other = {}
    for col in ("prev_symbol", "alias_symbol"):
        for syms, e in zip(h[col], h["entrez_id"]):
            if isinstance(syms, str):
                for s in syms.split("|"):
                    other.setdefault(s.strip(), set()).add(e)
    fb = {s: next(iter(v)) for s, v in other.items() if len(v) == 1 and s not in appr}
    return appr, fb


def load_reactome(depmap_sym2ent):
    z = zipfile.ZipFile(ROOT / "data" / "raw_external" / "reactome_gmt" / "ReactomePathways.gmt.zip")
    txt = z.read(z.namelist()[0]).decode("utf-8")
    appr, fb = hgnc_symbol_map()
    sets, n_sym, n_mapped = {}, 0, 0
    for line in txt.splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        name = f"{parts[0]} ({parts[1]})"
        ents = set()
        for s in parts[2:]:
            s = s.strip()
            if not s:
                continue
            n_sym += 1
            e = depmap_sym2ent.get(s) or appr.get(s) or fb.get(s)
            if e:
                n_mapped += 1
                ents.add(str(e))
        sets[name] = ents
    return sets, {"n_sets": len(sets), "n_symbol_entries": n_sym, "n_mapped": n_mapped,
                  "fraction_mapped": n_mapped / max(n_sym, 1)}


def enrichment(fg: set, bg: set, sets: dict, lo=10, hi=500) -> pd.DataFrame:
    rows = []
    N, K = len(bg), len(fg)
    for name, s in sets.items():
        s = s & bg
        n = len(s)
        if n < lo or n > hi:
            continue
        k = len(s & fg)
        p = stats.hypergeom.sf(k - 1, N, n, K)
        a, b, c, d = k, n - k, K - k, N - n - K + k
        orr = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
        rows.append({"set": name, "set_size_in_background": n, "overlap": k, "expected": n * K / N,
                     "fold_enrichment": (k / (n * K / N)) if K else np.nan, "odds_ratio": orr, "p": p})
    df = pd.DataFrame(rows)
    if len(df):
        df["q_bh"] = bh(df["p"].to_numpy())
        df = df.sort_values(["q_bh", "p", "fold_enrichment"], ascending=[True, True, False])
    return df


# ============================================================================ main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--table", default=str(ROOT / "configs" / "predictable_genes.csv"))
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--skip-anchors", action="store_true")
    a = ap.parse_args(argv)
    cfg = PR.load_config()["characterisation"]
    rng = np.random.default_rng(SEED)
    OUT.mkdir(parents=True, exist_ok=True)

    tab = PR.load_table(a.table)
    tier = tab["tier"]
    mean_r = tab["mean_r"].astype(np.float64)
    sp = pd.read_csv(ROOT / "configs" / "splits.csv", index_col=0)
    sp.index = sp.index.astype(str)
    dev = sp[sp["role"] == "dev"]
    lineage = dev["OncotreeLineage"].astype(str)
    lib = pd.read_csv(ROOT / "configs" / "library_groups.csv", index_col=0)["lib_group"].astype(int)
    lib.index = lib.index.astype(str)
    lib = lib.reindex(dev.index)
    Y = guard.load_targets(list(dev.index)).astype(np.float64)   # dev IDs only
    Y.columns = Y.columns.astype(str)
    genes = pd.Index(tab.index)
    gc_raw = pd.read_csv(ROOT / "configs" / "gene_classes.csv", dtype={"entrez": str}).set_index("entrez")
    gc = genesets.load_gene_classes()
    gc.index = gc.index.astype(str)

    props = pd.DataFrame(index=genes)
    props["tier"] = tier
    props["mean_r"] = mean_r
    props["confidence"] = tab["confidence"]
    # gene effect (all dev lines)
    props["effect_sd"] = gc_raw["sd"].reindex(genes).astype(np.float64)
    props["effect_frac_below_-0.5"] = gc_raw["frac_below"].reindex(genes).astype(np.float64)
    props["effect_mean"] = Y.mean(axis=0).reindex(genes)
    for c in ["variable", "selective", "rare_selective", "common_essential", "nonessential"]:
        props[f"class_{c}"] = gc[c].reindex(genes).fillna(False).astype(bool)
    # own expression
    E = pd.read_parquet(ROOT / "cache" / "features" / "E.parquet")
    E.columns = E.columns.astype(str)
    E.index = E.index.astype(str)
    have = genes.intersection(E.columns).intersection(Y.columns)
    lines = dev.index.intersection(E.index)
    r_own = M.per_gene_pearson(E.loc[lines, have].to_numpy(np.float64), Y.loc[lines, have].to_numpy(np.float64), 40)
    props["own_expr_r"] = pd.Series(r_own, index=have).reindex(genes)
    props["own_expr_abs_r"] = props["own_expr_r"].abs()
    props["own_expr_mean_log2tpm"] = E.loc[lines].mean(axis=0).reindex(genes)
    del E
    # paralogs
    pt = pd.read_parquet(ROOT / "cache" / "priors" / "paralog_top.parquet")
    pt.index = pt.index.astype(str)
    props["paralog_has_id20"] = (pt["family_size_id20"].reindex(genes).fillna(0) > 0)
    props["paralog_max_identity"] = pt["max_identity"].reindex(genes).fillna(0.0).astype(np.float64)
    props["paralog_high_identity_ge50"] = props["paralog_max_identity"] >= 50
    # screen reliability (declared methods)
    rel = cfg["screen_reliability"]
    rel_info = {}
    for key in ("primary", "secondary"):
        spec = rel[key]
        r, n, elig = lineage_mean_reliability(Y, lineage, lib, spec["groups"], spec["min_lines_per_lineage_per_group"],
                                              spec["min_observed_lines_for_a_mean"], spec["min_lineages"])
        props[spec["name"]] = r.reindex(genes)
        props[spec["name"] + "_n_lineages"] = n.reindex(genes)
        rel_info[spec["name"]] = {"groups": spec["groups"], "eligible_lineages": elig,
                                  "n_genes_defined": int(r.notna().sum())}
    if not a.skip_anchors:
        def loader(ids):
            Yd = guard.load_targets(ids)
            Yd.columns = Yd.columns.astype(str)
            return Yd.astype(np.float64)
        ra, na, ainfo = anchor_agreement(loader, rel["anchor_agreement"]["min_anchor_screens"])
        props["external_anchor_agreement"] = ra.reindex(genes)
        props["external_anchor_n_screens"] = na.reindex(genes)
        rel_info["external_anchor_agreement"] = {**ainfo, "n_genes_defined": int(ra.notna().sum())}
    del Y

    ev = tier.isin(TIERS3)
    res = {"generated": dt.datetime.now().isoformat(timespec="seconds"), "git_commit": git_head(),
           "table": Path(a.table).relative_to(ROOT).as_posix(), "table_sha256": PR.sha256_file(a.table),
           "population": cfg["population"], "comparisons_method": cfg["comparisons"], "reliability_info": rel_info,
           "n_per_tier": tier.value_counts().to_dict(), "continuous": {}, "binary": {}}
    cont = ["effect_sd", "effect_mean", "effect_frac_below_-0.5", "own_expr_r", "own_expr_abs_r",
            "own_expr_mean_log2tpm", "paralog_max_identity", rel["primary"]["name"], rel["secondary"]["name"]]
    if "external_anchor_agreement" in props:
        cont.append("external_anchor_agreement")
    for c in cont:
        print(f"  comparing {c}")
        res["continuous"][c] = compare_continuous(props.loc[ev, c], tier[ev], mean_r[ev], rng)
    # reliability within variable genes as well (most reliability signal is in variable genes)
    var = props["class_variable"]
    for c in [rel["primary"]["name"], rel["secondary"]["name"]] + (
            ["external_anchor_agreement"] if "external_anchor_agreement" in props else []) + ["own_expr_abs_r"]:
        res["continuous"][c + "__variable_genes"] = compare_continuous(props.loc[ev & var, c], tier[ev & var],
                                                                       mean_r[ev & var], rng)
    for c in ["class_variable", "class_selective", "class_rare_selective", "class_common_essential",
              "class_nonessential", "paralog_has_id20", "paralog_high_identity_ge50"]:
        res["binary"][c] = compare_binary(props.loc[ev, c], tier[ev])
    res["binary"]["paralog_has_id20__variable_genes"] = compare_binary(props.loc[ev & var, "paralog_has_id20"],
                                                                      tier[ev & var])
    # gene class x tier table (all tiers incl. not evaluable)
    res["class_by_tier"] = {c: props.loc[props[f"class_{c}"], "tier"].value_counts().reindex(PR.TIERS, fill_value=0)
                            .astype(int).to_dict()
                            for c in ["variable", "selective", "rare_selective", "common_essential", "nonessential"]}
    res["class_by_tier"]["all"] = tier.value_counts().reindex(PR.TIERS, fill_value=0).astype(int).to_dict()

    # biomarker pairs (list committed before this comparison)
    bp = pd.read_csv(ROOT / "configs" / "biomarker_pairs.csv")
    sym2ent = {s: e for e, s in gc_raw["symbol"].items()}
    bp["dependency_entrez"] = bp["dependency"].map(sym2ent)
    bp = bp.join(props[["tier", "mean_r", "confidence", "effect_sd", "effect_frac_below_-0.5", "own_expr_r"]],
                 on="dependency_entrez")
    bp = bp.join(tab[[f"r_fold{k}" for k in range(5)] + ["min_r"]], on="dependency_entrez")
    bp.to_csv(OUT / "biomarker_pairs.csv", index=False, lineterminator="\n")
    res["biomarker_pairs"] = {
        "file": "results/predictability/biomarker_pairs.csv",
        "tier_counts": bp["tier"].value_counts().to_dict(),
        "tier_by_relation": pd.crosstab(bp["relation"], bp["tier"]).to_dict(orient="index"),
        "tier_by_headline_input": pd.crosstab(bp["headline_input"], bp["tier"]).to_dict(orient="index"),
        "median_mean_r_by_input": bp.groupby("headline_input")["mean_r"].median().to_dict(),
        "fraction_high_overall_evaluable": float((tier == "high").sum() / ev.sum())}

    # enrichment
    import json as _j
    hall = {k: set(map(str, v)) for k, v in _j.loads((ROOT / "cache" / "features" / "gene_sets.json").read_text())[
        "hallmark"].items()}
    reac, reac_cov = load_reactome(sym2ent)
    bg_all = set(genes[ev])
    bg_var = set(genes[ev & var])
    tests = {"high_vs_evaluable": (set(genes[tier == "high"]), bg_all),
             "low_among_variable": (set(genes[(tier == "low") & var]), bg_var)}
    enr_rows = []
    res["enrichment"] = {"reactome_mapping": reac_cov, "hallmark_n_sets": len(hall), "tests": {}}
    for tname, (fg, bg) in tests.items():
        for cname, sets in (("hallmark", hall), ("reactome", reac)):
            df = enrichment(fg, bg, sets)
            df.insert(0, "collection", cname)
            df.insert(0, "test", tname)
            enr_rows.append(df)
            sig = df[df["q_bh"] < 0.05]
            res["enrichment"]["tests"][f"{tname}__{cname}"] = {
                "n_foreground": len(fg), "n_background": len(bg), "n_sets_tested": int(len(df)),
                "n_sig_q05": int(len(sig)),
                "top": sig.head(15)[["set", "set_size_in_background", "overlap", "fold_enrichment", "q_bh"]]
                .to_dict(orient="records")}
    enr = pd.concat(enr_rows, ignore_index=True)
    enr.to_csv(OUT / "enrichment.csv", index=False, lineterminator="\n", float_format="%.6g")

    props.round(6).rename_axis("entrez").to_csv(OUT / "gene_properties.csv.gz", lineterminator="\n")
    (OUT / "characterisation.json").write_text(json.dumps(clean(res), indent=1))
    print(f"wrote {OUT.relative_to(ROOT)}/characterisation.json")
    return res


if __name__ == "__main__":
    main()
