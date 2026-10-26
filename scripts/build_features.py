"""Build and cache raw (UNFITTED) line-feature blocks and gene-level prior tables.

Outputs (git-ignored):
  cache/features/<block>.parquet   E, C, C_instab, M, M_burden, H, F_pair, F_gene, meta, meta_prov,
                                   composite, P   (rows ModelID; see essml.features docstring)
  cache/features/gene_sets.json    {'hallmark': {set: [entrez]}, 'p53_targets': {...}}
  cache/priors/paralog_pairs.parquet, paralog_top.parquet, partners.parquet
Summary (committed): results/features_build_report.json — shapes, line coverage per modality and
split role, prior coverage over CRISPR genes (gene NAMES only, from cache/gene_map.parquet), the
fusion-filter counts, proteomics matching and the deep-deletion threshold evidence.

Usage: python scripts/build_features.py
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from essml import features, priors  # noqa: E402
from essml.io import load_gene_map, load_matrix  # noqa: E402

FEAT_DIR = ROOT / "cache" / "features"
PRI_DIR = ROOT / "cache" / "priors"
REPORT = ROOT / "results" / "features_build_report.json"


def deep_del_evidence():
    """Density of relative CN by bin and fraction not expressed (TPM<1), lines/genes in CN∩E."""
    cn, ex = load_matrix("cn"), load_matrix("expression")
    lines, genes = cn.index.intersection(ex.index), cn.columns.intersection(ex.columns)
    c = cn.loc[lines, genes].to_numpy()
    e = ex.loc[lines, genes].to_numpy()
    del cn, ex
    edges = [0, .05, .1, .15, .2, .25, .3, .35, .4, .45, .5, .75, 1.0]
    ok = ~np.isnan(c)
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = ok & (c >= lo) & (c < hi)
        n = int(sel.sum())
        rows.append({"bin": f"[{lo:.2f},{hi:.2f})", "n": n, "n_per_0.05": round(n / ((hi - lo) / .05)),
                     "frac_not_expressed": round(float((e[sel] < 1).mean()), 3) if n else None})
    v = c[ok]
    return {"bins": rows, "frac_values_below_thr": float((v < priors.DEEP_DEL_THR).mean()),
            "frac_of_sub_thr_below_0.05": float((v < .05).sum() / (v < priors.DEEP_DEL_THR).sum()),
            "median_deep_del_genes_per_line": float(np.median((np.where(ok, c, 1) < priors.DEEP_DEL_THR).sum(1)))}


def main():
    t0 = time.time()
    rep = {}
    splits = pd.read_csv(ROOT / "configs" / "splits.csv", index_col=0)
    hgnc = priors.load_hgnc()

    # ---------------- raw line-feature blocks
    blocks = features.build_all(report=rep)
    features.save_raw(blocks, FEAT_DIR)
    gs = {"hallmark": features.hallmark_sets(hgnc), "p53_targets": features.p53_target_set(hgnc)}
    (FEAT_DIR / "gene_sets.json").write_text(json.dumps(gs))
    rep["gene_sets"] = {"hallmark_n_sets": len(gs["hallmark"]),
                        "hallmark_genes_median": int(np.median([len(v) for v in gs["hallmark"].values()])),
                        "hallmark_genes_in_E": int(len({g for v in gs["hallmark"].values() for g in v}
                                                     & set(blocks["E"].columns))),
                        "p53_targets": gs["p53_targets"]["p53_targets"]}
    shapes = {}
    for k, v in blocks.items():
        cov = {r: int(v.index.isin(splits.index[splits["role"] == r]).sum()) for r in ["dev", "test", "noncancerous"]}
        shapes[k] = {"rows": int(v.shape[0]), "cols": int(v.shape[1]), **{f"lines_{r}": n for r, n in cov.items()}}
        print(f"{k:10s} {v.shape[0]:5d} x {v.shape[1]:6d}   dev {cov['dev']:4d}  test {cov['test']:3d}  "
              f"noncancerous {cov['noncancerous']:2d}")
    rep["blocks"] = shapes
    comp = blocks["composite"]
    rep["composite_prevalence_dev"] = comp.reindex(splits.index[splits["role"] == "dev"]).mean().round(3).to_dict()
    rep["C_genes_dropped_all_nan"] = int(len(load_gene_map().query("in_matrices.str.contains('cn')"))
                                         - blocks["C"].shape[1])
    del blocks

    # ---------------- priors
    pr = {}
    pairs = priors.load_paralog_pairs(hgnc=hgnc)
    top = priors.top_paralogs(pairs)
    within = priors.top_paralogs(pairs[pairs["type"] == "within_species_paralog"])
    PRI_DIR.mkdir(parents=True, exist_ok=True)
    pairs.to_parquet(PRI_DIR / "paralog_pairs.parquet", index=False)
    top.to_parquet(PRI_DIR / "paralog_top.parquet")
    srcs = [priors.load_bioplex(), priors.load_huri(hgnc=hgnc), priors.load_biogrid(report=pr)]
    partners = priors.partner_table(srcs)
    partners.to_parquet(PRI_DIR / "partners.parquet", index=False)
    pr["paralog_pairs_directed"] = int(len(pairs))
    pr["paralog_pairs_by_type"] = pairs["type"].value_counts().to_dict()
    pr["genes_with_paralog"] = int(len(top))
    pr["interaction_pairs_undirected"] = {s: int(len(d) // 2) for s, d in partners.groupby("source")}

    gm = load_gene_map()
    crispr_genes = gm.loc[gm["in_matrices"].str.split(",").apply(lambda x: "crispr" in x), "entrez"]
    cov = priors.coverage(crispr_genes, top, partners)
    cov["paralog_within_species_only"] = float(crispr_genes.isin(within.index).mean())
    rep["priors"] = pr
    rep["prior_coverage_crispr_genes"] = cov
    print("prior coverage over CRISPR genes:", json.dumps(cov, indent=1))

    # ---------------- local-feature smoke test on a few lines/genes
    b = priors.LocalFeatureBuilder.from_cache()
    ids = splits.index[splits["role"] == "dev"][:5]
    lf = b.local_features(ids, crispr_genes.iloc[:200])
    rep["local_features_smoke"] = {k: [list(v.shape), float(np.nanmean(v))] for k, v in lf.items()}
    del b, lf

    rep["deep_deletion_evidence"] = deep_del_evidence()
    print("deep-deletion evidence:", json.dumps(rep["deep_deletion_evidence"], indent=1))
    print("fusion filter:", rep["fusion_filter"])
    print("proteomics:", {k: v for k, v in rep["proteomics"].items() if k != "unmatched_names"})
    rep["runtime_s"] = round(time.time() - t0, 1)
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(rep, indent=1, default=str))
    print(f"done in {rep['runtime_s']} s -> {REPORT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
