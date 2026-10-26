"""Published gene sets for the "loosely comparable" block (DESIGN §4) -> configs/published_genesets.csv.

Sources (authors' own code repositories; plain GitHub raw downloads, no bot check):
- DeepDEP 1,298 genes: chenlabgccri/Prep4DeepDEP @ f713c84, inst/extdata/default_dep_genes_1298.RData
  (data.frame dep.data$Gene; converted to CSV with R 4.4.1). sha256 of the RData:
  38e185e3bc166b683c022316768c062f10a7cd7ea1b1b6ae23dd012552c566b3
- LEAP 1,539 genes (task 1, gene essentiality): owkin/leap @ 7c709df, data/list_of_perturbations_task_1.csv
  (lower-case symbols). sha256 c349c8a89caac3fa73ab1cf8889687eda43d5c74d741fa16e01d343806c5b632
Files live in data/raw_external/published_gene_sets/ (git-ignored; DATA.md row).

Symbols -> Entrez: (1) the CRISPR matrix symbol (cache/gene_map.parquet), case-insensitive;
(2) HGNC approved symbol; (3) HGNC previous/alias symbol if it maps to exactly one HGNC entry.
"""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from essml import io  # noqa: E402

SRC = ROOT / "data" / "raw_external" / "published_gene_sets"
FILES = {
    "deepdep_1298": (SRC / "deepdep" / "default_dep_genes_1298.csv",
                     "default_dep_genes_1298.RData",
                     "38e185e3bc166b683c022316768c062f10a7cd7ea1b1b6ae23dd012552c566b3"),
    "leap_1539": (SRC / "leap" / "list_of_perturbations_task_1.csv", None,
                  "c349c8a89caac3fa73ab1cf8889687eda43d5c74d741fa16e01d343806c5b632"),
}


def _symbols(name, path):
    if name == "deepdep_1298":
        return pd.read_csv(path)["Gene"].astype(str).tolist()
    return pd.read_csv(path, header=None)[0].astype(str).str.strip().tolist()


def main():
    gm = io.load_gene_map()
    crispr = gm[gm["in_matrices"].str.contains("crispr")]
    by_sym = {s.upper(): str(e) for s, e in zip(crispr["symbol"], crispr["entrez"])}
    h = pd.read_csv(ROOT / "data" / "raw_external" / "hgnc" / "hgnc_complete_set.txt", sep="\t",
                    dtype=str, usecols=["symbol", "alias_symbol", "prev_symbol", "entrez_id"])
    h = h.dropna(subset=["entrez_id"])
    appr = {s.upper(): e for s, e in zip(h["symbol"], h["entrez_id"])}
    other: dict[str, set] = {}
    for col in ("prev_symbol", "alias_symbol"):
        for syms, e in zip(h[col], h["entrez_id"]):
            if isinstance(syms, str):
                for s in syms.split("|"):
                    other.setdefault(s.strip().upper(), set()).add(e)
    rows = []
    for name, (path, raw, sha) in FILES.items():
        rawp = path.parent / raw if raw else path
        got = hashlib.sha256(rawp.read_bytes()).hexdigest()
        if got != sha:
            raise ValueError(f"{rawp}: sha256 {got} != {sha}")
        for s in _symbols(name, path):
            u = s.upper()
            if u in by_sym:
                e, via = by_sym[u], "crispr_symbol"
            elif u in appr:
                e, via = appr[u], "hgnc_approved"
            elif len(other.get(u, ())) == 1:
                e, via = next(iter(other[u])), "hgnc_prev_alias"
            else:
                e, via = "", "unmapped"
            rows.append({"set": name, "symbol": s, "entrez": e, "mapped_via": via})
    out = pd.DataFrame(rows)
    out.to_csv(ROOT / "configs" / "published_genesets.csv", index=False)
    print(out.groupby(["set", "mapped_via"]).size())


if __name__ == "__main__":
    main()
