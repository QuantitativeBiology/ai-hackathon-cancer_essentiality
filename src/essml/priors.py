"""Gene-level prior knowledge keyed by Entrez ID strings (DESIGN.md §2 "Gene-level features").

Sources (all open, non-CRISPR, non-text-mined; checksums in configs/external_manifest.json, DATA.md):
- HGNC complete set: the ID hub. Ensembl gene -> Entrez, symbol -> Entrez, chromosome, locus group.
- Ensembl 116 Compara homologies (homo_sapiens, protein_default): human-human paralogs.
- BioPlex 3.0 293T network (Entrez already), HuRI (Ensembl gene IDs), BioGRID physical
  (Homo sapiens on both sides, Experimental System Type == 'physical'). BioGRID ORCS (CRISPR screen
  database) is a separate download and is never read; rows whose Source Database is not a
  curation database are not used either (see ``load_biogrid``).

Paralog homology types
----------------------
Ensembl splits human-human paralogs into ``within_species_paralog`` (duplication node inside the
human lineage per the gene tree) and ``other_paralog`` (older duplications, e.g. SMARCA2-SMARCA4,
identity 78%). ``within_species_paralog`` alone covers only 4,460 human genes and misses the
canonical synthetic-lethal paralog pairs, so BOTH types are used by default (``PARALOG_TYPES``);
``gene_split`` rows (one gene split in two models) are excluded. Pair score =
max(identity, reverse identity) (% identical residues, query- and target-normalised).

Local (gene x line) features
----------------------------
``LocalFeatureBuilder.local_features(line_ids, genes)`` returns raw, unfitted arrays; any scaling is
done by the model inside a fold. Conventions:
- Missing paralog / partner (gene has none, or it is absent from the expression matrix) counts as
  expression 0 = "not expressed": a missing buffer and a silent buffer are treated alike.
- Lines with no CN / mutation data get 0 (= diploid, wild-type), the project-wide convention.
- Lines with no expression get NaN in expression-derived features (callers decide).
"""
from __future__ import annotations

import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw_external"
CACHE_DIR = ROOT / "cache"
PRIORS_DIR = CACHE_DIR / "priors"
HGNC_PATH = RAW / "hgnc" / "hgnc_complete_set.txt"
COMPARA_PATH = RAW / "ensembl_paralogs" / "Compara.116.protein_default.homologies.tsv.gz"
BIOPLEX_PATH = RAW / "bioplex_293t" / "BioPlex_293T_Network_10K_Dec_2019.tsv"
HURI_PATH = RAW / "huri" / "HuRI.tsv"
BIOGRID_PATH = RAW / "biogrid" / "BIOGRID-ALL-LATEST.tab3.zip"

PARALOG_TYPES = ("within_species_paralog", "other_paralog")
TOP_K = 3
# Deep deletion: relative CN < 0.25. Justification (scripts/build_features.py prints the evidence):
# the relative-CN distribution has a homozygous-loss mode at < 0.05 (97% of those genes have
# TPM < 1), a valley over 0.05-0.20, and the heterozygous-loss shoulder rises from 0.25 (density
# x3.3 from [0.25,0.30) to [0.30,0.35)). 0.25 is also < 1 retained copy at ploidy 4, so values
# below it cannot be single-copy retention in near-tetraploid lines.
DEEP_DEL_THR = 0.25
NOT_EXPRESSED_LOG = 1.0  # log2(TPM+1) < 1  <=>  TPM < 1
BIOGRID_DBS = ("BIOGRID",)  # human rows come from BioGRID curation; other curation DBs are model-organism


# ----------------------------------------------------------------------------- HGNC / ID mapping
def load_hgnc(path=HGNC_PATH) -> pd.DataFrame:
    cols = ["hgnc_id", "symbol", "locus_group", "location", "entrez_id", "ensembl_gene_id",
            "prev_symbol", "alias_symbol"]
    return pd.read_csv(path, sep="\t", dtype=str, usecols=cols)


def ensembl_to_entrez(hgnc: pd.DataFrame) -> dict:
    h = hgnc.dropna(subset=["ensembl_gene_id", "entrez_id"])
    h = h[~h["ensembl_gene_id"].duplicated(keep=False)]  # drop ambiguous ENSG (3 in HGNC)
    return dict(zip(h["ensembl_gene_id"], h["entrez_id"]))


def symbol_maps(hgnc: pd.DataFrame):
    """(approved symbol -> entrez, previous/alias symbol -> entrez if unambiguous)."""
    h = hgnc.dropna(subset=["entrez_id"])
    approved = dict(zip(h["symbol"], h["entrez_id"]))
    other: dict[str, set] = {}
    for col in ["prev_symbol", "alias_symbol"]:
        for syms, e in zip(h[col], h["entrez_id"]):
            if isinstance(syms, str):
                for s in syms.split("|"):
                    other.setdefault(s.strip(), set()).add(e)
    fallback = {s: next(iter(v)) for s, v in other.items() if len(v) == 1 and s not in approved}
    return approved, fallback


def map_symbols(symbols, hgnc: pd.DataFrame) -> pd.Series:
    """Symbols -> Entrez strings (approved first, then unique previous/alias); NaN if unmapped."""
    approved, fallback = symbol_maps(hgnc)
    out = [approved.get(s, fallback.get(s)) if isinstance(s, str) else None for s in symbols]
    return pd.Series(out, index=pd.Index(list(symbols)), dtype=object)


def chromosome_of(hgnc: pd.DataFrame) -> pd.Series:
    """Entrez -> chromosome ('1'..'22', 'X', 'Y'; NaN for mitochondria/unplaced) from HGNC location."""
    h = hgnc.dropna(subset=["entrez_id"])
    chrom = h["location"].fillna("").str.extract(r"^(\d{1,2}|X|Y)(?=[pq]|$| |cen)")[0]
    return pd.Series(chrom.to_numpy(), index=h["entrez_id"].to_numpy(), name="chrom")


def protein_coding_entrez(hgnc: pd.DataFrame) -> set:
    h = hgnc[(hgnc["locus_group"] == "protein-coding gene") & hgnc["entrez_id"].notna()]
    return set(h["entrez_id"])


# ----------------------------------------------------------------------------- paralogs
def load_paralog_pairs(path=COMPARA_PATH, hgnc=None, types=PARALOG_TYPES,
                       chunksize=500_000) -> pd.DataFrame:
    """Human-human paralog pairs in Entrez space, symmetric (a,b) and (b,a), one row per pair.

    Columns: a, b, identity (a->b, % of a), rev_identity (% of b), score = max of the two, type.
    Duplicate pairs (both directions listed, or several ENSG mapping to one Entrez) keep the max.
    """
    hgnc = load_hgnc() if hgnc is None else hgnc
    e2e = ensembl_to_entrez(hgnc)
    use = ["gene_stable_id", "identity", "homology_type", "homology_gene_stable_id",
           "homology_species", "homology_identity", "species"]
    parts = []
    for ch in pd.read_csv(path, sep="\t", usecols=use, chunksize=chunksize,
                          dtype={"identity": float, "homology_identity": float}):
        ch = ch[ch["homology_type"].isin(types) & (ch["species"] == "homo_sapiens")
                & (ch["homology_species"] == "homo_sapiens")]
        parts.append(ch)
    d = pd.concat(parts, ignore_index=True)
    d["a"] = d["gene_stable_id"].map(e2e)
    d["b"] = d["homology_gene_stable_id"].map(e2e)
    d = d.dropna(subset=["a", "b"])
    d = d[d["a"] != d["b"]]
    fwd = pd.DataFrame({"a": d["a"], "b": d["b"], "identity": d["identity"],
                        "rev_identity": d["homology_identity"], "type": d["homology_type"]})
    rev = pd.DataFrame({"a": d["b"], "b": d["a"], "identity": d["homology_identity"],
                        "rev_identity": d["identity"], "type": d["homology_type"]})
    p = pd.concat([fwd, rev], ignore_index=True)
    p["score"] = p[["identity", "rev_identity"]].max(axis=1)
    p = p.sort_values(["a", "b", "score"], ascending=[True, True, False])
    p = p.drop_duplicates(["a", "b"], keep="first").reset_index(drop=True)
    return p[["a", "b", "identity", "rev_identity", "score", "type"]]


def top_paralogs(pairs: pd.DataFrame, k: int = TOP_K) -> pd.DataFrame:
    """Per gene: top-k paralogs by score (ties: smaller Entrez number first), family size, max identity.

    Returns index entrez; columns par1..park (Entrez or NaN), par1_score..park_score,
    family_size (all paralogs), family_size_id20 (score >= 20%), max_identity.
    """
    p = pairs.assign(_bn=pd.to_numeric(pairs["b"]))
    p = p.sort_values(["a", "score", "_bn"], ascending=[True, False, True])
    p["rank"] = p.groupby("a").cumcount()
    g = p.groupby("a")
    out = pd.DataFrame({"family_size": g.size(),
                        "family_size_id20": g["score"].apply(lambda s: int((s >= 20).sum())),
                        "max_identity": g["score"].max()})
    top = p[p["rank"] < k]
    for r in range(k):
        t = top[top["rank"] == r].set_index("a")
        out[f"par{r + 1}"] = t["b"].reindex(out.index)
        out[f"par{r + 1}_score"] = t["score"].reindex(out.index)
    out.index.name = "entrez"
    cols = [f"par{r + 1}" for r in range(k)] + [f"par{r + 1}_score" for r in range(k)]
    return out[cols + ["family_size", "family_size_id20", "max_identity"]]


# ----------------------------------------------------------------------------- interactions
def _sym_pairs(a, b, source) -> pd.DataFrame:
    a = pd.Series(a, dtype=object).astype(str)
    b = pd.Series(b, dtype=object).astype(str)
    ok = a.notna() & b.notna() & (a != b) & (a != "nan") & (b != "nan")
    a, b = a[ok].to_numpy(), b[ok].to_numpy()
    d = pd.DataFrame({"a": np.concatenate([a, b]), "b": np.concatenate([b, a])})
    d = d.drop_duplicates().reset_index(drop=True)
    d["source"] = source
    return d


def load_bioplex(path=BIOPLEX_PATH) -> pd.DataFrame:
    d = pd.read_csv(path, sep="\t", dtype=str, usecols=["GeneA", "GeneB"])
    return _sym_pairs(d["GeneA"], d["GeneB"], "bioplex_293t")


def load_huri(path=HURI_PATH, hgnc=None) -> pd.DataFrame:
    hgnc = load_hgnc() if hgnc is None else hgnc
    e2e = ensembl_to_entrez(hgnc)
    d = pd.read_csv(path, sep="\t", header=None, names=["A", "B"], dtype=str)
    strip = lambda s: s.str.split(".").str[0]
    return _sym_pairs(strip(d["A"]).map(e2e), strip(d["B"]).map(e2e), "huri")


def load_biogrid(path=BIOGRID_PATH, chunksize=500_000, report=None) -> pd.DataFrame:
    """BioGRID physical interactions, human-human, curated by BioGRID; ORCS never read."""
    use = ["Entrez Gene Interactor A", "Entrez Gene Interactor B", "Experimental System Type",
           "Organism ID Interactor A", "Organism ID Interactor B", "Source Database"]
    z = zipfile.ZipFile(path)
    names = [n for n in z.namelist() if n.endswith(".tab3.txt")]
    if len(names) != 1 or "ORCS" in names[0].upper():
        raise ValueError(f"unexpected BioGRID archive content: {z.namelist()}")
    parts, n_rows = [], 0
    with z.open(names[0]) as f:
        for ch in pd.read_csv(f, sep="\t", usecols=use, dtype=str, chunksize=chunksize):
            n_rows += len(ch)
            ch = ch[(ch["Organism ID Interactor A"] == "9606") & (ch["Organism ID Interactor B"] == "9606")
                    & (ch["Experimental System Type"] == "physical")
                    & ch["Source Database"].isin(BIOGRID_DBS)]
            parts.append(ch[use[:2]])
    d = pd.concat(parts, ignore_index=True)
    if report is not None:
        report.update({"biogrid_file": names[0], "biogrid_rows": n_rows,
                       "biogrid_human_physical_rows": len(d)})
    return _sym_pairs(d.iloc[:, 0], d.iloc[:, 1], "biogrid_physical")


def partner_table(sources: list[pd.DataFrame]) -> pd.DataFrame:
    """Long table (a, b, source) of all per-source pairs plus source='union' (deduplicated)."""
    per = pd.concat(sources, ignore_index=True)
    union = per[["a", "b"]].drop_duplicates().assign(source="union")
    return pd.concat([per, union], ignore_index=True)


# ----------------------------------------------------------------------------- coverage
def coverage(genes, top: pd.DataFrame, partners: pd.DataFrame) -> dict:
    genes = pd.Index([str(g) for g in genes])
    cov = {"n_genes": len(genes),
           "paralog_any": float(genes.isin(top.index).mean()),
           "paralog_id20": float(genes.isin(top.index[top["family_size_id20"] > 0]).mean())}
    for src, d in partners.groupby("source"):
        cov[f"partner_{src}"] = float(genes.isin(pd.Index(d["a"].unique())).mean())
    return cov


# ----------------------------------------------------------------------------- local features
LOCAL_FEATURES = ["own_expr", "not_expressed", "expr_hinge", "own_cn", "deep_del", "own_mut",
                  "par1_expr", "par2_expr", "par3_expr", "par_min_expr", "par_max_expr",
                  "par_deep_del", "par_mut", "partner_mean_expr"]


class LocalFeatureBuilder:
    """Gene x line features from raw omics plus paralog/partner tables.

    expression: log2(TPM+1), lines x Entrez. cn_rel: raw relative CN (NOT log). damaging: 0/1/2.
    top: output of ``top_paralogs``. partners: long (a, b, source) table; ``partner_source``
    selects which source feeds partner_mean_expr ('union' by default).
    """

    def __init__(self, expression: pd.DataFrame, cn_rel: pd.DataFrame, damaging: pd.DataFrame,
                 top: pd.DataFrame, partners: pd.DataFrame, partner_source: str = "union"):
        self.E, self.CN, self.M = expression, cn_rel, damaging
        self.top = top
        self.k = sum(1 for c in top.columns if c.startswith("par") and c[3:].isdigit())
        p = partners[partners["source"] == partner_source]
        self.partners = p.groupby("a")["b"].apply(list).to_dict()

    @classmethod
    def from_cache(cls, cache_dir=CACHE_DIR, partner_source="union"):
        from essml.io import load_matrix
        pdir = Path(cache_dir) / "priors"
        return cls(load_matrix("expression", cache_dir), load_matrix("cn", cache_dir),
                   load_matrix("damaging", cache_dir),
                   pd.read_parquet(pdir / "paralog_top.parquet"),
                   pd.read_parquet(pdir / "partners.parquet"), partner_source)

    @staticmethod
    def _take(mat: pd.DataFrame, line_ids, genes, fill_line, fill_gene):
        """mat[lines, genes] as float32; absent lines -> fill_line, absent genes -> fill_gene."""
        rows = mat.index.get_indexer(line_ids)
        cols = mat.columns.get_indexer(genes)
        X = mat.to_numpy(dtype=np.float32, copy=False)
        out = np.full((len(line_ids), len(genes)), fill_gene, dtype=np.float32)
        r_ok, c_ok = rows >= 0, cols >= 0
        if r_ok.any() and c_ok.any():
            out[np.ix_(r_ok, c_ok)] = X[np.ix_(rows[r_ok], cols[c_ok])]
        out[~r_ok, :] = fill_line
        return out

    def local_features(self, line_ids, genes) -> dict[str, np.ndarray]:
        """Dict feature name -> float32 array (n_lines, n_genes). See module docstring."""
        line_ids = [str(i) for i in line_ids]
        genes = [str(g) for g in genes]
        nan = np.float32(np.nan)
        f: dict[str, np.ndarray] = {}
        E = lambda gs: self._take(self.E, line_ids, gs, nan, 0.0)
        cn = self._take(self.CN, line_ids, genes, 1.0, 1.0)
        cn = np.where(np.isnan(cn), 1.0, cn).astype(np.float32)  # NaN CN -> diploid
        f["own_expr"] = E(genes)
        f["not_expressed"] = np.where(np.isnan(f["own_expr"]), nan,
                                      (f["own_expr"] < NOT_EXPRESSED_LOG).astype(np.float32))
        f["expr_hinge"] = np.maximum(f["own_expr"] - NOT_EXPRESSED_LOG, 0)  # piecewise-linear knot at TPM=1
        f["own_cn"] = np.clip(np.log2((cn + 0.1) / 1.1), -3, 3).astype(np.float32)
        f["deep_del"] = (cn < DEEP_DEL_THR).astype(np.float32)
        f["own_mut"] = self._take(self.M, line_ids, genes, 0.0, 0.0)
        top = self.top.reindex(genes)
        par_e, par_dd, par_m = [], [], []
        for r in range(self.k):
            pg = top[f"par{r + 1}"].fillna("__none__").astype(str).tolist()
            has = np.array([g != "__none__" for g in pg])
            e = E(pg)
            par_e.append(e)
            c = self._take(self.CN, line_ids, pg, 1.0, 1.0)
            par_dd.append(((c < DEEP_DEL_THR) & has[None, :]).astype(np.float32))
            par_m.append(self._take(self.M, line_ids, pg, 0.0, 0.0) * has[None, :])
        for r in range(self.k):
            f[f"par{r + 1}_expr"] = par_e[r]
        n_par = top[[f"par{r + 1}" for r in range(self.k)]].notna().sum(axis=1).to_numpy()
        stack = np.stack(par_e)  # k x lines x genes
        # min over EXISTING paralogs only; genes with no paralog -> 0 (no buffer)
        exists = (np.arange(self.k)[:, None] < n_par[None, :])[:, None, :]
        f["par_min_expr"] = np.where(n_par[None, :] > 0,
                                     np.min(np.where(exists, stack, np.inf), axis=0), 0).astype(np.float32)
        f["par_max_expr"] = np.max(stack, axis=0).astype(np.float32)
        f["par_deep_del"] = np.max(np.stack(par_dd), axis=0)
        f["par_mut"] = np.max(np.stack(par_m), axis=0)
        f["partner_mean_expr"] = self._partner_mean(line_ids, genes)
        return f

    def _partner_mean(self, line_ids, genes) -> np.ndarray:
        from scipy import sparse
        ecols = {g: j for j, g in enumerate(self.E.columns)}
        ri, ci = [], []
        for i, g in enumerate(genes):
            js = [ecols[p] for p in self.partners.get(g, ()) if p in ecols]
            ri += [i] * len(js)
            ci += js
        if not ri:
            return np.zeros((len(line_ids), len(genes)), np.float32)
        W = sparse.csr_matrix((np.ones(len(ri), np.float32), (ri, ci)),
                              shape=(len(genes), len(ecols)))
        cnt = np.asarray(W.sum(axis=1)).ravel()
        W = sparse.diags(1.0 / np.maximum(cnt, 1)) @ W
        rows = self.E.index.get_indexer(line_ids)
        out = np.full((len(line_ids), len(genes)), np.nan, np.float32)
        ok = rows >= 0
        X = self.E.to_numpy(dtype=np.float32, copy=False)[rows[ok]]
        out[ok] = np.asarray((W @ X.T).T, dtype=np.float32)
        return out

    def static_features(self, genes) -> pd.DataFrame:
        """Line-independent gene features: family size, max identity, number of partners."""
        genes = [str(g) for g in genes]
        t = self.top.reindex(genes)
        return pd.DataFrame({
            "family_size": t["family_size"].fillna(0).to_numpy(np.float32),
            "family_size_id20": t["family_size_id20"].fillna(0).to_numpy(np.float32),
            "max_identity": t["max_identity"].fillna(0).to_numpy(np.float32),
            "n_partners": np.array([len(self.partners.get(g, ())) for g in genes], np.float32),
        }, index=pd.Index(genes, name="entrez"))


_DEFAULT_BUILDER = None


def local_features(line_ids, genes, builder: LocalFeatureBuilder | None = None):
    """Module-level convenience wrapper; loads matrices and prior tables from cache once."""
    global _DEFAULT_BUILDER
    if builder is None:
        if _DEFAULT_BUILDER is None:
            _DEFAULT_BUILDER = LocalFeatureBuilder.from_cache()
        builder = _DEFAULT_BUILDER
    return builder.local_features(line_ids, genes)
