"""External expression -> DepMap expression scale (DESIGN v3 section 7.3).

Pipeline per source
-------------------
1. **Parse** the deposited table (counts / TPM / FPKM / 3' tag counts) into a gene x sample frame keyed
   by Ensembl gene ID (unversioned) or gene symbol. Only expression tables are read; nothing here touches
   screen data or any prediction model.
2. **TPM over all genes**:
   - counts: rate = count / effective length (Ensembl 116 exon-union length - fragment length + 1),
     TPM = rate / sum(rate) * 1e6 over every gene that has a length;
   - FPKM: TPM = FPKM / sum(FPKM) * 1e6 (all rows of the table);
   - TPM (RSEM): spike-ins (non-ENSG rows) removed, renormalised to 1e6;
   - 3' tag counts (Quant-seq): CPM, *no* length correction; flagged ``tpm_compatible = False``;
   - sources whose table is pre-filtered ("restricted" universe) cannot give TPM over all genes: their
     protein-coding subset is rescaled to the DepMap dev median protein-coding TPM sum (~0.94M).
3. **Subset** to DepMap protein-coding Entrez columns (Ensembl -> Entrez via HGNC; symbol -> Entrez via
   HGNC approved, then unique previous/alias symbols); several Ensembl genes mapping to one Entrez are
   summed; then **log2(TPM + 1)**.
4. **Harmonise** with one of the candidates: ``none``, ``quantile`` (to the mean sorted DepMap dev
   profile), ``quantile_location`` (quantile + gene-wise offset fitted on matched DepMap dev lines),
   ``rank`` (within-profile ranks; both sides). Celligner-style alignment is deferred (DESIGN v3).
5. **Evaluate / report**: same-line nearest-neighbour recall among DepMap dev lines on the top-5000
   variable genes, profile agreement, and an applicability-domain distance (mean 1 - r to the 10
   nearest dev lines).

Missing genes stay NaN (the feature code fills absent expression with the training mean), except for
restricted-universe sources where an absent gene was filtered for low expression and is set to 0.
"""
from __future__ import annotations

import gzip
import io
import re
import tarfile
import zipfile
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
GTF_PATH = ROOT / "data" / "raw_external" / "ensembl116_gtf" / "Homo_sapiens.GRCh38.116.gtf.gz"
HGNC_PATH = ROOT / "data" / "raw_external" / "hgnc" / "hgnc_complete_set.txt"
OUT_DIR = ROOT / "cache" / "external_expr"

_ENSG = re.compile(r"^(ENSG\d{11})(?:\.\d+)?(?:_PAR_Y)?$")


# ============================================================================ annotation
def merged_length(starts, ends) -> int:
    """Total length of the union of closed intervals [start, end] (1-based, inclusive)."""
    s = np.asarray(starts, dtype=np.int64)
    e = np.asarray(ends, dtype=np.int64)
    if s.size == 0:
        return 0
    o = np.argsort(s, kind="stable")
    s, e = s[o], e[o]
    total, cur_s, cur_e = 0, s[0], e[0]
    for a, b in zip(s[1:], e[1:]):
        if a > cur_e + 1:          # disjoint (adjacent intervals are merged, no double count anyway)
            total += cur_e - cur_s + 1
            cur_s, cur_e = a, b
        else:
            cur_e = max(cur_e, b)
    return int(total + cur_e - cur_s + 1)


def _attr(s: str, key: str):
    m = re.search(key + r' "([^"]*)"', s)
    return m.group(1) if m else None


def gene_table_from_gtf(path=GTF_PATH) -> pd.DataFrame:
    """Ensembl GTF -> DataFrame indexed by unversioned gene_id with gene_name, gene_biotype, chrom and
    exon-union length (all transcripts' exons merged)."""
    exons: dict[str, tuple[list, list]] = {}
    meta: dict[str, tuple] = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if len(f) < 9:
                continue
            if f[2] == "gene":
                gid = _attr(f[8], "gene_id")
                meta[gid] = (_attr(f[8], "gene_name"), _attr(f[8], "gene_biotype"), f[0])
            elif f[2] == "exon":
                gid = _attr(f[8], "gene_id")
                st, en = exons.setdefault(gid, ([], []))
                st.append(int(f[3]))
                en.append(int(f[4]))
    rows = []
    for gid, (name, biotype, chrom) in meta.items():
        st, en = exons.get(gid, ([], []))
        rows.append((gid, name, biotype, chrom, merged_length(st, en)))
    t = pd.DataFrame(rows, columns=["gene_id", "gene_name", "gene_biotype", "chrom", "length"])
    return t.set_index("gene_id")


def load_gene_table(cache_dir=OUT_DIR, gtf=GTF_PATH) -> pd.DataFrame:
    p = Path(cache_dir) / "ensembl116_genes.parquet"
    if p.exists():
        return pd.read_parquet(p)
    t = gene_table_from_gtf(gtf)
    p.parent.mkdir(parents=True, exist_ok=True)
    t.to_parquet(p)
    return t


# ============================================================================ ID mapping
def strip_version(ids) -> pd.Index:
    """'ENSG00000000003.14' -> 'ENSG00000000003'; non-Ensembl IDs unchanged."""
    out = []
    for x in ids:
        m = _ENSG.match(str(x))
        out.append(m.group(1) if m else str(x))
    return pd.Index(out)


@dataclass
class GeneMapper:
    """Ensembl / symbol -> Entrez (HGNC) and symbol -> Ensembl (Ensembl GTF names, then HGNC)."""
    hgnc: pd.DataFrame
    genes: pd.DataFrame                       # gene_table_from_gtf output
    depmap_symbols: dict = field(default_factory=dict)   # entrez -> DepMap column symbol

    def __post_init__(self):
        h = self.hgnc
        he = h.dropna(subset=["ensembl_gene_id", "entrez_id"])
        he = he[~he["ensembl_gene_id"].duplicated(keep=False)]
        self.ens2entrez = dict(zip(he["ensembl_gene_id"], he["entrez_id"]))
        hs = h.dropna(subset=["entrez_id"])
        self.sym_approved = dict(zip(hs["symbol"], hs["entrez_id"]))
        other: dict[str, set] = {}
        for col in ["prev_symbol", "alias_symbol"]:
            for syms, e in zip(hs[col], hs["entrez_id"]):
                if isinstance(syms, str):
                    for s in syms.split("|"):
                        other.setdefault(s.strip(), set()).add(e)
        self.sym_fallback = {s: next(iter(v)) for s, v in other.items()
                             if len(v) == 1 and s not in self.sym_approved}
        g = self.genes.dropna(subset=["gene_name"])
        uniq = g[~g["gene_name"].duplicated(keep=False)]
        self.sym2ens_gtf = dict(zip(uniq["gene_name"], uniq.index))
        hens = h.dropna(subset=["ensembl_gene_id"])
        self.sym2ens_hgnc = dict(zip(hens["symbol"], hens["ensembl_gene_id"]))
        self.entrez2ens = {}
        for ens, ent in self.ens2entrez.items():
            self.entrez2ens.setdefault(ent, ens)

    def symbol_to_entrez(self, s):
        return self.sym_approved.get(s, self.sym_fallback.get(s))

    def symbol_to_ensembl(self, s):
        """GTF unique gene_name, else HGNC approved symbol's Ensembl, else via Entrez of a previous/alias."""
        e = self.sym2ens_gtf.get(s) or self.sym2ens_hgnc.get(s)
        if e is None:
            ent = self.symbol_to_entrez(s)
            e = self.entrez2ens.get(ent) if ent else None
        return e

    def ensembl_projection(self, depmap_cols) -> dict:
        """Ensembl gene -> DepMap Entrez column (HGNC; DepMap symbol via GTF name as fallback)."""
        cols = set(map(str, depmap_cols))
        proj = {ens: ent for ens, ent in self.ens2entrez.items() if ent in cols}
        have = set(proj.values())
        for ent in cols - have:
            sym = self.depmap_symbols.get(ent)
            ens = self.sym2ens_gtf.get(sym) if sym else None
            if ens and ens not in proj:
                proj[ens] = ent
        return proj

    def symbol_projection(self, symbols, depmap_cols) -> dict:
        cols = set(map(str, depmap_cols))
        out = {}
        for s in symbols:
            ent = self.symbol_to_entrez(s)
            if ent in cols:
                out[s] = ent
        sym2ent_dm = {v: k for k, v in self.depmap_symbols.items()}
        for s in symbols:
            if s not in out and s in sym2ent_dm:
                out[s] = sym2ent_dm[s]
        return out


def collapse_rows(df: pd.DataFrame, keys) -> pd.DataFrame:
    """Sum rows sharing a key (e.g. Ensembl versions, duplicated symbols); NaN keys dropped."""
    keys = pd.Index(keys)
    ok = ~pd.isna(keys)
    d = df.loc[ok].copy()
    d.index = keys[ok]
    return d.groupby(level=0, sort=False).sum(min_count=1)


# ============================================================================ unit conversion
def counts_to_tpm(counts: pd.DataFrame, lengths: pd.Series, fragment_length: float = 0.0):
    """Counts (genes x samples) -> TPM over all genes with a positive effective length.

    effective length = length - fragment_length + 1 if fragment_length > 0 else length. Genes with no
    length or effective length <= 0 get no TPM (as RSEM). Returns (tpm, fraction of counts used).
    """
    L = lengths.reindex(counts.index).astype(float)
    eff = L - fragment_length + 1.0 if fragment_length > 0 else L
    use = eff.notna() & (eff > 0)
    c = counts.loc[use].astype(float)
    rate = c.div(eff[use], axis=0)
    tpm = rate / rate.sum(axis=0) * 1e6
    frac_used = c.sum(axis=0) / counts.astype(float).sum(axis=0)
    return tpm, frac_used


def scale_to_million(x: pd.DataFrame) -> pd.DataFrame:
    """FPKM -> TPM, CPM from counts, or re-normalising a TPM table after dropping rows."""
    x = x.astype(float)
    return x / x.sum(axis=0) * 1e6


def project_to_depmap(tpm: pd.DataFrame, proj: dict, depmap_cols, absent: str = "nan",
                      pc_target: float | None = None) -> pd.DataFrame:
    """TPM (genes x samples; index = Ensembl or symbol) -> log2(TPM+1), samples x DepMap Entrez columns.

    proj: row key -> Entrez column. Rows sharing an Entrez are summed. Columns without any mapped row:
    NaN (absent='nan') or 0 TPM (absent='zero', pre-filtered tables). pc_target: if given, the mapped
    protein-coding subset is rescaled to sum to pc_target per sample (restricted-universe sources).
    """
    cols = pd.Index(list(map(str, depmap_cols)))
    keys = [proj.get(k) for k in tpm.index]
    sub = collapse_rows(tpm, keys)
    sub = sub.reindex(cols)
    if absent == "zero":
        sub = sub.fillna(0.0)
    if pc_target is not None:
        sub = sub * (pc_target / sub.sum(axis=0, min_count=1))
    return np.log2(sub.T.astype(float) + 1.0)


# ============================================================================ parsers
@dataclass
class RawTable:
    values: pd.DataFrame        # genes x samples (raw unit)
    id_type: str                # 'ensembl' | 'symbol'
    unit: str                   # 'counts' | 'tpm' | 'fpkm' | 'tags'
    extra: dict = field(default_factory=dict)


def _read_htseq(fh) -> pd.Series:
    s = pd.read_csv(fh, sep="\t", header=None, index_col=0, compression="gzip").iloc[:, 0]
    return s[~s.index.astype(str).str.startswith("__")]


def parse_htseq_tar(path, samples) -> RawTable:
    """GEO RAW.tar of per-GSM HTSeq count files; ``samples`` = GSM IDs to keep (file prefix)."""
    cols = {}
    with tarfile.open(path) as tf:
        for m in tf.getmembers():
            gsm = m.name.split("_")[0]
            if gsm in samples:
                cols[gsm] = _read_htseq(io.BytesIO(tf.extractfile(m).read()))
    missing = set(samples) - set(cols)
    if missing:
        raise KeyError(f"samples not in {path}: {sorted(missing)}")
    v = pd.DataFrame(cols)[list(samples)]
    return RawTable(collapse_rows(v, strip_version(v.index)), "ensembl", "counts")


def parse_htseq_files(directory, samples, unit="tags") -> RawTable:
    cols = {}
    for p in sorted(Path(directory).glob("GSM*.txt.gz")):
        gsm = p.name.split("_")[0]
        if gsm in samples:
            with open(p, "rb") as fh:
                cols[gsm] = _read_htseq(fh)
    v = pd.DataFrame(cols)[list(samples)]
    return RawTable(collapse_rows(v, strip_version(v.index)), "ensembl", unit)


def parse_csv_counts(path, samples) -> RawTable:
    v = pd.read_csv(path, index_col=0)[list(samples)]
    return RawTable(collapse_rows(v, strip_version(v.index)), "ensembl", "counts")


def parse_tsv_counts(directory, files, samples) -> RawTable:
    parts = []
    for f in files:
        t = pd.read_csv(Path(directory) / f, sep="\t", index_col=0)
        keep = [c for c in t.columns if c in samples]
        if keep:
            parts.append(t[keep].groupby(level=0, sort=False).sum(min_count=1))  # duplicated symbols summed
    v = pd.concat(parts, axis=1)   # genes absent from one file's (pre-filtered) table -> NaN -> 0 below
    return RawTable(v[list(samples)].fillna(0.0), "symbol", "counts")


def parse_calderon(path, samples) -> RawTable:
    v = pd.read_csv(path, sep="\t")      # header has one field fewer than rows: first column = index
    v = v[list(samples)]
    return RawTable(collapse_rows(v, strip_version(v.index)), "ensembl", "counts")


def parse_sack_tar(path, samples) -> RawTable:
    cols = {}
    with tarfile.open(path) as tf:
        for m in tf.getmembers():
            name = re.sub(r"^GSM\d+_", "", m.name).replace(".txt.gz", "")
            if name in samples:
                t = pd.read_csv(io.BytesIO(tf.extractfile(m).read()), sep="\t", index_col=0,
                                compression="gzip")
                cols[name] = t.iloc[:, 0]
    v = pd.DataFrame(cols)[list(samples)].fillna(0.0)
    return RawTable(collapse_rows(v, strip_version(v.index)), "ensembl", "counts")


def parse_drager_xlsx(path, sheet, samples) -> RawTable:
    import openpyxl
    wb = openpyxl.load_workbook(path, read_only=True)
    ws = wb[sheet]                                  # only the expression sheet is read
    it = ws.iter_rows(values_only=True)
    header = next(it)
    df = pd.DataFrame(list(it), columns=header)
    wb.close()
    df = df.dropna(subset=["ensembl_id"]).set_index("ensembl_id")
    v = df[list(samples)].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    return RawTable(collapse_rows(v, strip_version(v.index)), "ensembl", "tpm")


def parse_tzelepis_fpkm(path, member, samples) -> RawTable:
    outer_name, inner_name = member.split("/", 1)
    with zipfile.ZipFile(path) as z:
        inner = io.BytesIO(z.read(outer_name))
    with zipfile.ZipFile(inner) as z2:
        t = pd.read_csv(io.BytesIO(z2.read(inner_name)), sep="\t", index_col=0)
    t = t.groupby(level=0, sort=False).sum(min_count=1)
    return RawTable(t[list(samples)], "symbol", "fpkm")


def parse_encode_rsem(directory, samples) -> RawTable:
    """ENCODE RSEM gene quantifications: TPM (spike-ins dropped later) + expected_count."""
    tpm, cnt = {}, {}
    for s in samples:
        p = Path(directory) / f"{s}_gene_quantifications_GRCh38.tsv.gz"
        t = pd.read_csv(p, sep="\t", index_col=0, usecols=["gene_id", "expected_count", "TPM"])
        tpm[s], cnt[s] = t["TPM"], t["expected_count"]
    v = pd.DataFrame(tpm)
    c = pd.DataFrame(cnt)
    return RawTable(collapse_rows(v, strip_version(v.index)), "ensembl", "tpm",
                    extra={"expected_count": collapse_rows(c, strip_version(c.index))})


def parse_source(spec: dict, samples, root=ROOT) -> RawTable:
    p = Path(root) / spec["path"]
    kind = spec["parser"]
    if kind == "encode_rsem":
        return parse_encode_rsem(p, samples)
    if kind == "tzelepis_fpkm":
        return parse_tzelepis_fpkm(p, spec["member"], samples)
    if kind == "htseq_tar":
        return parse_htseq_tar(p, samples)
    if kind == "htseq_files":
        return parse_htseq_files(p, samples)
    if kind == "csv_counts":
        return parse_csv_counts(p, samples)
    if kind == "tsv_counts":
        return parse_tsv_counts(p, spec["files"], samples)
    if kind == "calderon":
        return parse_calderon(p, samples)
    if kind == "sack_tar":
        return parse_sack_tar(p, samples)
    if kind == "drager_xlsx":
        return parse_drager_xlsx(p, spec["sheet"], samples)
    raise ValueError(kind)


def to_tpm(raw: RawTable, genes: pd.DataFrame, mapper: GeneMapper, fragment_length: float = 0.0):
    """RawTable -> (TPM over all genes, genes x samples, keyed by Ensembl where possible; info)."""
    v = raw.values
    info = {}
    if raw.id_type == "symbol" and raw.unit == "counts":
        keys = [mapper.symbol_to_ensembl(s) for s in v.index]
        info["frac_rows_mapped_to_ensembl"] = float(np.mean([k is not None for k in keys]))
        v = collapse_rows(v, keys)
    if raw.unit == "counts":
        tpm, used = counts_to_tpm(v, genes["length"], fragment_length)
        info["frac_counts_with_length"] = used.round(4).to_dict()
    elif raw.unit == "tpm":
        is_gene = v.index.str.match(r"^ENSG\d{11}$") if raw.id_type == "ensembl" else np.ones(len(v), bool)
        info["frac_tpm_spikein_or_nongene"] = (1 - v.loc[is_gene].sum() / v.sum()).round(4).to_dict()
        tpm = scale_to_million(v.loc[is_gene])
    elif raw.unit in ("fpkm", "tags"):
        tpm = scale_to_million(v)
    else:
        raise ValueError(raw.unit)
    return tpm, info


# ============================================================================ harmonisation
def reference_quantiles(ref: np.ndarray) -> np.ndarray:
    """Mean of the row-wise sorted reference matrix (lines x genes, no NaN)."""
    return np.sort(np.asarray(ref, dtype=float), axis=1).mean(axis=0)


def _avg_rank(x: np.ndarray) -> np.ndarray:
    from scipy.stats import rankdata
    return rankdata(x, method="average")


def quantile_normalise(x: pd.DataFrame, ref_q: np.ndarray) -> pd.DataFrame:
    """Map each row's observed values onto the reference distribution (ties get the mean rank's value).

    A row with n observed genes uses plotting positions (rank - 0.5) / n, interpolated on the
    reference quantiles at (i + 0.5) / N; with n = N and no ties the row becomes exactly ref_q
    in its own order.
    """
    ref_q = np.asarray(ref_q, dtype=float)
    N = ref_q.size
    grid = (np.arange(N) + 0.5) / N
    out = np.full(x.shape, np.nan)
    X = x.to_numpy(dtype=float)
    for i in range(X.shape[0]):
        ok = ~np.isnan(X[i])
        n = ok.sum()
        if n == 0:
            continue
        p = (_avg_rank(X[i, ok]) - 0.5) / n
        out[i, ok] = np.interp(p, grid, ref_q)
    return pd.DataFrame(out, index=x.index, columns=x.columns)


def rank_transform(x: pd.DataFrame) -> pd.DataFrame:
    """Within-row average rank / n of observed genes (in (0, 1])."""
    X = x.to_numpy(dtype=float)
    out = np.full(X.shape, np.nan)
    for i in range(X.shape[0]):
        ok = ~np.isnan(X[i])
        if ok.any():
            out[i, ok] = _avg_rank(X[i, ok]) / ok.sum()
    return pd.DataFrame(out, index=x.index, columns=x.columns)


def fit_location(ext: pd.DataFrame, dep: pd.DataFrame) -> pd.Series:
    """Gene-wise mean of (DepMap - external) over matched rows (row i of ext <-> row i of dep)."""
    d = dep.to_numpy(dtype=float) - ext.to_numpy(dtype=float)
    n = (~np.isnan(d)).sum(axis=0)
    off = np.where(n > 0, np.nansum(d, axis=0) / np.maximum(n, 1), 0.0)
    return pd.Series(off, index=ext.columns)


def apply_location(ext: pd.DataFrame, offsets: pd.Series) -> pd.DataFrame:
    return ext + offsets.reindex(ext.columns).fillna(0.0).to_numpy()[None, :]


def heldout_location(ext: pd.DataFrame, dep: pd.DataFrame, groups, fit_rows=None,
                     audit: list | None = None) -> pd.DataFrame:
    """Each row corrected with offsets fitted only on rows of OTHER groups (leave-one-group-out).

    ext, dep: aligned (row i of ext <-> row i of dep). groups: label per row (matched DepMap line for
    LOO-line, source for LOSO). fit_rows: boolean mask of rows eligible for fitting (e.g. exact dev
    matches; default all). ``audit`` (optional list) receives (left-out group, training row labels).
    """
    groups = np.asarray(list(groups), dtype=object)
    fit_rows = np.ones(len(ext), bool) if fit_rows is None else np.asarray(fit_rows, bool)
    out = ext.copy()
    for g in pd.unique(groups):
        test = groups == g
        train = fit_rows & ~test
        if audit is not None:
            audit.append((g, list(ext.index[train])))
        off = fit_location(ext.loc[train], dep.loc[train]) if train.any() else pd.Series(0.0, index=ext.columns)
        out.loc[test] = apply_location(ext.loc[test], off).to_numpy()
    return out


# ============================================================================ metrics
def top_variable_genes(ref: pd.DataFrame, n: int = 5000) -> pd.Index:
    v = ref.var(axis=0)
    return v.sort_values(ascending=False, kind="stable").index[:n]


def corr_to_pool(x: pd.DataFrame, pool: pd.DataFrame) -> pd.DataFrame:
    """Pearson r of each row of x with each row of pool over genes observed in that x row."""
    P = pool.to_numpy(dtype=float)
    out = np.full((x.shape[0], P.shape[0]), np.nan)
    X = x.to_numpy(dtype=float)
    for i in range(X.shape[0]):
        ok = ~np.isnan(X[i]) & ~np.isnan(P).any(axis=0)
        if ok.sum() < 3:
            continue
        a = X[i, ok] - X[i, ok].mean()
        B = P[:, ok] - P[:, ok].mean(axis=1, keepdims=True)
        out[i] = (B @ a) / (np.sqrt((B ** 2).sum(axis=1)) * np.sqrt((a ** 2).sum()) + 1e-300)
    return pd.DataFrame(out, index=x.index, columns=pool.index)


def nn_recall(r: pd.DataFrame, truth: pd.Series) -> pd.DataFrame:
    """Rank of the true line among pool columns (1 = nearest) and recall@1/5/10 per row."""
    rows = []
    for pid, t in truth.items():
        rr = r.loc[pid]
        rank = int((rr > rr[t]).sum()) + 1
        rows.append(dict(profile_id=pid, true_rank=rank, r_true=float(rr[t]), nearest=rr.idxmax(),
                         r_nearest=float(rr.max()), recall1=rank <= 1, recall5=rank <= 5, recall10=rank <= 10))
    return pd.DataFrame(rows).set_index("profile_id")


def pearson_rows(a: pd.DataFrame, b: pd.DataFrame) -> pd.Series:
    """Row-wise Pearson between aligned rows of a and b over jointly observed genes."""
    out = {}
    for (ia, ra), (_, rb) in zip(a.iterrows(), b.iterrows()):
        ok = ra.notna() & rb.notna()
        out[ia] = float(np.corrcoef(ra[ok], rb[ok])[0, 1]) if ok.sum() > 2 else np.nan
    return pd.Series(out)


def applicability(r: pd.DataFrame, lineage: pd.Series, k: int = 10) -> pd.DataFrame:
    """Mean correlation distance (1 - r) to the k nearest pool lines, nearest lineage (majority among k;
    ties broken by mean distance) and nearest line."""
    rows = []
    for pid, rr in r.iterrows():
        top = rr.sort_values(ascending=False).iloc[:k]
        d = 1.0 - top
        lin = lineage.reindex(top.index)
        cnt = lin.value_counts()
        best = cnt[cnt == cnt.max()].index
        if len(best) > 1:
            best = [d[lin == b].mean() for b in best], list(best)
            best = best[1][int(np.argmin(best[0]))]
        else:
            best = best[0]
        rows.append(dict(profile_id=pid, ad_distance=float(d.mean()), nearest_lineage=best,
                         nearest_lineage_n=int(cnt.max()), nearest_line=top.index[0],
                         nearest_line_lineage=lineage.get(top.index[0]), r_nearest=float(top.iloc[0])))
    return pd.DataFrame(rows).set_index("profile_id")


def pool_self_ad(pool_r: pd.DataFrame, k: int = 10) -> pd.Series:
    """AD distance of each pool line to its k nearest OTHER pool lines (reference distribution)."""
    R = pool_r.to_numpy(dtype=float).copy()
    np.fill_diagonal(R, -np.inf)
    top = -np.sort(-R, axis=1)[:, :k]
    return pd.Series((1 - top).mean(axis=1), index=pool_r.index)
