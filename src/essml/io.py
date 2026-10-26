"""Fast-load cache for DepMap matrices.

Column-naming convention
------------------------
Raw DepMap gene columns look like ``"SYMBOL (EntrezID)"``. Cached matrices use
the **Entrez ID as a string** (e.g. ``"7105"``) as column name; the symbol is
kept only in ``cache/gene_map.parquet`` (entrez, symbol, in_matrices).
Rows are indexed by ``ModelID`` (``ACH-xxxxxx``); omics matrices use only the
row with ``IsDefaultEntryForModel == "Yes"`` for each model. Values are float32.

Matrix names: crispr, expression, cn, damaging, hotspot (model x gene matrices)
and fusion (long table, default-entry rows only, not a matrix).
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data"
CACHE_DIR = ROOT / "cache"

FILES = {
    "crispr": "CRISPRGeneEffect.csv",
    "expression": "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
    "cn": "OmicsCNGeneWGS.csv",
    "damaging": "OmicsSomaticMutationsMatrixDamaging.csv",
    "hotspot": "OmicsSomaticMutationsMatrixHotspot.csv",
    "fusion": "OmicsFusionFiltered.csv",
}
MATRICES = ["crispr", "expression", "cn", "damaging", "hotspot"]
META_COLS = ["SequencingID", "ModelConditionID", "ModelID",
             "IsDefaultEntryForMC", "IsDefaultEntryForModel"]
_COL_RE = re.compile(r"^(.*) \((\d+)\)$")


def parse_gene_column(col: str):
    """'TSPAN6 (7105)' -> ('TSPAN6', '7105'); None if not parseable."""
    m = _COL_RE.match(str(col))
    return (m.group(1), m.group(2)) if m else None


def is_default(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin(["yes", "true"])


def select_default_entries(df: pd.DataFrame) -> pd.DataFrame:
    """Keep rows with IsDefaultEntryForModel == 'Yes', indexed by ModelID.

    Raises if a model has more than one default row.
    """
    d = df[is_default(df["IsDefaultEntryForModel"])]
    if d["ModelID"].duplicated().any():
        raise ValueError("more than one default entry for a ModelID")
    return d.set_index("ModelID")


def read_header(path) -> list[str]:
    return list(pd.read_csv(path, nrows=0).columns)


def split_columns(columns):
    """Return (meta_cols, parsed [(raw_col, symbol, entrez)], unparseable_cols)."""
    meta, parsed, bad = [], [], []
    for i, c in enumerate(columns):
        if i == 0 or c in META_COLS or str(c).startswith("Unnamed"):
            meta.append(c)
            continue
        p = parse_gene_column(c)
        if p is None:
            bad.append(c)
        else:
            parsed.append((c, p[0], p[1]))
    return meta, parsed, bad


def convert_matrix(path, name: str, chunksize: int = 200):
    """Raw matrix CSV -> (float32 DataFrame [ModelID x Entrez str], info dict).

    CRISPR: first column is ModelID. Omics: default entries only.
    Duplicated Entrez columns: first kept, count in info['n_dup_entrez'].
    """
    header = read_header(path)
    _, parsed, bad = split_columns(header)
    gene_cols = [p[0] for p in parsed]
    dtypes = {c: np.float32 for c in gene_cols}
    is_crispr = name == "crispr"
    if is_crispr:
        usecols = [header[0]] + gene_cols
    else:
        usecols = [c for c in header if c in META_COLS] + gene_cols
    parts = []
    for ch in pd.read_csv(path, usecols=usecols, dtype=dtypes, chunksize=chunksize):
        if is_crispr:
            ch = ch.rename(columns={header[0]: "ModelID"}).set_index("ModelID")[gene_cols]
        else:
            ch = select_default_entries(ch)[gene_cols]
        parts.append(ch.astype(np.float32))
    mat = pd.concat(parts)
    ent = np.array([p[2] for p in parsed])
    keep = ~pd.Index(ent).duplicated()
    mat = mat.iloc[:, np.flatnonzero(keep)]
    mat.columns = pd.Index(ent[keep], name="entrez")
    mat.index.name = "ModelID"
    sym = {p[2]: p[1] for p in parsed[::-1]}  # first occurrence wins
    info = {"n_dup_entrez": int((~keep).sum()), "unparseable": bad,
            "symbols": {e: sym[e] for e in mat.columns}}
    return mat, info


def build_cache(data_dir=DATA_DIR, cache_dir=CACHE_DIR, names=None, force=False):
    """Convert raw CSVs to parquet under cache_dir and (re)write gene_map.parquet."""
    data_dir, cache_dir = Path(data_dir), Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    for name in (names or MATRICES):
        out = cache_dir / f"{name}.parquet"
        if out.exists() and not force:
            continue
        mat, info = convert_matrix(data_dir / FILES[name], name)
        mat.to_parquet(out)
        pd.DataFrame({"entrez": list(info["symbols"]), "symbol": list(info["symbols"].values())}
                     ).to_parquet(cache_dir / f"_symbols_{name}.parquet", index=False)
        del mat
    fp = cache_dir / "fusion.parquet"
    if (names is None or "fusion" in names) and (force or not fp.exists()):
        f = pd.read_csv(data_dir / FILES["fusion"], index_col=0)
        f[is_default(f["IsDefaultEntryForModel"])].reset_index(drop=True).to_parquet(fp)
    write_gene_map(cache_dir)
    return cache_dir


def write_gene_map(cache_dir=CACHE_DIR):
    """gene_map.parquet: entrez, symbol, in_matrices (comma-separated names)."""
    cache_dir = Path(cache_dir)
    sym, present = {}, {}
    for name in MATRICES:
        p = cache_dir / f"_symbols_{name}.parquet"
        if p.exists():
            t = pd.read_parquet(p)
            for e, s in zip(t["entrez"], t["symbol"]):
                sym.setdefault(e, s)
                present.setdefault(e, []).append(name)
    gm = pd.DataFrame({"entrez": list(present), "symbol": [sym[e] for e in present],
                       "in_matrices": [",".join(v) for v in present.values()]})
    gm.to_parquet(cache_dir / "gene_map.parquet", index=False)
    return gm


def load_matrix(name: str, cache_dir=CACHE_DIR, unguarded: bool = False) -> pd.DataFrame:
    """Load cached matrix: float32, index=ModelID (default entries), columns=Entrez ID strings.

    name in {'crispr','expression','cn','damaging','hotspot'}; 'fusion' returns the long
    default-entry fusion table. Builds the cache on first use.

    CRISPR targets are guarded: use ``essml.guard.load_targets(ids)``. ``unguarded=True`` is
    reserved for the pre-split data audit and split construction (scripts/audit_data.py,
    scripts/make_splits.py), which ran before the test split existed.
    """
    if name == "crispr" and not unguarded:
        raise PermissionError("CRISPR targets must be loaded via essml.guard.load_targets(ids)")
    p = Path(cache_dir) / f"{name}.parquet"
    if not p.exists():
        build_cache(cache_dir=cache_dir, names=[name])
    return pd.read_parquet(p)


def load_gene_map(cache_dir=CACHE_DIR) -> pd.DataFrame:
    return pd.read_parquet(Path(cache_dir) / "gene_map.parquet")
