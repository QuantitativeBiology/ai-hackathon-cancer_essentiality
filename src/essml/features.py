"""Line-level feature matrices (DESIGN.md §2 "Inputs").

Two layers, so that nothing fitted across lines can leak:

1. ``build_*`` functions: RAW, UNFITTED, row-wise transforms. Each value depends only on that line's
   own data plus fixed external knowledge (HGNC, MSigDB), never on statistics over other lines.
   They are built for every model with data (dev, test, set A, unscreened lines) and cached by
   ``scripts/build_features.py`` under ``cache/features/``. External cells can be pushed through
   the same functions.
2. ``FoldTransformer``: every cross-line statistic (SD filter, min-count filter, standardisation,
   imputation means, one-hot levels, gene-set z-scores) is fitted on TRAINING IDs only and then
   applied to any IDs.

Convention: 0 = diploid / wild-type
-----------------------------------
Genomic raw features are encoded so that 0 means "diploid, wild-type, no event":
C = log2((CN + 0.1) / 1.1) is exactly 0 at relative CN 1; M/H/F indicators are 0 when no event is
called; composite flags are 0 when intact. Missing genomic values (NaN CN of a gene in a line,
or a line with no CN / mutation / fusion data) are filled with 0, i.e. assumed diploid wild-type,
and the FoldTransformer adds one line-level ``<block>:absent`` flag for lines without that
modality. For external cells, engineered genotypes (e.g. TP53-/-) are set explicitly; everything
else follows the convention. Expression is never zero-filled (absent lines get the training mean).

Raw blocks (name -> content)
----------------------------
E           log2(TPM+1), DepMap protein-coding genes (as supplied)
C           clip(log2((CN+0.1)/1.1), -3, 3), NaN -> 0, genes NaN in every line dropped
C_instab    fraction of autosomal genes with |relative CN - 1| > 0.3 (raw relative CN)
M           damaging-mutation dosage 0/1/2 (all genes; >= 5 training lines filter at fit time)
M_burden    log1p(number of genes with a damaging mutation)
H           hotspot indicator (value > 0)
F_pair      fusion gene-pair indicators (unordered Entrez pair), after filtering
F_gene      gene-level indicator "gene is a partner in >= 1 retained fusion"
meta        OncotreeLineage, OncotreePrimaryDisease, Sex (categorical), Age (numeric)
meta_prov   GrowthPattern, ModelType, PrimaryOrMetastasis (ABLATION ONLY: proxy the screen batch)
composite   p53 / Rb pathway flags (see ``build_composites``)
P           Nusinow 2020 TMT proteomics (normalised log2 relative abundance; NaN kept)
Gene-set scores (p53-target score, Hallmarks) are fitted blocks on E: kind='score'.
"""
from __future__ import annotations

import json
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from essml import priors

ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data"
CACHE_DIR = ROOT / "cache"
FEATURE_DIR = CACHE_DIR / "features"
SPLITS = ROOT / "configs" / "splits.csv"
HALLMARK_PATH = ROOT / "data" / "raw_external" / "msigdb_hallmark" / "h.all.v2024.1.Hs.symbols.gmt"
PROTEOMICS_PATH = ROOT / "data" / "external_omics" / "nusinow2020" / "protein_quant_current_normalized.csv.gz"

CN_OFFSET = 0.1
CN_CLIP = 3.0
INSTAB_THR = 0.3          # |relative CN - 1| > 0.3 (R2 data review)
AMP_THR = 2.0             # relative CN > 2 (>= 2x ploidy) = amplification
LOSS_THR = 0.75           # relative CN < 0.75 = copy loss (1 of 2 or 1 of 3-4 copies lost)
DEEP_DEL_THR = priors.DEEP_DEL_THR
P53_TARGETS = ["CDKN1A", "MDM2", "RRM2B", "ZMAT3", "TNFRSF10B", "FDXR", "SESN1", "DDB2"]
COMPOSITE_GENES = ["TP53", "MDM2", "MDM4", "PPM1D", "ATM", "CHEK2", "RB1", "CDKN2A"]
META_COLS = ["OncotreeLineage", "OncotreePrimaryDisease", "Sex", "Age"]
META_PROV_COLS = ["GrowthPattern", "ModelType", "PrimaryOrMetastasis"]
FORBIDDEN_FIT_ROLES = ("test", "noncancerous")


# ============================================================================ raw transforms
def cn_transform(cn: pd.DataFrame) -> pd.DataFrame:
    """Relative CN -> clip(log2((CN+0.1)/1.1), -3, 3); all-NaN genes dropped; other NaN -> 0."""
    keep = cn.notna().any(axis=0)
    x = cn.loc[:, keep].to_numpy(dtype=np.float32)
    y = np.clip(np.log2((x + CN_OFFSET) / (1.0 + CN_OFFSET)), -CN_CLIP, CN_CLIP)
    y = np.nan_to_num(y, nan=0.0).astype(np.float32)
    return pd.DataFrame(y, index=cn.index, columns=cn.columns[keep.to_numpy()])


def cn_instability(cn: pd.DataFrame, chrom: pd.Series, thr: float = INSTAB_THR) -> pd.Series:
    """Fraction of autosomal genes (non-NaN) with |relative CN - 1| > thr, per line."""
    auto = {str(i) for i in range(1, 23)}
    cols = [c for c in cn.columns if chrom.get(c) in auto]
    x = cn[cols].to_numpy(dtype=np.float32)
    ok = ~np.isnan(x)
    dev = (np.abs(np.where(ok, x, 1.0) - 1.0) > thr) & ok
    return pd.Series((dev.sum(1) / np.maximum(ok.sum(1), 1)).astype(np.float32), index=cn.index,
                     name="C_instab")


def mutation_burden(damaging: pd.DataFrame) -> pd.Series:
    n = (damaging.to_numpy() > 0).sum(axis=1)
    return pd.Series(np.log1p(n).astype(np.float32), index=damaging.index, name="M_burden")


def hotspot_binary(hotspot: pd.DataFrame) -> pd.DataFrame:
    return (hotspot > 0).astype(np.float32)


def _parse_fusion_gene(s: pd.Series):
    """'RARG (ENSG00000172819.17)' -> (symbol, ENSG without version)."""
    m = s.astype(str).str.extract(r"^(.*?) \((ENSG\d+)(?:\.\d+)?(?:_PAR_Y)?\)$")
    return m[0], m[1]


def filter_fusions(fus: pd.DataFrame, coding_entrez: set, ens2entrez: dict,
                   sym2entrez: dict | None = None) -> tuple[pd.DataFrame, dict]:
    """Drop self-fusions and fusions with a non-coding partner.

    A partner is protein-coding if its Entrez (via Ensembl ID, else approved symbol) is in
    ``coding_entrez`` (expression-matrix genes or HGNC locus_group 'protein-coding gene').
    Unmapped partners (RP11-*, lncRNA, pseudogenes without Entrez) count as non-coding.
    Self-fusion: same Ensembl ID, same symbol, or same Entrez after mapping.
    Returns long table (ModelID, a, b) with a < b (unordered pair) and a report dict.
    """
    sym2entrez = sym2entrez or {}
    s1, e1 = _parse_fusion_gene(fus["Gene1"])
    s2, e2 = _parse_fusion_gene(fus["Gene2"])
    z1 = e1.map(ens2entrez).fillna(s1.map(sym2entrez))
    z2 = e2.map(ens2entrez).fillna(s2.map(sym2entrez))
    self_f = (e1 == e2) | (s1 == s2) | (z1.notna() & (z1 == z2))
    coding = z1.isin(coding_entrez) & z2.isin(coding_entrez)
    keep = ~self_f & coding
    a = np.where(pd.to_numeric(z1, errors="coerce") <= pd.to_numeric(z2, errors="coerce"), z1, z2)
    b = np.where(pd.to_numeric(z1, errors="coerce") <= pd.to_numeric(z2, errors="coerce"), z2, z1)
    out = pd.DataFrame({"ModelID": fus["ModelID"].to_numpy(), "a": a, "b": b})[keep.to_numpy()]
    out = out.drop_duplicates().reset_index(drop=True)
    rep = {"rows_in": int(len(fus)), "self_fusion_rows": int(self_f.sum()),
           "noncoding_rows": int((~self_f & ~coding).sum()), "rows_kept": int(keep.sum()),
           "unique_pairs_kept": int(out[["a", "b"]].drop_duplicates().shape[0])}
    return out, rep


def fusion_matrices(long: pd.DataFrame, index, min_lines_total: int = 5):
    """(pair matrix, gene matrix) of 0/1 over ``index``.

    Columns present in < ``min_lines_total`` lines overall are not cached: they can never pass
    the fit-time filter (>= 5 TRAINING lines), so this pre-filter is a superset and uses no labels.
    """
    index = pd.Index(index, name="ModelID")
    long = long[long["ModelID"].isin(index)]
    pair = long.assign(col=long["a"] + "--" + long["b"])[["ModelID", "col"]].drop_duplicates()
    gene = pd.concat([long[["ModelID", "a"]].rename(columns={"a": "col"}),
                      long[["ModelID", "b"]].rename(columns={"b": "col"})]).drop_duplicates()

    def to_mat(t):
        n = t.groupby("col")["ModelID"].nunique()
        cols = sorted(n.index[n >= min_lines_total])
        t = t[t["col"].isin(cols)]
        m = np.zeros((len(index), len(cols)), np.float32)
        m[index.get_indexer(t["ModelID"]), pd.Index(cols).get_indexer(t["col"])] = 1.0
        return pd.DataFrame(m, index=index, columns=pd.Index(cols, name="feature"))

    return to_mat(pair), to_mat(gene)


def build_meta(model: pd.DataFrame, cols=META_COLS) -> pd.DataFrame:
    m = model.set_index("ModelID")[list(cols)].copy()
    for c in cols:
        if c == "Age":
            m[c] = pd.to_numeric(m[c], errors="coerce").astype(np.float32)
        else:
            m[c] = m[c].astype("string").str.strip().replace({"": pd.NA, "Unknown": pd.NA})
    return m


def build_composites(cn_rel: pd.DataFrame, damaging: pd.DataFrame, hotspot: pd.DataFrame,
                     sym2entrez: dict, index=None) -> pd.DataFrame:
    """p53 / Rb pathway flags from raw relative CN, damaging and hotspot matrices (row-wise).

    0 = intact (wild-type / diploid); lines lacking a modality get 0 for flags from it.
    Columns: TP53_dam (0/1/2), TP53_hot, TP53_mut (dam or hot), TP53_cnloss (CN < 0.75),
    MDM2_amp, MDM4_amp (CN > 2), PPM1D_dam, ATM_dam, CHEK2_dam (binary), p53_altered (any of the
    binary p53 flags), RB1_dam, RB1_deepdel, CDKN2A_deepdel (CN < 0.25), RB_altered (any).
    The p53-target expression score is a fitted 'score' block (needs cross-line z-scores).
    """
    g = {s: sym2entrez[s] for s in COMPOSITE_GENES}
    idx = pd.Index(index if index is not None else
                   cn_rel.index.union(damaging.index).union(hotspot.index), name="ModelID")

    def col(mat, sym, fill=0.0):
        e = g[sym]
        if e not in mat.columns:
            return pd.Series(fill, index=idx, dtype=np.float32)
        return mat[e].reindex(idx).astype(np.float32).fillna(fill)

    cn = lambda s: col(cn_rel, s, fill=1.0)
    out = pd.DataFrame(index=idx)
    out["TP53_dam"] = col(damaging, "TP53")
    out["TP53_hot"] = (col(hotspot, "TP53") > 0).astype(np.float32)
    out["TP53_mut"] = ((out["TP53_dam"] > 0) | (out["TP53_hot"] > 0)).astype(np.float32)
    out["TP53_cnloss"] = (cn("TP53") < LOSS_THR).astype(np.float32)
    out["MDM2_amp"] = (cn("MDM2") > AMP_THR).astype(np.float32)
    out["MDM4_amp"] = (cn("MDM4") > AMP_THR).astype(np.float32)
    for s in ["PPM1D", "ATM", "CHEK2"]:
        out[f"{s}_dam"] = (col(damaging, s) > 0).astype(np.float32)
    out["p53_altered"] = out[["TP53_mut", "TP53_cnloss", "MDM2_amp", "MDM4_amp", "PPM1D_dam",
                              "ATM_dam", "CHEK2_dam"]].max(axis=1)
    out["RB1_dam"] = (col(damaging, "RB1") > 0).astype(np.float32)
    out["RB1_deepdel"] = (cn("RB1") < DEEP_DEL_THR).astype(np.float32)
    out["CDKN2A_deepdel"] = (cn("CDKN2A") < DEEP_DEL_THR).astype(np.float32)
    out["RB_altered"] = out[["RB1_dam", "RB1_deepdel"]].max(axis=1)
    return out.astype(np.float32)


def read_gmt(path, sym2entrez: dict) -> dict[str, list[str]]:
    sets = {}
    for line in Path(path).read_text().splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) > 2:
            ents = [sym2entrez.get(s) for s in parts[2:]]
            sets[parts[0]] = sorted({e for e in ents if e is not None}, key=int)
    return sets


def _sym_lookup(hgnc) -> dict:
    approved, fallback = priors.symbol_maps(hgnc)
    return {**fallback, **approved}


def hallmark_sets(hgnc=None, path=HALLMARK_PATH) -> dict[str, list[str]]:
    hgnc = priors.load_hgnc() if hgnc is None else hgnc
    return read_gmt(path, _sym_lookup(hgnc))


def p53_target_set(hgnc=None) -> dict[str, list[str]]:
    hgnc = priors.load_hgnc() if hgnc is None else hgnc
    lk = _sym_lookup(hgnc)
    return {"p53_targets": [lk[s] for s in P53_TARGETS]}


def build_proteomics(path=PROTEOMICS_PATH, model: pd.DataFrame | None = None, hgnc=None,
                     report: dict | None = None) -> pd.DataFrame:
    """Nusinow et al. 2020 CCLE TMT proteomics -> lines (ModelID) x Entrez, float32, NaN kept.

    Sample columns '<CCLEName>_TenPxNN' are matched to Model.csv CCLEName; lines measured in
    two plexes are averaged (NaN-aware). Proteins map to Entrez via HGNC (approved symbol, then
    unique previous/alias symbol); several proteins per Entrez -> the one quantified in most
    lines (ties: Swiss-Prot before TrEMBL, then file order).
    """
    model = pd.read_csv(DATA_DIR / "Model.csv", usecols=["ModelID", "CCLEName"]) if model is None else model
    hgnc = priors.load_hgnc() if hgnc is None else hgnc
    p = pd.read_csv(path)
    samp = [c for c in p.columns if "_TenPx" in c and not c.endswith("_Peptides")]
    names = pd.Index([c.rsplit("_TenPx", 1)[0] for c in samp])
    cc = model.dropna(subset=["CCLEName"])
    if cc["CCLEName"].duplicated().any():
        raise ValueError("duplicated CCLEName in Model.csv")
    name2id = dict(zip(cc["CCLEName"], cc["ModelID"]))
    ids = names.map(lambda n: name2id.get(n))
    X = p[samp].to_numpy(np.float32)
    ent = priors.map_symbols(p["Gene_Symbol"].tolist(), hgnc).to_numpy()
    nobs = (~np.isnan(X)).sum(1)
    order = pd.DataFrame({"ent": ent, "nobs": -nobs, "tr": p["Protein_Id"].str.startswith("tr").to_numpy(),
                          "pos": np.arange(len(p))}).dropna(subset=["ent"])
    order = order.sort_values(["ent", "nobs", "tr", "pos"]).drop_duplicates("ent")
    order = order.sort_values("pos")
    Xg = X[order["pos"].to_numpy()].T  # samples x proteins
    df = pd.DataFrame(Xg, index=pd.Index(ids, name="ModelID"),
                      columns=pd.Index(order["ent"].astype(str).to_numpy(), name="entrez"))
    matched = df.index.notna()
    df = df[matched].groupby(level=0).mean().astype(np.float32)
    if report is not None:
        report.update({"samples": len(samp), "unique_ccle_names": int(names.nunique()),
                       "lines_matched": int(df.shape[0]), "unmatched_names": sorted(set(names[~matched])),
                       "replicated_lines": int(names.duplicated().sum()), "proteins": int(len(p)),
                       "proteins_mapped": int(pd.notna(ent).sum()), "entrez_genes": int(df.shape[1])})
    return df


# ============================================================================ build all (raw)
def build_all(cache_dir=CACHE_DIR, report: dict | None = None) -> dict[str, pd.DataFrame]:
    """Build every raw block from the cached DepMap matrices (generator-friendly, one at a time)."""
    from essml.io import load_matrix
    rep = report if report is not None else {}
    hgnc = priors.load_hgnc()
    lk = _sym_lookup(hgnc)
    out = {}
    E = load_matrix("expression", cache_dir)
    out["E"] = E
    cn = load_matrix("cn", cache_dir)
    out["C"] = cn_transform(cn)
    out["C_instab"] = cn_instability(cn, priors.chromosome_of(hgnc)).to_frame()
    dam = load_matrix("damaging", cache_dir)
    hot = load_matrix("hotspot", cache_dir)
    out["composite"] = build_composites(cn, dam, hot, lk)
    del cn
    out["M"] = dam
    out["M_burden"] = mutation_burden(dam).to_frame()
    out["H"] = hotspot_binary(hot)
    fus = load_matrix("fusion", cache_dir)
    coding = set(E.columns) | priors.protein_coding_entrez(hgnc)
    long, frep = filter_fusions(fus, coding, priors.ensembl_to_entrez(hgnc), lk)
    f_index = E.index.union(pd.Index(fus["ModelID"].unique()))
    out["F_pair"], out["F_gene"] = fusion_matrices(long, f_index)
    rep["fusion_filter"] = frep
    model = pd.read_csv(DATA_DIR / "Model.csv", dtype=str)
    out["meta"] = build_meta(model, META_COLS)
    out["meta_prov"] = build_meta(model, META_PROV_COLS)
    prep = {}
    out["P"] = build_proteomics(model=model[["ModelID", "CCLEName"]], hgnc=hgnc, report=prep)
    rep["proteomics"] = prep
    return out


def save_raw(blocks: dict[str, pd.DataFrame], feature_dir=FEATURE_DIR):
    feature_dir = Path(feature_dir)
    feature_dir.mkdir(parents=True, exist_ok=True)
    for k, v in blocks.items():
        v = v.copy()
        v.columns = v.columns.astype(str)
        v.index.name = "ModelID"
        v.to_parquet(feature_dir / f"{k}.parquet")


def load_raw(names, feature_dir=FEATURE_DIR) -> dict[str, pd.DataFrame]:
    return {n: pd.read_parquet(Path(feature_dir) / f"{n}.parquet") for n in names}


def load_gene_sets(feature_dir=FEATURE_DIR) -> dict[str, dict[str, list[str]]]:
    return json.loads((Path(feature_dir) / "gene_sets.json").read_text())


# ============================================================================ fitted transforms
def _f64(df: pd.DataFrame) -> np.ndarray:
    """C-contiguous float64 copy: reductions then do not depend on pandas' block layout."""
    return np.ascontiguousarray(df.to_numpy(dtype=np.float64))


@dataclass
class BlockSpec:
    """How a raw block is turned into model features inside a fold.

    kind: 'continuous' | 'count' | 'categorical' | 'score'
    min_sd: drop columns with training SD < min_sd (columns with SD == 0 are always dropped).
    min_count: ('count') keep columns with value != 0 in >= min_count training lines.
    max_missing: drop columns whose NaN fraction among present training lines > max_missing;
        remaining NaNs get the training column mean.
    absent_fill: lines not in the raw block: 'zero' (raw 0 = diploid WT) or 'mean' (training mean).
    absent_flag: add one '<block>:absent' 0/1 column.
    standardize: z-score with training mean/SD (categorical: numeric columns only).
    source / gene_sets / min_set_genes: ('score') mean z-score over each set's genes in raw block
        ``source``; z uses training mean/SD per gene; sets with < min_set_genes usable genes dropped.
    """
    kind: str = "continuous"
    min_sd: float = 0.0
    min_count: int = 0
    max_missing: float = 1.0
    absent_fill: str = "zero"
    absent_flag: bool = False
    standardize: bool = True
    min_level_count: int = 1
    source: str | None = None
    gene_sets: dict | None = None
    min_set_genes: int = 3


def default_specs(gene_sets: dict | None = None, include_provenance: bool = False) -> dict[str, BlockSpec]:
    """DESIGN §2 defaults. gene_sets: {'hallmark': {...}, 'p53_targets': {...}} for score blocks.

    meta_prov (GrowthPattern/ModelType/PrimaryOrMetastasis) is ablation-only and is included only
    with include_provenance=True. Drop keys for blocks a configuration does not use.
    """
    s = {
        "E": BlockSpec(min_sd=0.1, absent_fill="mean"),
        "C": BlockSpec(absent_flag=True),
        "C_instab": BlockSpec(),
        "M": BlockSpec(kind="count", min_count=5, absent_flag=True),
        "M_burden": BlockSpec(absent_fill="mean"),
        "H": BlockSpec(kind="count", min_count=5, absent_flag=True),
        "F_pair": BlockSpec(kind="count", min_count=5, absent_flag=True),
        "F_gene": BlockSpec(kind="count", min_count=5, absent_flag=True),
        "meta": BlockSpec(kind="categorical"),
        "composite": BlockSpec(),
        "P": BlockSpec(max_missing=0.2, absent_fill="mean", absent_flag=True),
    }
    if include_provenance:
        s["meta_prov"] = BlockSpec(kind="categorical")
    if gene_sets:
        if "p53_targets" in gene_sets:
            s["p53score"] = BlockSpec(kind="score", source="E", gene_sets=gene_sets["p53_targets"])
        if "hallmark" in gene_sets:
            s["hallmark"] = BlockSpec(kind="score", source="E", gene_sets=gene_sets["hallmark"],
                                      min_set_genes=10)
    return s


def check_fit_ids(train_ids, splits_path=SPLITS, forbidden=FORBIDDEN_FIT_ROLES):
    """Leakage assertion (DESIGN §8): fit IDs must not be test / set-A lines."""
    ids = pd.Index([str(i) for i in train_ids])
    if ids.duplicated().any():
        raise ValueError("duplicated training IDs")
    p = Path(splits_path)
    if p.exists():
        role = pd.read_csv(p, index_col=0)["role"].reindex(ids)
        bad = role[role.isin(forbidden)]
        if len(bad):
            raise ValueError(f"fit called with {len(bad)} forbidden IDs (roles {sorted(set(bad))})")
    return ids


class FoldTransformer:
    """Fit cross-line statistics on training IDs only; transform any IDs.

    >>> ft = FoldTransformer({"E": spec_E, "C": spec_C}).fit(raw, train_ids)
    >>> X_val = ft.transform(raw, val_ids)          # DataFrame, columns '<block>:<feature>'

    ``raw`` maps block name -> raw DataFrame (rows ModelID). Only rows in ``train_ids`` are read
    during ``fit`` (rows are selected before any statistic is computed). ``params_`` holds every
    fitted quantity; ``fitted_params()`` flattens it for comparison in tests.
    """

    def __init__(self, specs: dict[str, BlockSpec], check_roles: bool = True, splits_path=SPLITS):
        self.specs = specs
        self.check_roles = check_roles
        self.splits_path = splits_path
        self.params_: dict[str, dict] = {}

    # -- helpers
    @staticmethod
    def _rows(df: pd.DataFrame, ids):
        present = df.index.intersection(ids)
        return df.loc[present], present

    def fit(self, raw: dict[str, pd.DataFrame], train_ids):
        ids = (check_fit_ids(train_ids, self.splits_path) if self.check_roles
               else pd.Index([str(i) for i in train_ids]))
        self.train_ids_ = ids
        self.params_ = {}
        for name, spec in self.specs.items():
            src = raw[spec.source if spec.kind == "score" else name]
            X, present = self._rows(src, ids)
            fit = {"score": self._fit_score, "categorical": self._fit_cat}.get(spec.kind, self._fit_num)
            self.params_[name] = fit(X, spec)
            self.params_[name]["n_train_present"] = len(present)
        return self

    def _fit_num(self, X: pd.DataFrame, spec: BlockSpec) -> dict:
        A = _f64(X)
        nanfrac = np.isnan(A).mean(0) if len(A) else np.zeros(A.shape[1])
        keep = nanfrac <= spec.max_missing
        with np.errstate(all="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)  # all-NaN column -> mean 0 below
            mean = np.nanmean(A, axis=0) if len(A) else np.zeros(A.shape[1])
        mean = np.nan_to_num(mean)
        Af = np.where(np.isnan(A), mean, A)
        sd = Af.std(axis=0) if len(A) else np.zeros(A.shape[1])
        if spec.kind == "count":
            keep &= (Af != 0).sum(axis=0) >= spec.min_count
        keep &= (sd > 0) & (sd >= spec.min_sd)
        return {"columns": X.columns[keep].tolist(), "mean": mean[keep], "sd": sd[keep]}

    def _fit_cat(self, X: pd.DataFrame, spec: BlockSpec) -> dict:
        p = {"levels": {}, "numeric": {}}
        for c in X.columns:
            if pd.api.types.is_numeric_dtype(X[c]):
                v = X[c].astype(float)
                med = float(v.median()) if v.notna().any() else 0.0
                f = v.fillna(med)
                p["numeric"][c] = {"median": med, "mean": float(f.mean()) if len(f) else 0.0,
                                   "sd": float(f.std(ddof=0)) if len(f) > 1 else 0.0}
            else:
                vc = X[c].dropna().astype(str).value_counts()
                p["levels"][c] = sorted(vc.index[vc >= spec.min_level_count])
        return p

    def _fit_score(self, X: pd.DataFrame, spec: BlockSpec) -> dict:
        A = _f64(X)
        mean = A.mean(0) if len(A) else np.zeros(A.shape[1])
        sd = A.std(0) if len(A) else np.zeros(A.shape[1])
        col = {c: j for j, c in enumerate(X.columns)}
        sets = {}
        for sname, genes in spec.gene_sets.items():
            js = [col[g] for g in genes if g in col and sd[col[g]] > 0]
            if len(js) >= spec.min_set_genes:
                sets[sname] = {"genes": [X.columns[j] for j in js], "mean": mean[js], "sd": sd[js]}
        p = {"sets": sets}
        if spec.standardize and sets:  # second-level standardisation of the scores (training)
            S = self._scores(X, sets)
            p["score_mean"], p["score_sd"] = S.mean(0), np.where(S.std(0) > 0, S.std(0), 1.0)
        return p

    @staticmethod
    def _scores(X: pd.DataFrame, sets: dict) -> np.ndarray:
        out = np.zeros((len(X), len(sets)))
        for k, s in enumerate(sets.values()):
            Z = (_f64(X[s["genes"]]) - s["mean"]) / s["sd"]
            out[:, k] = Z.mean(1)
        return out

    # -- transform
    def transform(self, raw: dict[str, pd.DataFrame], ids) -> pd.DataFrame:
        ids = pd.Index([str(i) for i in ids], name="ModelID")
        blocks = []
        for name, spec in self.specs.items():
            p = self.params_[name]
            src = raw[spec.source if spec.kind == "score" else name]
            absent = ~ids.isin(src.index)
            if spec.kind == "score":
                B = self._tr_score(src, ids, absent, p, spec)
            elif spec.kind == "categorical":
                B = self._tr_cat(src, ids, p, spec)
            else:
                B = self._tr_num(src, ids, absent, p, spec)
            B.columns = [f"{name}:{c}" for c in B.columns]
            if spec.absent_flag:
                B[f"{name}:absent"] = absent.astype(np.float32)
            blocks.append(B)
        return pd.concat(blocks, axis=1).astype(np.float32)

    def _tr_num(self, src, ids, absent, p, spec):
        cols = p["columns"]
        A = _f64(src.reindex(index=ids, columns=cols))
        fill_absent = np.zeros(len(cols)) if spec.absent_fill == "zero" else p["mean"]
        A[absent] = fill_absent
        A = np.where(np.isnan(A), p["mean"], A)
        if spec.standardize:
            A = (A - p["mean"]) / p["sd"]
        return pd.DataFrame(A.astype(np.float32), index=ids, columns=cols)

    def _tr_cat(self, src, ids, p, spec):
        X = src.reindex(ids)
        parts = {}
        for c, levels in p["levels"].items():
            v = X[c].astype("string")
            for lv in levels:
                parts[f"{c}={lv}"] = (v == lv).fillna(False).to_numpy(np.float32)
        for c, s in p["numeric"].items():
            v = X[c].astype(float)
            parts[f"{c}_missing"] = v.isna().to_numpy(np.float32)
            f = v.fillna(s["median"]).to_numpy()
            if spec.standardize and s["sd"] > 0:
                f = (f - s["mean"]) / s["sd"]
            parts[c] = f.astype(np.float32)
        return pd.DataFrame(parts, index=ids)

    def _tr_score(self, src, ids, absent, p, spec):
        sets = p["sets"]
        genes = sorted({g for s in sets.values() for g in s["genes"]})
        X = src.reindex(index=ids, columns=genes)
        S = self._scores(X, sets)
        if "score_mean" in p:
            S = (S - p["score_mean"]) / p["score_sd"]
        S[absent] = 0.0  # absent lines -> training average
        S = np.nan_to_num(S)
        return pd.DataFrame(S.astype(np.float32), index=ids, columns=list(sets))

    def fit_transform(self, raw, train_ids):
        return self.fit(raw, train_ids).transform(raw, train_ids)

    def fitted_params(self) -> dict[str, np.ndarray]:
        flat = {}

        def walk(prefix, obj):
            if isinstance(obj, dict):
                for k, v in obj.items():
                    walk(f"{prefix}/{k}", v)
            elif isinstance(obj, (list, tuple)):
                flat[prefix] = np.asarray(obj, dtype=object)
            else:
                flat[prefix] = np.asarray(obj)
        walk("", self.params_)
        return flat
