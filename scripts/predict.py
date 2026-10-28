"""Predict CRISPR gene effects (Chronos scale) for new cell lines from their expression profile, with a
per-gene confidence score (DESIGN §4 "Per-gene predictability", D22).

    python scripts/predict.py --expression my_lines.csv --out predictions.csv
    python scripts/predict.py --expression ext.tsv --harmonise quantile --out preds.parquet   # external profiles
    python scripts/predict.py --expression x.csv --save-models          # persist the fitted model for next time
    python scripts/predict.py --expression x.csv --run-id <cv run id>   # use a run's persisted fold models
    python scripts/predict.py --expression x.csv --cn cn.csv --damaging dam.csv --hotspot hot.csv \
        --metadata meta.csv --out preds.csv                              # m5_avg with user-supplied genomics

Input (``--expression``): a csv / tsv / parquet table of log2(TPM + 1) on DepMap's protein-coding genes, one
row per sample (first column = sample id). Gene columns are Entrez IDs ("7105") or DepMap headers
("TSPAN6 (7105)"); ``--genes-in-rows`` reads the transposed layout. Genes absent from the input are set to
the training mean by the model's feature code. External (non-DepMap) profiles should be harmonised with
``--harmonise quantile`` (the frozen D9 choice: essml.expr_harmonise.quantile_normalise to the mean sorted
DepMap dev profile). TPM must be computed over all genes before subsetting (DESIGN §7.3; see
essml.expr_harmonise for count / FPKM conversion). DepMap release files (OmicsExpressionTPMLogp1...,
OmicsCNGeneWGS, OmicsSomaticMutationsMatrixDamaging / Hotspot) can be passed as they are: their identifier columns
(SequencingID, ModelID, ...) are recognised, rows are restricted to IsDefaultEntryForModel and indexed by ModelID.

Genomic inputs (optional; m5_avg only; added post hoc on 2026-10-28, REPORT §7 item 4). Same table layout as
``--expression``, one row per sample, sample IDs matching the expression table:
- ``--cn``: gene-level relative copy number, linear scale (DepMap OmicsCNGeneWGS convention, 1 = the line's
  ploidy). Transformed exactly as for DepMap lines (essml.features.cn_transform).
- ``--damaging``: damaging-mutation dosage 0 / 1 / 2 (DepMap OmicsSomaticMutationsMatrixDamaging: LikelyLoF;
  2 = summed allele fraction > 0.95).
- ``--hotspot``: hotspot-mutation matrix 0 / 1 / 2 (DepMap OmicsSomaticMutationsMatrixHotspot: OncoKB / COSMIC
  tier 1 hotspots); any value > 0 counts.
- ``--metadata``: OncotreeLineage, OncotreePrimaryDisease, Sex, Age (Model.csv column names and vocabulary;
  levels unseen in training match no level).
A modality that is not supplied, or a sample missing from its table, is absent = diploid wild type / unknown
metadata (DESIGN §2), as without these options. A gene missing from a supplied CN / mutation table is diploid /
wild type. Fusions are not read: none of m5_avg's three inputs uses them. The calls must follow DepMap's pipelines
(Mutect2 + LikelyLoF / hotspot annotation, WGS relative CN); other pipelines are not harmonised.

Models (``--model``):
- ``m5_avg`` (default; the G4 headline, D23f): the plain average of m1_E_resB0, M2 (per-gene paralog stage,
  nested local lambda) and M3 (LightGBM on M1 residuals, nested early stopping), each refitted on the dev
  lines exactly as in CV (essml.models.m5_stack.StackM5.predict_from_blocks). Only expression is read from the
  input unless the genomic options above are given. Without them M2/M3's genomic inputs (copy number,
  damaging / hotspot mutations) are absent = diploid wild type and metadata (lineage, sex, age) unknown: the
  DESIGN §2 convention used for every external cell (and for model T, §5). A warning says so: mutation-driven
  dependencies (D24) then rest on expression alone. With all four genomic inputs in DepMap's conventions this is
  the evaluated headline (dev-CV P1 0.3175, test 0.319): DepMap release rows of fold-0 dev lines fed through this
  script reproduce the registered dev-CV predictions (scripts/predict_genomics_check.py ->
  results/predict_genomics_check.json).
  Fitting takes several minutes (three models); ``--save-models`` keeps the fitted inputs for next time.
  Tiers: configs/predictable_genes.csv (primary).
  Evaluation of THIS use (expression only): post hoc, dev CV only (G7, after the REPORT red team;
  scripts/predict_expression_only_cv.py -> results/predict_expression_only_cv.json, docs/predict_expression_only.md).
  The G6 fold models applied to their held-out dev lines with expression only give P1 (variable genes) 0.3135,
  against 0.3081 for m1_E_resB0 (paired fold difference +0.0054, SE 0.0001) and 0.3175 for m5_avg with its real
  genomic inputs (the G4 test number 0.319 is with real genomics). The default rule, fixed by the PI before the
  number was seen (DESIGN §4: > delta 0.005 and > 1 SE over m1_E_resB0), is met, narrowly, so m5_avg stays the
  default. NRAS / KRAS / BRAF / PIK3CA hotspot AUPRC (dev CV) is 0.22 / 0.77 / 0.55 / 0.60 expression-only, about
  m1_E_resB0's 0.21 / 0.76 / 0.56 / 0.60 and below m5_avg with genomics 0.63 / 0.81 / 0.66 / 0.65.
- ``m1_E_resB0`` (expression-only kernel ridge, B0 residual; the evaluated expression-only G4 finalist, test P1
  0.310; about 1.5 min). Tiers: configs/predictable_genes_m1.csv (secondary). Genomic options are refused.
Either model is fitted on the DepMap dev lines (minus ``--exclude-lines``) with the launcher's sentinel, or
loaded from cache/models/<tag>/ if it was persisted earlier with ``--save-models``; ``--run-id`` instead
averages a CV run's persisted fold models (scripts/run.py --save-models) via essml.modelstore.predict_new.

Output (``--out``, csv or parquet): one row per (sample, gene): sample, entrez, symbol, predicted_effect,
confidence, tier, ad_distance, ad_percentile_vs_dev, nearest_lineage, note. ``note`` flags NRAS / KRAS / BRAF /
PIK3CA for both models (RAS/RAF/PIK3CA hotspot dependencies, underpredicted by both expression-only configurations;
D24), with a one-line warning; for m5_avg the flag is dropped when ``--hotspot`` is supplied. ``confidence`` / ``tier`` come from
the dev-CV per-gene predictability table of the model used (m1_E_resB0: configs/predictable_genes_m1.csv,
the secondary table; the frozen G4 tiers in configs/predictable_genes.csv are those of m5_avg). The
applicability-domain distance is the mean (1 - Pearson r) to the 10 nearest DepMap dev lines on the top-5000
variable genes (essml.expr_harmonise.applicability), with its percentile among dev lines. A probability of
effect < -0.5 is not produced: the models are point regressors and no calibrated probability exists in the code.
"""
from __future__ import annotations

import os
import sys

if __name__ == "__main__":
    _n = sys.argv[sys.argv.index("--threads") + 1] if "--threads" in sys.argv else "4"
    for _v in ["OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"]:
        os.environ[_v] = _n

import argparse  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import re  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
from essml import expr_harmonise as XH  # noqa: E402
from essml import guard, modelstore  # noqa: E402
from essml import predictability as PR  # noqa: E402

SUPPORTED = {"m5_avg": {"config": ROOT / "configs" / "runs" / "m5_avg.json",
                        "table": ROOT / "configs" / "predictable_genes.csv"},
             "m1_E_resB0": {"config": ROOT / "configs" / "runs" / "m1_E_resB0.json",
                            "table": ROOT / "configs" / "predictable_genes_m1.csv"}}
DEFAULT_MODEL = "m5_avg"
GENOMIC_INPUT_MODELS = {"m5_avg"}  # read genomics when available; here only expression is supplied (DESIGN §2)
# Post-hoc dev-CV evaluation of m5_avg applied with expression only (results/predict_expression_only_cv.json; G7).
EXPRESSION_ONLY_CV = {"m5_avg_P1_variable": 0.3135, "m1_E_resB0_P1_variable": 0.3081, "dP1_vs_m1": 0.0054,
                      "dP1_vs_m1_se": 0.0001, "m5_avg_with_genomics_P1_variable": 0.3175,
                      "nras_auprc_m5_avg_expression_only": 0.22, "nras_auprc_m1_E_resB0": 0.21,
                      "nras_auprc_m5_avg_with_genomics": 0.63}
NO_GENOMICS_NOTE = ("mutation-driven dependency; no mutation / copy-number input (diploid wild type assumed); "
                    "underpredicted by both expression-only configurations")
NO_GENOMICS_WARNING = ("WARNING: m5_avg reads copy number, mutations, hotspots, fusions and metadata when it has them; this "
                       "script supplies expression only, so they are set to diploid wild type / unknown (DESIGN section 2 "
                       "convention). This expression-only use was evaluated only post hoc, in dev CV "
                       "(results/predict_expression_only_cv.json): P1 on variable genes 0.3135, vs 0.3081 for "
                       "m1_E_resB0 and 0.3175 for m5_avg with its real genomic inputs (the G4 test number 0.319 used real "
                       "genomics). Mutation-driven dependencies (D24: NRAS, KRAS, BRAF, PIK3CA hotspot addiction) are "
                       "underpredicted by both expression-only configurations: NRAS AUPRC 0.22 here vs 0.21 for "
                       "m1_E_resB0 and 0.63 for m5_avg with genomics (see the 'note' column).")
# PI note (D24, docs/mutation_dependencies.md): expression-only prediction misses RAS/RAF/PI3K hotspot addiction within
# lineages. With genomics m5_avg recovers most of it (NRAS AUPRC 0.63), but this script supplies no genomics, so both
# models here are expression-only (NRAS AUPRC 0.21 m1_E_resB0, 0.22 m5_avg). Reporting only; tiers are unaffected.
MUTATION_NOTE_GENES = {"NRAS", "KRAS", "BRAF", "PIK3CA"}
MUTATION_NOTE = "mutation-driven dependency; underpredicted by both expression-only configurations"
M1_WARNING = ("WARNING: m1_E_resB0 is expression-only. RAS/RAF/PIK3CA hotspot dependencies (NRAS, KRAS, BRAF, PIK3CA; "
              "D24) are underpredicted by both expression-only configurations (this model and m5_avg as run by this "
              "script, which supplies no genomics); see the 'note' column.")
EXPRESSION_ONLY_MODELS = {"m1_E_resB0"}
AD_K = 10
N_HVG = 5000
X_BUILDER = (lambda E: {"E": E})  # raw expression block for the kernel models
_HDR = re.compile(r"^.*\((\d+)\)\s*$")
DEPMAP_ID_COLS = {"SequencingID", "ModelConditionID", "ModelID", "ProfileID", "IsDefaultEntryForMC",
                  "IsDefaultEntryForModel"}
GENOMIC_OPTIONS = ("cn", "damaging", "hotspot", "metadata")


# ============================================================================ input
def _read_table(path, genes_in_rows=False) -> pd.DataFrame:
    p = Path(path)
    if p.suffix == ".parquet":
        df = pd.read_parquet(p)
    else:
        df = pd.read_csv(p, sep="\t" if p.suffix in (".tsv", ".txt") else ",", index_col=0, low_memory=False)
    if genes_in_rows:
        df = df.T
    if "IsDefaultEntryForModel" in df.columns:  # DepMap release layout -> default entries, indexed by ModelID
        df = df[df["IsDefaultEntryForModel"].astype(str).str.lower().isin(["true", "yes", "1"])]
        df = df.set_index("ModelID")
    return df.drop(columns=[c for c in df.columns if c in DEPMAP_ID_COLS])


def read_expression(path, genes_in_rows=False) -> pd.DataFrame:
    """Samples x Entrez-string columns, float64."""
    df = _read_table(path, genes_in_rows)
    cols = []
    for c in map(str, df.columns):
        m = _HDR.match(c)
        cols.append(m.group(1) if m else c.strip())
    df.columns = cols
    if not all(c.isdigit() for c in cols):
        bad = [c for c in cols if not c.isdigit()][:5]
        raise ValueError(f"gene columns must be Entrez IDs or 'SYMBOL (ENTREZ)'; got e.g. {bad}")
    df = df.T.groupby(level=0, sort=False).mean().T  # duplicated Entrez columns -> mean
    df.index = pd.Index([str(i) for i in df.index], name="sample")
    return df.astype(np.float64)


def read_genomic(path, genes_in_rows=False, log=print) -> pd.DataFrame:
    """CN / mutation table -> samples x Entrez-string columns, parsed as the DepMap cache does (essml.io
    convert_matrix): columns without an Entrez ID (e.g. 'ENSG..._PAR_Y') are dropped, duplicated Entrez -> first."""
    df = _read_table(path, genes_in_rows)
    ent = [(m.group(1) if (m := _HDR.match(str(c))) else str(c).strip()) for c in df.columns]
    ok = np.array([e.isdigit() for e in ent])
    if (~ok).any():
        log(f"{Path(path).name}: {int((~ok).sum())} column(s) without an Entrez ID dropped, "
            f"e.g. {list(df.columns[~ok][:3])}")
    df = df.loc[:, ok]
    df.columns = [e for e, k in zip(ent, ok) if k]
    df = df.loc[:, ~df.columns.duplicated()]
    df.index = pd.Index([str(i) for i in df.index], name="sample")
    return df.apply(pd.to_numeric, errors="coerce").astype(np.float64)


def read_metadata(path) -> pd.DataFrame:
    """Sample x META_COLS (essml.features.build_meta encoding); missing columns -> unknown."""
    from essml import features as F
    df = _read_table(path)
    m = df.reindex(columns=F.META_COLS).rename_axis("ModelID").reset_index()
    m["ModelID"] = m["ModelID"].astype(str)
    out = F.build_meta(m.astype({c: "string" for c in F.META_COLS if c != "Age"}), F.META_COLS)
    out.index = pd.Index(out.index.astype(str), name="sample")
    return out


def _cached_columns(name) -> list[str]:
    import pyarrow.parquet as pq
    path = ROOT / "cache" / (f"features/{name}.parquet" if name in ("C", "M", "H") else f"{name}.parquet")
    return [c for c in pq.read_schema(path).names if c != "ModelID"]


def genomic_blocks(samples, cn=None, damaging=None, hotspot=None, metadata=None, log=print) -> dict:
    """m5_avg's non-expression raw blocks for user samples, built as for DepMap lines (essml.features):
    C = cn_transform(CN) and '_cn' = CN (M2 gene-level), M = '_damaging' = dosage, H = hotspot > 0, meta.
    Columns follow the DepMap caches; a gene missing from a supplied table is diploid (C 0; '_cn' NaN, which M2
    fills as for DepMap NaN) / wild type (0). Samples missing from a table are left out of that block = absent."""
    from essml import features as F
    samples = pd.Index(samples)
    out = {}

    def rows(df, what):
        have = samples.intersection(df.index)
        miss = samples.difference(df.index)
        if len(miss):
            log(f"{what}: {len(miss)} sample(s) not in the table -> absent (diploid wild type), e.g. {list(miss[:3])}")
        return df.loc[have]

    if cn is not None:
        cn = rows(cn, "copy number")
        cols = _cached_columns("cn")
        out["_cn"] = cn.reindex(columns=cols).astype(np.float32)
        x = np.clip(np.log2((cn.to_numpy(np.float32) + F.CN_OFFSET) / (1.0 + F.CN_OFFSET)), -F.CN_CLIP, F.CN_CLIP)
        C = pd.DataFrame(np.nan_to_num(x, nan=0.0).astype(np.float32), index=cn.index, columns=cn.columns)
        out["C"] = C.reindex(columns=_cached_columns("C"), fill_value=0.0)
    if damaging is not None:
        dam = rows(damaging, "damaging mutations")
        out["M"] = dam.reindex(columns=_cached_columns("M"), fill_value=0.0).fillna(0.0).astype(np.float32)
        out["_damaging"] = out["M"]
    if hotspot is not None:
        hot = rows(hotspot, "hotspot mutations")
        out["H"] = F.hotspot_binary(hot.reindex(columns=_cached_columns("H"), fill_value=0.0).fillna(0.0))
    if metadata is not None:
        out["meta"] = rows(metadata, "metadata")
    return out


def with_input_ids(blocks: dict) -> dict:
    """Prefix the sample index with 'input::' (the IDs predict_table / predict_new use for new samples)."""
    return {k: v.set_axis([f"input::{s}" for s in v.index], axis=0) for k, v in blocks.items()}


def dev_reference():
    """DepMap dev expression (cache/features/E.parquet), dev lineages."""
    sp = pd.read_csv(ROOT / "configs" / "splits.csv", index_col=0)
    sp.index = sp.index.astype(str)
    dev = sp.index[sp["role"] == "dev"]
    E = pd.read_parquet(ROOT / "cache" / "features" / "E.parquet")
    E.columns, E.index = E.columns.astype(str), E.index.astype(str)
    return E.loc[E.index.intersection(dev)].astype(np.float64), sp["OncotreeLineage"].astype(str)


def harmonise(expr: pd.DataFrame, pool: pd.DataFrame, method: str) -> pd.DataFrame:
    if method == "none":
        return expr
    if method == "quantile":
        q = XH.reference_quantiles(pool.to_numpy())
        cols = pool.columns
        x = expr.reindex(columns=cols)
        return XH.quantile_normalise(x, q)
    raise ValueError("harmonise must be 'none' or 'quantile'")


def applicability_domain(expr: pd.DataFrame, pool: pd.DataFrame, lineage: pd.Series, k=None, n_hvg=None):
    k, n_hvg = k or AD_K, n_hvg or N_HVG
    hvg = XH.top_variable_genes(pool, n_hvg)
    r = XH.corr_to_pool(expr.reindex(columns=hvg), pool[hvg])
    ad = XH.applicability(r, lineage, k)
    pr = pd.DataFrame(np.corrcoef(pool[hvg].to_numpy()), index=pool.index, columns=pool.index)
    self_ad = XH.pool_self_ad(pr, k)
    ad["ad_percentile_vs_dev"] = [float((self_ad < d).mean() * 100) for d in ad["ad_distance"]]
    return ad


# ============================================================================ model
def train_tag(model_name, cfg, train_ids):
    h = hashlib.sha256(json.dumps([cfg, sorted(map(str, train_ids))], sort_keys=True).encode()).hexdigest()
    return f"predict_{model_name}_{h[:12]}"


def _keep_inputs():
    from essml.models.m5_stack import StackM5
    StackM5.KEEP_INPUTS = True  # M5 input refits serve every sample of the run


def fit_model(model_name, train_ids, log=print):
    """Fit on dev lines ``train_ids`` with the launcher's guarded fit (sentinel + FIT_LOG); dev targets via guard."""
    import run as RUN
    cfg = json.loads(SUPPORTED[model_name]["config"].read_text())
    if cfg.get("grid"):
        raise NotImplementedError("configs with a tuning grid need the launcher's inner CV (not needed for M1)")
    roles = guard.roles_of(train_ids)
    if (roles != "dev").any():
        raise guard.BlindingError("predict.py trains on DepMap dev lines only")
    sp = pd.read_csv(ROOT / "configs" / "splits.csv", index_col=0)
    sp.index = sp.index.astype(str)
    Y = guard.load_targets(list(train_ids)).astype(np.float64)
    excluded = sp.index.difference(pd.Index(train_ids))
    sent = RUN.Sentinel(excluded, "predict.py")
    t0 = time.perf_counter()
    model = RUN.guarded_fit(RUN.MODELS[cfg["model"]](**cfg.get("params", {})), sent, "predict", Y,
                            sp["OncotreeLineage"].astype(str), None, sp["PatientID"].astype(str))
    if cfg["model"] == "M5":  # average mode stores the data; its three inputs are refitted on the dev lines here
        _keep_inputs()
        model.fitted_inputs()
    log(f"fitted {model_name} on {len(Y)} dev lines in {time.perf_counter() - t0:.0f} s")
    return model, cfg


def get_models(model_name, train_ids, save=False, model_dir=modelstore.MODEL_DIR, log=print):
    """[(fold, model, meta)] loaded from cache/models/<tag>/ if persisted, else fitted (and saved if asked)."""
    cfg = json.loads(SUPPORTED[model_name]["config"].read_text())
    tag = train_tag(model_name, cfg, train_ids)
    if (Path(model_dir) / tag / "fold0.pkl").exists():
        log(f"loading persisted model cache/models/{tag}")
        _keep_inputs()
        return modelstore.load_fold_models(tag, model_dir), tag
    model, _ = fit_model(model_name, train_ids, log)
    if save:
        if hasattr(model, "compact"):
            model.compact()
        modelstore.save_fold_model(tag, 0, model, train_ids, [], model_dir)
        log(f"saved cache/models/{tag}")
        return modelstore.load_fold_models(tag, model_dir), tag
    return [(0, model, {"train_ids": list(train_ids)})], tag


# ============================================================================ prediction
def add_notes(long: pd.DataFrame, model_name: str, genomics_given: bool = False) -> pd.DataFrame:
    """``note`` column: flags mutation-driven dependencies (D24; not a tier change) for expression-only models, and
    for m5_avg when no genomic input was supplied (always, in this script)."""
    flagged = long["symbol"].isin(MUTATION_NOTE_GENES)
    if model_name in EXPRESSION_ONLY_MODELS:
        note = np.where(flagged, MUTATION_NOTE, "")
    elif model_name in GENOMIC_INPUT_MODELS and not genomics_given:
        note = np.where(flagged, NO_GENOMICS_NOTE, "")
    else:
        note = np.full(len(long), "", dtype=object)
    return long.assign(note=note)


def predict_table(expr: pd.DataFrame, models, table: pd.DataFrame, x_builder=None, ad: pd.DataFrame | None = None):
    """Mean prediction of ``models`` for the samples of ``expr`` -> tidy table with confidence / tier (+ AD)."""
    x_builder = x_builder or X_BUILDER
    ids = pd.Index([f"input::{s}" for s in expr.index])
    X = x_builder(expr.set_axis(ids, axis=0))
    preds = [m.predict(ids, lineage=None, X=X) for _, m, _ in models]
    wide = sum(preds) / len(preds)
    if not np.isfinite(wide.to_numpy()).all():
        raise ValueError("non-finite predictions (check the input expression)")
    wide.index = expr.index
    wide.index.name = "sample"
    long = PR.annotate_long(wide.rename_axis("ModelID"), table).rename(columns={"ModelID": "sample"})
    if ad is not None:
        long = long.join(ad[["ad_distance", "ad_percentile_vs_dev", "nearest_lineage"]], on="sample")
    return long


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--expression", required=True)
    ap.add_argument("--out", required=True, help=".csv or .parquet")
    ap.add_argument("--model", default=DEFAULT_MODEL, choices=sorted(SUPPORTED))
    ap.add_argument("--harmonise", default="none", choices=["none", "quantile"])
    ap.add_argument("--genes-in-rows", action="store_true")
    ap.add_argument("--exclude-lines", default="", help="comma-separated dev ModelIDs left out of training")
    ap.add_argument("--save-models", action="store_true")
    ap.add_argument("--run-id", default=None, help="average the persisted fold models of this CV run instead")
    ap.add_argument("--predictable-genes", default=None, help="tier table (default: the model's own)")
    ap.add_argument("--no-ad", action="store_true")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--cn", default=None, help="m5_avg: gene-level relative copy number (linear; DepMap layout ok)")
    ap.add_argument("--damaging", default=None, help="m5_avg: damaging-mutation dosage 0/1/2")
    ap.add_argument("--hotspot", default=None, help="m5_avg: hotspot-mutation matrix 0/1/2")
    ap.add_argument("--metadata", default=None, help="m5_avg: OncotreeLineage, OncotreePrimaryDisease, Sex, Age")
    a = ap.parse_args(argv)
    given = [o for o in GENOMIC_OPTIONS if getattr(a, o)]
    if given and a.model not in GENOMIC_INPUT_MODELS:
        ap.error(f"{a.model} reads expression only; --{', --'.join(given)} apply to m5_avg")
    try:
        from threadpoolctl import threadpool_limits
        threadpool_limits(a.threads)
    except ImportError:
        pass
    expr = read_expression(a.expression, a.genes_in_rows)
    pool, lineage = dev_reference()
    excl = [s.strip() for s in a.exclude_lines.split(",") if s.strip()]
    train_ids = pool.index.difference(pd.Index(excl))  # dev lines with expression (= all 1022)
    pool_ad = pool.loc[train_ids]                      # reference pool excludes held-out lines too
    expr_h = harmonise(expr, pool_ad, a.harmonise)
    cov = float(expr.columns.isin(pool.columns).mean())
    n_cov = int(pool.columns.isin(expr.columns).sum())
    print(f"{len(expr)} samples; {n_cov}/{pool.shape[1]} model genes present in the input "
          f"(absent genes -> training mean); harmonise={a.harmonise}")
    ad = None if a.no_ad else applicability_domain(expr_h, pool_ad, lineage)
    G = {}
    if given:  # genomics are not harmonised: they must follow DepMap's calling conventions (docstring)
        G = with_input_ids(genomic_blocks(
            expr.index, *[read_genomic(getattr(a, o), a.genes_in_rows) if getattr(a, o) else None
                          for o in ("cn", "damaging", "hotspot")],
            read_metadata(a.metadata) if a.metadata else None))
    x_builder = (lambda E: {**X_BUILDER(E), **G}) if G else X_BUILDER
    table_path = a.predictable_genes or (a.run_id and PR.table_for_run(a.run_id)) or SUPPORTED[a.model]["table"]
    table = PR.load_table(table_path)
    if a.model in EXPRESSION_ONLY_MODELS:
        print(M1_WARNING)
    if a.model in GENOMIC_INPUT_MODELS and not given:
        print(NO_GENOMICS_WARNING)
    elif given:
        absent = [o for o in GENOMIC_OPTIONS if o not in given]
        print(f"m5_avg with user-supplied {', '.join(given)}" + (
            f"; not supplied (diploid wild type / unknown): {', '.join(absent)}" if absent else ""))
    if a.run_id:
        wide = modelstore.predict_new(a.run_id, [f"input::{s}" for s in expr.index],
                                      X=x_builder(expr_h.set_axis([f"input::{s}" for s in expr.index], axis=0)),
                                      output="wide")
        wide.index = expr.index
        long = PR.annotate_long(wide.rename_axis("ModelID"), table).rename(columns={"ModelID": "sample"})
        if ad is not None:
            long = long.join(ad[["ad_distance", "ad_percentile_vs_dev", "nearest_lineage"]], on="sample")
        tag = a.run_id
    else:
        models, tag = get_models(a.model, train_ids, a.save_models)
        long = predict_table(expr_h, models, table, x_builder=x_builder, ad=ad)
    long = add_notes(long, a.model, genomics_given=bool(a.hotspot))  # D24 recovery runs through hotspots
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.suffix == ".parquet":
        long.to_parquet(out)
    else:
        long.to_csv(out, index=False)
    meta = {"model": a.model, "model_tag": tag,
            "genomic_inputs": ({o: getattr(a, o) for o in given} if given else
                               "none supplied: diploid wild type / unknown metadata (DESIGN §2 convention)"
                               if a.model in GENOMIC_INPUT_MODELS else "not used by this model"),
            "genomic_samples": {k[1:] if k.startswith("_") else k: int(len(v)) for k, v in G.items()},
            "evaluation": ("all four genomic inputs: the evaluated headline (dev CV, G4 test); block path checked in "
                           "results/predict_genomics_check.json" if len(given) == len(GENOMIC_OPTIONS) else
                           "post hoc, dev CV only: results/predict_expression_only_cv.json" if not given else
                           "partial genomic inputs: not evaluated as such" if a.model in GENOMIC_INPUT_MODELS
                           else "G4 test set (expression-only model)"),
            "n_training_lines": int(len(train_ids)), "excluded_lines": excl,
            "harmonise": a.harmonise, "input_gene_fraction_recognised": cov, "model_genes_present": n_cov,
            "predictable_genes": str(table_path), "predictable_genes_sha256": PR.sha256_file(table_path),
            "applicability": None if ad is None else ad.reset_index().to_dict(orient="records")}
    out.with_suffix(".meta.json").write_text(json.dumps(meta, indent=1, default=str))
    print(f"wrote {out} ({len(long)} rows) and {out.with_suffix('.meta.json')}")
    return long


if __name__ == "__main__":
    main()
