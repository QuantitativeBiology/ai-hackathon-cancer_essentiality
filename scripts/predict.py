"""Predict CRISPR gene effects (Chronos scale) for new cell lines from their expression profile, with a
per-gene confidence score (DESIGN §4 "Per-gene predictability", D22).

    python scripts/predict.py --expression my_lines.csv --out predictions.csv
    python scripts/predict.py --expression ext.tsv --harmonise quantile --out preds.parquet   # external profiles
    python scripts/predict.py --expression x.csv --save-models          # persist the fitted model for next time
    python scripts/predict.py --expression x.csv --run-id <cv run id>   # use a run's persisted fold models

Input (``--expression``): a csv / tsv / parquet table of log2(TPM + 1) on DepMap's protein-coding genes, one
row per sample (first column = sample id). Gene columns are Entrez IDs ("7105") or DepMap headers
("TSPAN6 (7105)"); ``--genes-in-rows`` reads the transposed layout. Genes absent from the input are set to
the training mean by the model's feature code. External (non-DepMap) profiles should be harmonised with
``--harmonise quantile`` (the frozen D9 choice: essml.expr_harmonise.quantile_normalise to the mean sorted
DepMap dev profile). TPM must be computed over all genes before subsetting (DESIGN §7.3; see
essml.expr_harmonise for count / FPKM conversion).

Model: ``m1_E_resB0`` (expression-only kernel ridge, B0 residual; an additional G4 finalist). It is fitted on
the DepMap dev lines (minus ``--exclude-lines``) with the launcher's sentinel, or loaded from
cache/models/<tag>/ if it was persisted earlier with ``--save-models``; ``--run-id`` instead averages a CV
run's persisted fold models (scripts/run.py --save-models) via essml.modelstore.predict_new.
TODO(m5_avg): the G4 headline m5_avg averages m1_E_resB0, M2 (paralog priors) and M3 (LightGBM on residuals);
M2/M3 need gene-level features built from the new line's expression and a nested fit protocol, which this
script does not implement yet. Until then predictions come from m1_E_resB0 (dev-CV P1 0.308 vs 0.3175).

Output (``--out``, csv or parquet): one row per (sample, gene): sample, entrez, symbol, predicted_effect,
confidence, tier, ad_distance, ad_percentile_vs_dev, nearest_lineage. ``confidence`` / ``tier`` come from
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

SUPPORTED = {"m1_E_resB0": {"config": ROOT / "configs" / "runs" / "m1_E_resB0.json",
                            "table": ROOT / "configs" / "predictable_genes_m1.csv"}}
AD_K = 10
N_HVG = 5000
X_BUILDER = (lambda E: {"E": E})  # raw expression block for the kernel models
_HDR = re.compile(r"^.*\((\d+)\)\s*$")


# ============================================================================ input
def read_expression(path, genes_in_rows=False) -> pd.DataFrame:
    """Samples x Entrez-string columns, float64."""
    p = Path(path)
    if p.suffix == ".parquet":
        df = pd.read_parquet(p)
    else:
        df = pd.read_csv(p, sep="\t" if p.suffix in (".tsv", ".txt") else ",", index_col=0)
    if genes_in_rows:
        df = df.T
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
    log(f"fitted {model_name} on {len(Y)} dev lines in {time.perf_counter() - t0:.0f} s")
    return model, cfg


def get_models(model_name, train_ids, save=False, model_dir=modelstore.MODEL_DIR, log=print):
    """[(fold, model, meta)] loaded from cache/models/<tag>/ if persisted, else fitted (and saved if asked)."""
    cfg = json.loads(SUPPORTED[model_name]["config"].read_text())
    tag = train_tag(model_name, cfg, train_ids)
    if (Path(model_dir) / tag / "fold0.pkl").exists():
        log(f"loading persisted model cache/models/{tag}")
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
def predict_table(expr: pd.DataFrame, models, table: pd.DataFrame, x_builder=None, ad: pd.DataFrame | None = None):
    """Mean prediction of ``models`` for the samples of ``expr`` -> tidy table with confidence / tier (+ AD)."""
    x_builder = x_builder or X_BUILDER
    ids = pd.Index([f"input::{s}" for s in expr.index])
    X = x_builder(expr.set_axis(ids, axis=0))
    preds = [m.predict(ids, lineage=None, X=X) for _, m, _ in models]
    wide = sum(preds) / len(preds)
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
    ap.add_argument("--model", default="m1_E_resB0", choices=sorted(SUPPORTED) + ["m5_avg"])
    ap.add_argument("--harmonise", default="none", choices=["none", "quantile"])
    ap.add_argument("--genes-in-rows", action="store_true")
    ap.add_argument("--exclude-lines", default="", help="comma-separated dev ModelIDs left out of training")
    ap.add_argument("--save-models", action="store_true")
    ap.add_argument("--run-id", default=None, help="average the persisted fold models of this CV run instead")
    ap.add_argument("--predictable-genes", default=None, help="tier table (default: the model's own)")
    ap.add_argument("--no-ad", action="store_true")
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args(argv)
    if a.model == "m5_avg":
        raise SystemExit("m5_avg is not supported by predict.py yet (TODO in the module docstring); use m1_E_resB0")
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
    table = PR.load_table(a.predictable_genes or SUPPORTED[a.model]["table"])
    if a.run_id:
        wide = modelstore.predict_new(a.run_id, [f"input::{s}" for s in expr.index],
                                      X=X_BUILDER(expr_h.set_axis([f"input::{s}" for s in expr.index], axis=0)),
                                      output="wide")
        wide.index = expr.index
        long = PR.annotate_long(wide.rename_axis("ModelID"), table).rename(columns={"ModelID": "sample"})
        if ad is not None:
            long = long.join(ad[["ad_distance", "ad_percentile_vs_dev", "nearest_lineage"]], on="sample")
        tag = a.run_id
    else:
        models, tag = get_models(a.model, train_ids, a.save_models)
        long = predict_table(expr_h, models, table, ad=ad)
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.suffix == ".parquet":
        long.to_parquet(out)
    else:
        long.to_csv(out, index=False)
    meta = {"model": a.model, "model_tag": tag, "n_training_lines": int(len(train_ids)), "excluded_lines": excl,
            "harmonise": a.harmonise, "input_gene_fraction_recognised": cov, "model_genes_present": n_cov,
            "predictable_genes": str(a.predictable_genes or SUPPORTED[a.model]["table"]),
            "predictable_genes_sha256": PR.sha256_file(a.predictable_genes or SUPPORTED[a.model]["table"]),
            "applicability": None if ad is None else ad.reset_index().to_dict(orient="records")}
    out.with_suffix(".meta.json").write_text(json.dumps(meta, indent=1, default=str))
    print(f"wrote {out} ({len(long)} rows) and {out.with_suffix('.meta.json')}")
    return long


if __name__ == "__main__":
    main()
