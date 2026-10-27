"""Persist the fitted fold models of a run so they can predict new feature rows (set A / set B).

Layout: cache/models/<run_id>/fold<k>.pkl (gitignored), one pickled dict per outer dev fold:
    {"model": fitted model, "fold": k, "train_ids": [...], "features": [...], "run_id": ...}
Large training-target arrays (e.g. B2's neighbour profiles, attribute ``Y_``) are not pickled; they
are re-read from dev targets through ``guard.load_targets`` (dev IDs only) on load.

``predict_new(run_id, index, lineage, X)`` averages the fold models' predictions (DESIGN v3 §5:
set A/B predictions use the mean of the 5 dev fold-models). ``predict_anchor(run_id, model_id, ...)``
uses only the fold models whose training set excludes that anchor line.
"""
from __future__ import annotations

import pickle
from pathlib import Path

import pandas as pd

from . import guard

ROOT = Path(__file__).resolve().parents[2]
MODEL_DIR = ROOT / "cache" / "models"


def save_fold_model(run_id, fold, model, train_ids, features, model_dir=MODEL_DIR):
    d = Path(model_dir) / run_id
    d.mkdir(parents=True, exist_ok=True)
    state = dict(model.__dict__)
    dropped = []
    if "Y_" in state:  # re-loadable from dev targets; avoid ~100 MB per fold on disk
        state["Y_"] = None
        dropped.append("Y_")
    obj = {"cls": type(model), "state": state, "dropped": dropped, "fold": fold,
           "train_ids": list(map(str, train_ids)), "features": list(features), "run_id": run_id}
    with open(d / f"fold{fold}.pkl", "wb") as fh:
        pickle.dump(obj, fh, protocol=pickle.HIGHEST_PROTOCOL)


def load_fold_models(run_id, model_dir=MODEL_DIR, target_loader=None):
    """Return list of (fold, model, meta) for a run, with dropped arrays restored."""
    target_loader = target_loader or guard.load_targets
    out = []
    for p in sorted((Path(model_dir) / run_id).glob("fold*.pkl")):
        with open(p, "rb") as fh:
            obj = pickle.load(fh)
        m = obj["cls"].__new__(obj["cls"])
        m.__dict__.update(obj["state"])
        if "Y_" in obj["dropped"]:
            Y = target_loader(obj["train_ids"])
            m.Y_ = Y.loc[:, m.genes_].to_numpy("float64")
        out.append((obj["fold"], m, obj))
    if not out:
        raise FileNotFoundError(f"no fold models for run {run_id}")
    return out


def predict_new(run_id, index, lineage=None, X=None, model_dir=MODEL_DIR, target_loader=None,
                per_model=False, output="tidy", predictable_genes=None):
    """Mean prediction of the run's fold models for new lines (index) given lineage / features.

    output="tidy" (default, user-facing): one row per (line, gene) with the per-gene ``confidence`` and
    ``tier`` of configs/predictable_genes.csv (DESIGN §4 "Per-gene predictability"; another table via
    ``predictable_genes``); columns ModelID, entrez, symbol, predicted_effect, confidence, tier.
    output="wide": the [line x gene] matrix (internal use, e.g. ensembles).
    per_model=True returns {fold: wide prediction} instead (for ensemble-size-matched comparisons).
    """
    if output not in ("tidy", "wide"):
        raise ValueError("output must be 'tidy' or 'wide'")
    preds = {f: m.predict(index, lineage=lineage, X=X)
             for f, m, _ in load_fold_models(run_id, model_dir, target_loader)}
    if per_model:
        return preds
    mean = sum(preds.values()) / len(preds)
    if output == "wide":
        return mean
    from . import predictability as PR
    return PR.annotate_long(mean, PR.load_table(predictable_genes or PR.TABLE))


def predict_anchor(run_id, model_id, lineage=None, X=None, model_dir=MODEL_DIR,
                   target_loader=None):
    """Prediction for an anchor (a dev line) from the fold models NOT trained on it, averaged.

    With 5-fold CV a dev line is in the training set of 4 fold models and held out of exactly 1,
    so this is the single out-of-fold model (the line's own CV prediction). Returns
    (prediction DataFrame [1 x gene], list of folds used). Raises if every model saw the line.
    """
    mid = str(model_id)
    use = [(f, m) for f, m, meta in load_fold_models(run_id, model_dir, target_loader)
           if mid not in set(meta["train_ids"])]
    if not use:
        raise ValueError(f"every fold model of {run_id} was trained on {mid}")
    preds = [m.predict([mid], lineage=lineage, X=X) for _, m in use]
    return sum(preds) / len(preds), [f for f, _ in use]
