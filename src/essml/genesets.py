"""Gene classes for evaluation (FROZEN: hashed in configs/frozen_hashes.json; DESIGN.md §4).

Classes are computed from TRAINING rows only (a CV training fold, or all dev lines for the test).
All statistics ignore NaN: SDs and prevalences are over observed (non-NaN) lines of each gene.

- all              every gene in Y_train
- variable         SD > 0.25 (ddof=1)
- selective        5% <= fraction of observed lines with effect < -0.5 <= 95%
- rare_selective   >= 3 observed lines < -0.5 and fraction < 5%
- common_essential fraction < -0.5 > 95%
- nonessential     Hart NEGv1, mapped to Entrez via HGNC (see load_hart_neg)

Genes with no observed training line belong only to 'all' (and possibly 'nonessential').
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
NEG_PATH = ROOT / "data" / "raw_external" / "hart_neg" / "NEGv1.txt"
HGNC_PATH = ROOT / "data" / "raw_external" / "hgnc" / "hgnc_complete_set.txt"
CLASSES = ["all", "variable", "selective", "rare_selective", "common_essential", "nonessential"]


def load_hart_neg(neg_path=NEG_PATH, hgnc_path=HGNC_PATH):
    """Hart NEGv1 -> (set of Entrez ID strings, report dict).

    Mapping order per NEGv1 row: (1) HGNC_ID -> HGNC entrez_id; (2) approved symbol -> entrez_id;
    (3) previous/alias symbol if it maps to exactly one HGNC entry. The Entrez column shipped in
    NEGv1 is used only as a cross-check (disagreements are counted in the report).
    """
    neg = pd.read_csv(neg_path, sep="\t", dtype=str)
    h = pd.read_csv(hgnc_path, sep="\t", dtype=str, usecols=["hgnc_id", "symbol", "entrez_id",
                                                          "prev_symbol", "alias_symbol"])
    h = h[h["entrez_id"].notna()]
    by_id = dict(zip(h["hgnc_id"], h["entrez_id"]))
    by_sym = dict(zip(h["symbol"], h["entrez_id"]))
    other = {}
    for col in ["prev_symbol", "alias_symbol"]:
        for syms, e in zip(h[col], h["entrez_id"]):
            if isinstance(syms, str):
                for s in syms.split("|"):
                    other.setdefault(s.strip(), set()).add(e)
    out, how, mismatch = set(), {"hgnc_id": 0, "symbol": 0, "prev_alias": 0, "unmapped": 0}, 0
    unmapped = []
    for _, r in neg.iterrows():
        e = by_id.get(r.get("HGNC_ID"))
        if e is not None:
            how["hgnc_id"] += 1
        elif r["GENE"] in by_sym:
            e = by_sym[r["GENE"]]
            how["symbol"] += 1
        elif len(other.get(r["GENE"], ())) == 1:
            e = next(iter(other[r["GENE"]]))
            how["prev_alias"] += 1
        else:
            how["unmapped"] += 1
            unmapped.append(r["GENE"])
            continue
        if isinstance(r.get("ENTREZ_ID"), str) and r["ENTREZ_ID"] != e:
            mismatch += 1
        out.add(str(e))
    rep = {"n_neg": len(neg), "n_mapped": len(out), "by": how, "unmapped": unmapped,
           "entrez_disagreements_with_file": mismatch}
    return out, rep


def gene_classes(Y_train: pd.DataFrame, neg_entrez=None, thr=-0.5, sd_thr=0.25) -> pd.DataFrame:
    """Boolean class table (index = genes of Y_train, columns = CLASSES) from training rows only.

    neg_entrez: iterable of Entrez strings for 'nonessential'; None -> load Hart NEGv1 via HGNC.
    """
    if neg_entrez is None:
        neg_entrez, _ = load_hart_neg()
    Y = np.asarray(Y_train, dtype=np.float64)
    obs = ~np.isnan(Y)
    n = obs.sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        n_below = ((Y < thr) & obs).sum(0)
        frac = n_below / n
        mu = np.where(obs, Y, 0).sum(0) / n
        sd = np.sqrt(np.where(obs, (Y - mu) ** 2, 0).sum(0) / (n - 1))
    has = n > 0
    sd_ok = n > 1
    out = pd.DataFrame(index=Y_train.columns)
    out["all"] = True
    out["variable"] = sd_ok & (sd > sd_thr)
    out["selective"] = has & (frac >= 0.05) & (frac <= 0.95)
    out["rare_selective"] = has & (n_below >= 3) & (frac < 0.05)
    out["common_essential"] = has & (frac > 0.95)
    out["nonessential"] = out.index.astype(str).isin(set(map(str, neg_entrez)))
    return out
