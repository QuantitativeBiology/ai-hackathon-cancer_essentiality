"""Minimal sanity check of the external screens in configs/external_screens.json (+ extracted ORCS screens).

For each screen arm: number of guides/genes, the sample columns present, and a positive/negative
control check - median gene score of Hart CEGv2 core-essential genes vs NEGv1 non-essential genes.
  * counts / RPM:  crude log2 fold change T_end vs T0 per guide:
        lfc = log2((end + 0.5)/(sum_end + 1) * 1e6) - log2((t0 + 0.5)/(sum_t0 + 1) * 1e6),
        replicate columns are averaged after per-column CPM normalisation; guides with mean raw
        T0 < MIN_T0 reads (or < MIN_T0_RPM for already-normalised tables) are dropped; gene = mean of guides.
        No centring, no scaling - this is a QC, not a harmonised score.
  * guide LFC tables: gene = mean guide LFC.  * gene-level tables: the published score.
"essential_is" says whether essential genes are expected at low (LFC, CS, phenotype) or high (BF) scores.
AUROC = P(random CEG gene is more essential than random NEG gene); pass = AUROC >= 0.75.
Data files are read as data only. Output: results/external_screens/sanity.csv.
"""
import argparse
import gzip
import io
import json
import pathlib
import re
import zipfile

import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
MIN_T0 = 30
MIN_T0_RPM = 5.0
ORCS_DIRECTION = {"Bayes Factor": "high", "CRISPR Score (CS)": "low", "Log2FC": "low", "Log2": "low",
                  "Z-score": "low", "ZLFC": "low", "mean fold change": "low", "Beta Score": "low",
                  "FS (Fitness Score)": "low", "Essentiality Score": "low", "STARS Score": "high",
                  "Mean Depletion": "low"}


# ---------------------------------------------------------------- gene symbol normalisation
def load_hgnc():
    h = pd.read_csv(ROOT / "data/raw_external/hgnc/hgnc_complete_set.txt", sep="\t", dtype=str,
                    usecols=["hgnc_id", "symbol", "alias_symbol", "prev_symbol", "entrez_id"])
    approved = set(h["symbol"])
    prev, alias = {}, {}
    for r in h.itertuples(index=False):
        for col, d in ((r.prev_symbol, prev), (r.alias_symbol, alias)):
            if isinstance(col, str):
                for s in col.split("|"):
                    d.setdefault(s, set()).add(r.symbol)
    by_id = dict(zip(h["hgnc_id"], h["symbol"]))
    by_entrez = {e: s for e, s in zip(h["entrez_id"], h["symbol"]) if isinstance(e, str)}

    def norm(sym):
        s = str(sym).strip()
        if s in approved:
            return s
        for d in (prev, alias):
            if s in d and len(d[s]) == 1:
                return next(iter(d[s]))
        return s
    return norm, by_id, by_entrez


def load_controls(norm, by_id):
    def read(p):
        t = pd.read_csv(p, sep="\t", dtype=str)
        syms = {by_id.get(i, norm(g)) for g, i in zip(t["GENE"], t["HGNC_ID"])}
        return syms, set(t["ENTREZ_ID"].dropna())
    ceg = read(ROOT / "data/raw_external/hart_ceg2/CEGv2.txt")
    neg = read(ROOT / "data/raw_external/hart_neg/NEGv1.txt")
    return ceg, neg


# ---------------------------------------------------------------- readers
def _bytes_from_zip(path, member, inner_zip=None):
    z = zipfile.ZipFile(path)
    if inner_zip:
        z = zipfile.ZipFile(io.BytesIO(z.read(inner_zip)))
    return z.read(member)


def _text_table(raw, sep):
    # universal newlines: some deposits use bare CR line endings
    return pd.read_csv(io.StringIO(raw.decode("latin1").replace("\r\n", "\n").replace("\r", "\n")), sep=sep)


def read_spec(d, s):
    f = d / s["file"]
    fmt = s["format"]
    if fmt == "table":
        raw = gzip.open(f).read() if f.name.endswith(".gz") else f.read_bytes()
        return _text_table(raw, s["sep"])
    if fmt in ("zip_table", "zip_lfc", "zip_gene"):
        return _text_table(_bytes_from_zip(f, s["member"], s.get("inner_zip")), s["sep"])
    if fmt in ("excel", "gene_excel", "mageck_sgrna_excel"):
        src = io.BytesIO(gzip.open(f).read()) if s.get("gzip") else f
        return pd.read_excel(src, sheet_name=s.get("sheet", 0), header=s.get("header", 0))
    if fmt == "excel_multi":
        return pd.concat([pd.read_excel(f, sheet_name=sh, header=s.get("header", 0)) for sh in s["sheets"]], ignore_index=True)
    if fmt == "zip_excel_lfc":
        return pd.read_excel(io.BytesIO(_bytes_from_zip(f, s["member"])), sheet_name=s["sheet"], header=s["header"])
    if fmt == "zip_excel_mapped":
        b = io.BytesIO(_bytes_from_zip(f, s["member"]))
        x = pd.read_excel(b, sheet_name=s["sheet"], header=s["header"])
        mp = pd.read_excel(b, sheet_name=s["map_sheet"], header=s["map_header"])
        x = x.merge(mp[[s["map_guide_col"], s["map_gene_col"]]], left_on=s["guide_col"],
                    right_on=s["map_guide_col"], how="left")
        return x.rename(columns={s["map_gene_col"]: "__gene"})
    raise ValueError(fmt)


# ---------------------------------------------------------------- scoring
def crude_lfc(x, t0_cols, end_cols, normalised):
    t0 = x[t0_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    en = x[end_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    cpm = lambda m: (m + 0.5).div(m.sum(axis=0) + 1, axis=1) * 1e6
    keep = t0.mean(axis=1) >= (MIN_T0_RPM if normalised else MIN_T0)
    lfc = np.log2(cpm(en).mean(axis=1)) - np.log2(cpm(t0).mean(axis=1))
    return lfc[keep], int(keep.sum())


def auroc(ess, non, direction):
    a, b = np.asarray(ess, float), np.asarray(non, float)
    a, b = a[~np.isnan(a)], b[~np.isnan(b)]
    if len(a) == 0 or len(b) == 0:
        return np.nan
    if direction == "high":
        a, b = -a, -b
    r = pd.Series(np.concatenate([a, b])).rank().to_numpy()
    u = r[: len(a)].sum() - len(a) * (len(a) + 1) / 2
    return 1 - u / (len(a) * len(b))  # P(ess < non) after orienting so that low = essential


def summarise(gene_scores, ceg, neg, direction, by="symbol"):
    gs = gene_scores.dropna()
    sets = (ceg[0], neg[0]) if by == "symbol" else (ceg[1], neg[1])
    e = gs[gs.index.isin(sets[0])]
    n = gs[gs.index.isin(sets[1])]
    au = auroc(e, n, direction) if direction in ("low", "high") else np.nan
    return {"n_genes_scored": len(gs), "n_ceg": len(e), "n_neg": len(n),
            "median_ceg": round(float(e.median()), 3) if len(e) else np.nan,
            "median_neg": round(float(n.median()), 3) if len(n) else np.nan,
            "auroc_ceg_vs_neg": round(float(au), 3) if not np.isnan(au) else np.nan,
            "pass": (bool(au >= 0.75) if not np.isnan(au) else None)}


def check_config_screen(sc, norm, ceg, neg):
    d = ROOT / "data/external_screens" / sc["deposit"]
    s = sc["sanity"]
    x = read_spec(d, s)
    x.columns = [str(c).strip() if not str(c).startswith("Unnamed") else str(c) for c in x.columns]
    base = {"screen_id": sc["id"], "category": sc["category"], "library": sc["library"],
            "data_level": sc["data_level"], "n_rows": len(x), "note": s.get("contrast_note", "")}
    rows = []
    if s["format"] == "zip_excel_mapped":
        genes = x["__gene"]
    elif "gene_from_guide_regex" in s:
        rx = re.compile(s["gene_from_guide_regex"])
        genes = x[s["guide_col"]].astype(str).map(lambda g: (rx.search(g).group(1) if rx.search(g) else np.nan))
    else:
        genes = x[s["gene_col"]]
    genes = genes.map(lambda g: norm(g) if isinstance(g, str) else np.nan)
    base["n_genes"] = int(genes.nunique())
    fmt = s["format"]
    if fmt in ("table", "zip_table", "excel", "excel_multi", "zip_excel_mapped"):
        sample_cols = [c for c in x.columns if c not in (s.get("guide_col"), s.get("gene_col"), "__gene", s.get("map_guide_col"))]
        base["columns"] = ";".join(map(str, sample_cols))[:400]
        for arm, (t0c, endc) in s["arms"].items():
            lfc, used = crude_lfc(x, t0c, endc, s.get("is_normalised", False))
            g = lfc.groupby(genes[lfc.index]).mean()
            rows.append({**base, "arm": arm, "t0": "+".join(t0c), "t_end": "+".join(endc),
                         "n_guides_used": used, "score": "crude LFC", **summarise(g, ceg, neg, "low")})
    elif fmt in ("zip_lfc", "zip_excel_lfc"):
        base["columns"] = ";".join(s["lfc_cols"])
        cols = s["lfc_cols"] if fmt == "zip_lfc" else [s["lfc_cols"]]
        for c in cols:
            v = x[c].apply(pd.to_numeric, errors="coerce")
            v = v.mean(axis=1) if isinstance(v, pd.DataFrame) else v
            g = v.groupby(genes).mean()
            rows.append({**base, "arm": c if isinstance(c, str) else "mean_reps", "t0": "pDNA/T0 (published LFC)",
                         "t_end": "", "n_guides_used": int(v.notna().sum()), "score": "published guide LFC",
                         **summarise(g, ceg, neg, "low")})
    elif fmt == "mageck_sgrna_excel":
        base["columns"] = "control_mean;treat_mean (MAGeCK normalised)"
        lfc = np.log2((x["treat_mean"].astype(float) + 0.5) / (x["control_mean"].astype(float) + 0.5))
        g = lfc.groupby(genes).mean()
        rows.append({**base, "arm": s["sheet"], "t0": "control (non-dividing)", "t_end": "treatment (dividing)",
                     "n_guides_used": int(lfc.notna().sum()), "score": "LFC of MAGeCK means",
                     **summarise(g, ceg, neg, "low")})
    elif fmt in ("gene_excel", "zip_gene"):
        if "sheets_each" in s:  # several sheets, one score column each
            for sh in s["sheets_each"]:
                y = pd.read_excel(d / s["file"], sheet_name=sh, header=s.get("header", 0))
                gn = y[s["gene_col"]].map(lambda q: norm(q) if isinstance(q, str) else np.nan)
                g = pd.to_numeric(y[s["score_col"]], errors="coerce").groupby(gn).mean()
                rows.append({**base, "n_rows": len(y), "n_genes": int(gn.nunique()), "columns": s["score_col"],
                             "arm": sh, "t0": "", "t_end": "", "n_guides_used": np.nan,
                             "score": s["score_col"], **summarise(g, ceg, neg, s["essential_is"])})
        else:
            for c in s.get("score_cols", [s.get("score_col")]):
                g = pd.to_numeric(x[c], errors="coerce").groupby(genes).mean()
                rows.append({**base, "columns": ";".join(s.get("score_cols", [s.get("score_col")])), "arm": c,
                             "t0": "", "t_end": "", "n_guides_used": np.nan, "score": c,
                             **summarise(g, ceg, neg, s["essential_is"])})
    return rows


def check_orcs(ceg, neg):
    res = ROOT / "results/external_screens"
    ext = ROOT / "data/external_screens/biogrid_orcs/extracted"
    rows = []
    for tab, cat_col in ((res / "orcs_nontransformed_screens.csv", "category"), (res / "orcs_anchor_screens.csv", "category")):
        if not tab.exists():
            continue
        meta = pd.read_csv(tab)
        for r in meta.to_dict("records"):
            fn = r["extracted_file"]
            if not isinstance(fn, str) or not (ext / fn).exists():
                continue
            x = pd.read_csv(ext / fn, sep="	", dtype={"IDENTIFIER_ID": str}, low_memory=False)
            sc_type = r["SCORE.1_TYPE"]
            v = pd.to_numeric(x["SCORE.1"], errors="coerce")
            g = v.groupby(x["IDENTIFIER_ID"]).mean()
            direction = ORCS_DIRECTION.get(sc_type, "unknown")
            rows.append({"screen_id": f"ORCS_{r['#SCREEN_ID']}", "category": f"ORCS {r[cat_col]}",
                         "library": r["LIBRARY"], "data_level": f"gene-level {sc_type}", "n_rows": len(x),
                         "note": f"{r['AUTHOR']}; {r['CELL_LINE']}; {r['PHENOTYPE']}; {r['EXPERIMENTAL_SETUP']}; {r['CONDITION_NAME']}; {r['DURATION']}",
                         "n_genes": int(x["IDENTIFIER_ID"].nunique()), "columns": "SCORE.1..5",
                         "arm": r["CELL_LINE"], "t0": "", "t_end": "", "n_guides_used": np.nan,
                         "score": f"SCORE.1 ({sc_type}; direction {direction})",
                         **summarise(g, ceg, neg, direction, by="entrez")})
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=str(ROOT / "configs/external_screens.json"))
    ap.add_argument("--out", default=str(ROOT / "results/external_screens/sanity.csv"))
    ap.add_argument("--no-orcs", action="store_true")
    a = ap.parse_args()
    norm, by_id, _ = load_hgnc()
    ceg, neg = load_controls(norm, by_id)
    cfg = json.load(open(a.config, encoding="utf-8"))
    rows = []
    for sc in cfg["screens"]:
        try:
            rows += check_config_screen(sc, norm, ceg, neg)
        except Exception as e:  # report and continue
            rows.append({"screen_id": sc["id"], "category": sc["category"], "note": f"ERROR {e!r}"})
            print("ERROR", sc["id"], repr(e))
    if not a.no_orcs:
        rows += check_orcs(ceg, neg)
    out = pd.DataFrame(rows)
    cols = ["screen_id", "arm", "category", "library", "data_level", "score", "n_rows", "n_genes",
            "n_guides_used", "t0", "t_end", "n_genes_scored", "n_ceg", "n_neg", "median_ceg", "median_neg",
            "auroc_ceg_vs_neg", "pass", "note", "columns"]
    out = out.reindex(columns=cols)
    pathlib.Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(a.out, index=False)
    with pd.option_context("display.width", 250, "display.max_rows", 500, "display.max_colwidth", 30):
        print(out[["screen_id", "arm", "n_rows", "n_genes", "n_guides_used", "n_ceg", "n_neg",
                   "median_ceg", "median_neg", "auroc_ceg_vs_neg", "pass"]].to_string())


if __name__ == "__main__":
    main()
