# predict.py with user-supplied genomics (post hoc, 2026-10-28, after G7)

**Status:** this implements REPORT §7 item 4, added after the report was final. It is not a new model and not a
new evaluation. With all four genomic inputs, the script computes the evaluated headline `m5_avg`: dev-CV P1
0.3175, test P1 0.319, NRAS hotspot AUPRC 0.631 in dev CV.

## Usage
```
python scripts/predict.py --expression expr.csv --cn cn.csv --damaging damaging.csv --hotspot hotspot.csv \
    --metadata meta.csv --out preds.csv
```
- **Table layout:** each table has one row per sample, and its IDs must match the expression table. Gene columns are
  Entrez IDs or DepMap headers ("NRAS (4893)").
- **DepMap release files work as they are.** Identifier columns are recognised, rows are restricted to
  `IsDefaultEntryForModel` and indexed by ModelID. Columns without an Entrez ID are dropped, as in the DepMap cache.
- **`--cn`:** gene-level relative copy number on a linear scale (OmicsCNGeneWGS; 1 = the line's ploidy).
- **`--damaging`:** damaging-mutation dosage 0 / 1 / 2 (OmicsSomaticMutationsMatrixDamaging; LikelyLoF).
- **`--hotspot`:** hotspot matrix 0 / 1 / 2 (OmicsSomaticMutationsMatrixHotspot); any value > 0 counts.
- **`--metadata`:** OncotreeLineage, OncotreePrimaryDisease, Sex and Age, using Model.csv column names and
  vocabulary.
- **Not supplied → absent.** A modality not supplied, or a sample missing from its table, is absent: diploid wild
  type, metadata unknown (DESIGN §2). A gene missing from a supplied table is diploid / wild type.
- **Fusions are not read:** none of `m5_avg`'s three inputs uses them.
- **Which model:** the options apply to `m5_avg` only. `--model m1_E_resB0` refuses them.
- **The RAS/RAF/PIK3CA note** in the output is dropped when `--hotspot` is given.

## Caveats
- **Pipelines must match.** Calls must follow DepMap's pipelines (Mutect2 with LikelyLoF / hotspot annotation;
  WGS relative CN). Genomics from other pipelines are not harmonised, and the effect of a mismatch was not measured.
- **Partial inputs were not evaluated as such.** The evaluated configurations are all four inputs (the headline)
  and none (expression only, post hoc in dev CV: P1 0.3135; docs/predict_expression_only.md).
- **No transfer result applies to this path.** The G6 transfer analysis used expression only for external cells.

## Reproduction check (dev lines only)
`scripts/predict_genomics_check.py` → `results/predict_genomics_check.json`.

**What it does:**
- The held-out lines of CV fold 0 (208 dev lines) are written in DepMap release layout from `data/`, into the
  git-ignored `cache/predict_genomics_check/`.
- **Block check:** the blocks the script builds must equal the cached blocks the CV run used (C, M, H, metadata,
  and M2's gene-level CN and damaging).
- **End-to-end check:** `predict.py` trains on fold 0's training lines (`--exclude-lines`) and must reproduce the
  registered out-of-fold predictions of the m5_avg CV run, within 1e-3, the float32 storage tolerance.

**Results:**
- Blocks: identical (max difference 0, the same NaN pattern).
- End-to-end: all 208 lines × 18,434 genes reproduced, max |difference| 2.4e-7 and mean 4.4e-9. The CLI took 389 s:
  the fit stores the data, and the three inputs are refitted at prediction time.
- **Why the match is meaningful:** the same fold models applied with expression only differ from these
  predictions (dev-CV P1 0.3135 against 0.3175), so the check exercises the genomic path.
