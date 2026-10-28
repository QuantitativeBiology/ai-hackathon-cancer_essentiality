# m5_avg applied with expression only: post-hoc dev-CV diagnostic (G7)

**Status:** post hoc. Added at G7 after the REPORT red team. It uses dev lines only; no test-set or external
target was read. The rule for the `scripts/predict.py` default (below) was fixed by the PI before the number was
seen.

- **Code:** `scripts/predict_expression_only_cv.py`.
- **Output:** `results/predict_expression_only_cv.json`, which holds every number below, the fold-model hashes
  and the commit.

## Why
`scripts/predict.py` reads expression only. With `--model m5_avg` (the default), M2 and M3's genomic inputs (copy
number, damaging and hotspot mutations, fusions) are therefore absent, which means diploid wild type, and the
metadata (lineage, sex, age) is unknown. This is `StackM5.predict_from_blocks`, the same path G6 used for external
cells. The G4 numbers for m5_avg (test P1 0.319) used the real genomic inputs, so the configuration users actually
run had not been evaluated on cancer lines.

## Method
1. **Fold models.** The persisted G6 fold models `cache/models/m5_avg_6dbebcf4` are used: fold k is trained on
   the other four dev folds. They were checked against CV run `20261027-161041_m5_avg_6dbebcf4`:
   - the config hash (`6dbebcf4…`) matches;
   - the fold signature matches `configs/splits.csv` cv_fold;
   - each fold's training lines are exactly the other folds;
   - the pickle sha256 values equal those in `results/transfer/g6_input_hashes.json`.

   m5_avg is in "average" mode, so it stores only its training data. Each fold model refits its three inputs
   (m1_E_resB0, m2n_E_paralogs, g3_m3n) on that fold's training lines, as in CV, G6 and predict.py
   (`StackM5.fitted_inputs`). No config, seed, weight or input was changed.
2. **Self-check.** Each held-out line was predicted with all of its own raw blocks: E, C, M, H, meta, and the
   gene-level expression, CN and damaging data, as in the G6 self-check. This was done for every held-out line,
   not just a sample. The result reproduces the registered out-of-fold predictions
   (`cache/runs/20261027-161041_m5_avg_6dbebcf4/pred.parquet`) with a maximum |diff| of 2.4e-7 in each of the 5
   folds (tolerance 1e-3).
3. **Expression only.** Each fold's held-out lines were predicted exactly as predict.py predicts a new line:
   - IDs `input::<line>`;
   - `X = {"E": the line's DepMap log2(TPM+1) row}`, which equals `cache/features/E.parquet` (`--harmonise none`);
   - lineage `None`.
4. **Scoring.** The frozen `essml.evaluate.evaluate` was called as `scripts/run.py` calls it: dev-wide
   `configs/gene_classes.csv`, `configs/eval_genes.csv` and min_obs 40. Applied to the registered m5_avg OOF, the
   same call reproduces its `per_fold.csv` P1 within 5e-7.
   - **Paired fold differences:** taken against the registered `per_fold.csv` of m1_E_resB0
     (`20261026-222117_m1_E_resB0_91e52f9d`) and of m5_avg, with `scripts/g2_report.paired`
     (SE = sd / √5).
   - **Hotspot AUPRC:** as in docs/mutation_dependencies.md, using `essml.mutdeps.fold_metrics` (all-line AP for
     effect < −0.5) with the declared metric_params, averaged over the 5 folds. The recomputed values for the
     registered runs equal `results/mutation_dependencies.json`.

## Results (dev CV, 5 folds)

### P1 on variable genes

| | fold 0 | fold 1 | fold 2 | fold 3 | fold 4 | mean |
|---|---|---|---|---|---|---|
| m5_avg, expression only | 0.3186 | 0.3145 | 0.3097 | 0.3170 | 0.3075 | **0.3135** |
| m1_E_resB0 (registered) | | | | | | 0.3081 |
| m5_avg with genomics (registered) | | | | | | 0.3175 |

Paired fold differences of the expression-only m5_avg:
- **vs m1_E_resB0:** +0.0054, SE 0.0001. Each fold is positive: +0.0056, +0.0056, +0.0053, +0.0050, +0.0053.
- **vs m5_avg with genomics:** −0.0041, SE 0.0001. All 5 folds are negative.

Removing genomics costs m5_avg about 0.004 of P1 on variable genes. That is about 40% of its CV gain over
m1_E_resB0 (0.0094).

### Other classes

| class | P1 expression-only | m1_E_resB0 | m5_avg with genomics |
|---|---|---|---|
| all | 0.1553 | 0.1547 | 0.1558 |
| selective | 0.2918 | 0.2881 | 0.2942 |
| common essential | 0.2527 | 0.2510 | 0.2557 |

On variable genes, the expression-only m5_avg has P2a 0.969 and P2b nAP (prevalence ≤ 0.5) 0.299.

### Hotspot dependencies (AUPRC, effect < −0.5, mean of 5 folds)

| dependency | m5_avg, expression only | m1_E_resB0 | m5_avg with genomics |
|---|---|---|---|
| NRAS | 0.219 | 0.212 | 0.631 |
| KRAS | 0.767 | 0.756 | 0.805 |
| BRAF | 0.554 | 0.558 | 0.661 |
| PIK3CA | 0.603 | 0.597 | 0.647 |

Without genomics, m5_avg loses the recovery of RAS/RAF/PIK3CA hotspot addiction. The NRAS difference vs m5_avg
with genomics is −0.41 (SE 0.05), which puts it at m1_E_resB0's level. These dependencies are underpredicted by
both expression-only configurations.

## predict.py default
- **Rule, fixed by the PI before the number was seen (DESIGN §4):** the expression-only m5_avg stays the default
  only if its dev-CV P1 on variable genes beats m1_E_resB0 by more than δ = 0.005 and by more than 1 SE of the
  paired fold differences.
- **Result:** +0.0054 > 0.005 and > SE 0.0001, so the rule is **met** and m5_avg stays the default. The margin
  over δ is small: 0.0004, though all 5 folds clear 0.005.
- **predict.py changes:**
  - its m5_avg warning now says that this use was evaluated only post hoc in dev CV, and gives the numbers;
  - the m1_E_resB0 warning and note no longer advise `--model m5_avg` for mutation-driven dependencies, because
    the script takes no genomics;
  - both models flag NRAS/KRAS/BRAF/PIK3CA as underpredicted by both expression-only configurations.
- **Not evaluated:** these are dev-CV numbers only. The test set was not used for the expression-only
  configuration.
