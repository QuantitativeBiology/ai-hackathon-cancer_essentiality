# Round 2 review: data (DESIGN v1)

## Verdict
Accepted with minor changes. All round-1 data items are addressed in v1. The remaining issues are small and concern the evaluation gene set, the CN fill value and the instability threshold.

## Checks of v1 (numbers computed from `configs/splits.csv` and the cache)

**Test and CV gene drops.**
- Test set (255 lines): the rule "at least 50 observed held-out lines" drops 519 of 18,434 genes (2.8%). These are 41 of 2,159 variable genes and 19 of 2,594 selective genes. None of the 817 dev common-essential genes is dropped.
- The dropped genes have a mean dev effect of +0.02, against −0.15 for kept genes, and are 82% NaN in dev. They are genes screened in only one library group.
- CV folds (rule: at least 40 observed): 504, 519, 519, 504 and 179 genes dropped in folds 0 to 4. The variable genes dropped are 34, 41, 41, 34 and 6.
- The evaluable gene set therefore differs by fold, and fold 4 is very different from the others.
- No data problem arises for fitting. Every gene has at least 135 observed training rows in every fold. The smallest of the 44 dev NaN patterns has 174 observed rows, so the per-pattern exact ridge is well-posed.

**P2b on the test set.** 27 of 2,594 selective genes have fewer than 3 positives or fewer than 3 negatives among the 255 test lines. The existing ≥3/≥3 rule already handles them.

**CN transform.** log2(CN+0.1) has quantiles of −3.32 (min), −1.61 (0.1%), 0.17 (median), 1.96 (99.9%) and 8.03 (max). Only 0.11% of values are clipped at ±3, which is acceptable.

## Remaining changes

**1. [MUST] §4: fix one evaluable gene set per fold type.**
- Define the evaluable set once before modelling: genes with at least 40 observed lines in every fold and at least 50 in test.
- Commit it as `configs/eval_genes.csv`, alongside `gene_classes.csv`.
- This removes fold-dependent gene sets from the paired comparisons. Fold 4 has only 179 dropped genes against about 510 in the others.
- Report the dropped 519 genes as "not evaluable" and count them in the summary. Their effect is near zero, so this is a coverage statement, not a bias.
- Without this, P1 differences between folds partly reflect which genes are scored.

**2. [MUST] §2 C: correct the fill value.**
- v1 says other NaNs are filled with 0 and calls that "diploid". In log2(CN+0.1) space diploid is log2(1.1) = 0.14, and the genome-wide median is 0.17. Zero corresponds to CN of about 0.9.
- Change: fill NaNs with the training-fold gene median in that space, or with log2(1.1). Standardisation then centres it.

**3. [MUST] §2 C: re-derive the instability threshold on the right scale.**
- The score is "fraction of genes with |log2CN| > 0.3". With the new offset, the median line already has 41% of genes above 0.3 (0.14 of the diploid baseline sits inside that), so the score would mostly measure the offset.
- Change: compute it on raw relative CN as the fraction of autosomal genes with |CN − 1| > 0.3. Its median per line is 0.23, which is a usable spread. Define the threshold in the config.

**4. [SHOULD] §2 and §4: CV-fold gene classes should keep gene-class counts stable.**
- Rare-selective genes are defined as at least 3 and fewer than 5% of lines below −0.5. For NaN-heavy genes, "lines" means observed lines only. State this in the config, otherwise the prevalence for those genes is computed on a different denominator than for fully observed genes.

**5. [COULD] §2 Chronos disclosure.**
- Add one line that Chronos NaNs (not filled) arise from the library's gene set, so "NaN group" is a library proxy but not a true screen label. This is already requested via ScreenMap and needs no further action.

## No change needed
- Exact ridge per NaN pattern: feasible and well-posed (44 dev patterns, at least 174 observed rows each).
- Fusion filter and F encoding: as specified in v1 (self-fusions, non-coding partners, at least 5 training lines). Nothing further needed.
- Meta restrictions and set-A tiers: accepted.
