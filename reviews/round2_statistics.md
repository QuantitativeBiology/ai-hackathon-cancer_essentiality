# Review round 2 — statistics (DESIGN v1)

## Verdict
**Accept after three must-fixes.**
- Round-1 items 1, 3, 4, 6, 7, 11, 13 and 14 are adequately addressed.
- Items 2 and 9 are partly addressed.
- Item 12's new headline does not measure context specificity.
- §7.6 never defines "transfers".

## Must
1. **§7.5 recall at fixed NEGv1 FPR.**
   - *Problem:* "not common-essential" admits genes that are dependencies in up to 95% of cancer lines. B0 ranks those from the gene mean alone.
   - *Change:*
     - Positives = observed < −0.5 **and** dev dependency fraction ≤ 0.5.
     - Headline = **Δrecall vs B0 and vs B1**.
     - FPR is computed on the **evaluation half** of NEGv1.
     - Ceiling = the same statistic with a same-cell-type replicate screen as the predictor.
     - Clopper–Pearson CIs.
2. **§7.6 decision rule.**
   - *Problem:* only "fails" is defined, with no uncertainty method or aggregation rule.
   - *Change:* commit the following.
     - **Transfers:** lower 95% bound of skill over B0/B1 > 0 **and** point estimate ≥ 50% of same-pipeline anchor skill.
     - **Fails:** the current conditions plus an upper bound below the MDE.
     - **Inconclusive:** everything else.
     - Uncertainty: gene block bootstrap, with blocks = complexes/paralog families.
     - Category verdict by a count rule, e.g. "≥ ⌈n/2⌉ screens transfer and none fails".
     - Primary verdicts are only Q3a and Q3b. H1–H4 are secondary counts.
3. **§3 A1 is one genetic background.**
   - *Change:* state that A1 is a single-background case study. Report per subline plus a background-level mean. Claims about hTERT cells in general also need set B.

## Should
4. **Holm (§4/§6).** Holm adjusts p-values, not CIs.
   - Apply Holm to paired patient-bootstrap p-values and show Bonferroni CIs (1 − 0.05/m).
   - Commit the contrast family to `configs/ablation_contrasts.json` before running.
   - Restore |Δ| > δ.
   - Note that within-fold bootstrap CIs hold the models fixed and are anti-conservative. Show the fold-wise t(4) interval too.
5. **Power statement (§4).** "Averages over ≥2,000 genes" implies 0.06/√2000. That is wrong, because gene errors are correlated through shared lines.
   - Replace it with an empirical MDE: patient-bootstrap ΔP1 (best vs B1) on dev out-of-fold predictions subsampled to 255 lines, MDE ≈ 2.8 × SD.
   - Finalist gaps below the MDE count as unresolved.
6. **Rare-selective (§4).** "≥3 lines < −0.5" admits many noise genes (SD ≈ 0.14, heavy tails).
   - Also require min < −1.0 and ≥3 lines below median − 4·MAD; exclude NEGv1.
   - Report P1 for this class only for genes with ≥3 held-out positives. Use per-gene AUROC as the main metric.
   - Use the same list in §7.5.
7. **Round-1 leftovers.** Add the skill 1 − SSE/SSE_B0. Calibrate "n genes r > 0.3" and "fraction beating B1" against the permuted-feature control.
8. **Profile-swap null (§7.5).**
   - Exclude same-cell-type profiles from the mismatches and add cancer-test profiles.
   - Pre-declare the statistic: rank of matched among mismatched, combined by a permutation test on summed ranks.
   - The expression-swap control reuses it.
9. **Set A ranking (§3).** Build the cancer-test reference distribution with the same model applied to set A (T).
10. **Selection wording (§4).** ΔP1 = mean of the 5 fold differences, SE = SD/√5. Add B2 to the order. The R rungs are comparators, not candidates.

## Could
11. Near-duplicates (§8): use the 99.9th percentile of between-patient correlation as the threshold, plus an exclusion sensitivity analysis.
12. Anchor skill comes from CV models trained on ~818 lines, while non-transformed cells are predicted by the frozen 1022-line model. Note the slight asymmetry.

## No change needed
- P1 fold aggregation and the B0 test.
- Frozen gene classes.
- AP normalisation, ties and the π ≤ 0.5 headline.
- Within-lineage rule.
- Observed-line minimum.
- Affine validation.
- Anchor restriction.
- Chronos disclosure.
- Test guard.
