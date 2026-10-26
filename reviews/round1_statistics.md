# Review round 1 — statistics

## 1. Verdict

**Accept with revisions.** The cancer half is sound. Splits are locked, patient-grouped and lineage-stratified (255 test lines from 247 patients; dev folds of 201–209 lines). Baselines are fair and the main leakage controls are present. There are three gaps:
- (a) Some metric rules are ambiguous, and one of them biases the baselines.
- (b) The 1-SE rule and the ablation claims lack paired inference.
- (c) The non-transformed half has no defined unit, null or estimand. Set A has about 11 independent backgrounds and each external cell type has one screen.

The harmonisation checks as written cannot detect the affine step they are meant to validate.

## 2. Proposed changes

1. **§4 Unit of computation (must).** Compute per-gene r within each held-out fold, then average over folds. Never pool out-of-fold predictions. *Why:* B0 in fold k is the training mean, which falls linearly as the held-out fold's mean rises. Pooled across folds it scores a spurious **negative** r, and any model with fold-specific intercepts inherits the artefact. Add a unit test that B0 scores exactly 0.

2. **§4 Scale-aware complement (should).** Add the per-gene skill 1 − SSE/SSE_B0 and the per-gene r of residuals from B1, on *variable* genes. *Why:* Pearson ignores offset and scale. §7 transfer needs calibrated magnitude, and "beats lineage" should be measured directly.

3. **§4 Gene classes frozen (must).** In CV, define the classes per training fold. For test, set A and external screens, use one list from the 1022 dev lines, committed and hashed as `configs/gene_classes.csv` before unblinding. *Why:* classes built from the evaluated lines would select genes for their chance variance there. That is selection on the outcome, and it inflates absolute levels. The P2b inclusion rule (≥3 positives and ≥3 negatives) uses held-out labels but never predictions, so it is identical across models. Keep it, and report the gene count.

4. **§4 P2b prevalence (must).** Replace mean(AUPRC − π) with (AP − π)/(1 − π), and also report AUROC. The headline P2b uses selective genes with held-out π ≤ 0.5; report genes with π > 0.5 separately. Use `sklearn.metrics.average_precision_score` and unit-test that a constant prediction gives exactly π (B1 is tied within lineage). *Why:* raw differences are bounded by 1 − π. When π is high, the "dependency" is the majority class.

5. **§4 P2a (should).** Report the per-line ΔAP versus B0, paired by line. *Why:* B0 already ranks common essentials, so raw P2a mainly measures the gene mean.

6. **§4 Selection rule (must).** Take the SE from **paired** fold-wise differences against the best configuration: SD(5 Δ)/√5. Pre-declare a complexity order (B0<B1<B2<M1<…<M5; within a rung, fewer modalities is simpler) and a minimum effect δ = 0.005 in P1. "Beats the previous rung" (§5) uses the same rule. *Why:* the unpaired fold SD adds fold difficulty to the noise term. CV SDs also understate the true variance, because training sets overlap (Bengio & Grandvalet 2004).

7. **§4 Test protocol (must).** Name one headline model before unblinding. The other finalists are secondary and are never chosen among on test.

8. **§6 Ablation multiplicity (should).** Treat the ablation as estimation over a pre-listed family of contrasts (add-one, leave-one-out, E-free, prior-knowledge ladder; about 25). Report paired Δ with CIs from a patient bootstrap within each held-out fold. To claim that X adds information, require a Holm-adjusted CI that excludes 0 and |Δ| > δ. Otherwise report "no detectable contribution".

9. **§4 Uncertainty (should).**
   - No gene-level tests that treat about 2k correlated genes as independent.
   - Calibrate "n genes r > 0.3" and "fraction beating lineage" against the permuted-label control.
   - The test bootstrap resamples patients stratified by lineage: 2000 resamples, percentile CIs, paired differences.
   - Label the fold SD as descriptive, not as a CI.

10. **§3/§4 Power (should).** State the following.
    - At n = 255, one per-gene r has SE ≈ 0.06, so the estimand is the mean over genes.
    - The probe gap (0.31 vs 0.18) is easily resolved.
    - A finalist-vs-finalist gap under about 0.01 counts as unresolved if its CI covers 0.
    - Report per-gene metrics by lineage only for lineages with ≥20 test lines (CNS 24, Lung 21). Other lineages get per-line metrics only.

11. **§3/§7.4 Set A (must).**
    - Effective n = 11 patients, because the five RPE1-ss lines come from one patient.
    - The mechanisms are mixed: hTERT (RPE1), SV40 large T (HA1E, BPH-1), virally transformed lymphoid lines (C8166, HG-3), an adenoma organoid and prostate lines. Stratify by mechanism.
    - The primary unit is the line. Rank each line's per-line metrics (P2a and item 12) as a percentile within the cancer-test distribution. Resample by patient.
    - Drop the per-gene correlation across the RPE1 sublines. It measures aneuploidy at n = 5, with per-gene SE ≈ 0.4.
    - Only 10 of the 15 have CN: use E-only models here, or report the reduced n.

12. **§7.4 Define "context-specific dependencies survive" (must).** With one screen per cell type, judge each screen against nulls.
    - (i) Specificity: compare the matched deviation correlation (Spearman of pred − μ_g against obs − μ_g over selective genes) with the correlations from every other line's prediction. Report the rank per screen and a permutation test across screens.
    - (ii) Survival: a predicted dependency (pred < −0.5 and pred − μ_g < −0.25) survives if obs < −0.5. Report precision and recall, each divided by the same metric computed between replicate screens of the same cell type (the ceiling).
    - (iii) Pre-specified genes: z-score obs − pred against the per-gene replicate noise SD from the harmonised anchors, and report counts with exact binomial CIs.
    - Results are descriptive, with no population claim.

13. **§7.2 step 5 — harmonisation validation (must).**
    - Affine scaling preserves ranks, so held-out CEG/NEG AUROC **cannot change** at step 2. Validate step 2 by the held-out CEG median ≈ −1, NEG median ≈ 0, their spreads and SSMD. Use AUROC only for step 3.
    - Repeat the half-splits over 20 seeded splits and report the spread.
    - DepMap scales each line on control lists that overlap CEGv2/NEGv1, so set A's CEG/NEG separation is partly built in. Also report transfer on dev-defined common essentials outside CEGv2.

14. **§7.2 step 3 — anchors (must).**
    - Measure agreement on deviations from the gene mean (selective genes) and as same-line recall rank among DepMap lines. Whole-profile r is about 0.9 from gene means alone.
    - Use empirical-Bayes shrinkage, with leave-one-anchor-out only for libraries with ≥3 anchors.
    - Exclude the evaluated cell type from the anchors. Anchoring external RPE1 on DepMap RPE1, which is in set A, makes the two evaluations dependent.
    - Assert that anchors are dev-only. All currently are, but only by chance.

15. **§8 Leakage additions (should).**
    - Flag test–dev near-duplicates: expression correlation above the 99.9th percentile of between-patient pairs, plus derivative notes in PublicComments. Run a sensitivity analysis without them.
    - Compute the fusion-prevalence filter and the §7.3 quantile reference on dev only.
    - Fill NaNs from the training fold only (the probe used full-dev means, which is acceptable for calibration only).
    - Require ≥20 observed held-out lines per gene and fold.

16. **§11 Chronos joint fit — disclose (should).**
    - 26Q3 Chronos fits guide efficacies and screen parameters jointly over all 1361 lines, test and set A included, and scales each line to its control-gene medians. Test targets are therefore not strictly independent of training lines.
    - The shared information is mainly gene-level: an approximately common gene-wise affine map, with each line contributing about 1/1361.
    - Per-gene across-line r (P1, P2b) is invariant to that map, so leakage there is negligible. RMSE, calibration, P2a and B0 may be slightly favoured.
    - Deployment is equivalent, since a newly screened line uses the same pretrained efficacies.
    - Refitting on dev alone needs the raw reads (could, if ScreenSequenceMap is obtained).
    - Note the asymmetry with the external screens, which are re-fitted on their own data.

## 3. No change needed
The thresholds (−0.5, with −1.0 as a sensitivity check), "variable" defined as SD > 0.25, and r = 0 for constant predictions (given item 1). The split construction and asserts (test lineage proportions deviate slightly, which is acceptable). Seeds, nested CV and the control set.
