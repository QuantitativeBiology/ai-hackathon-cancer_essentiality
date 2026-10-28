# Independent review of G6 (transfer to non-transformed cells): reproduction and red team

Reviewer: an independent agent, 2026-10-28. It reviewed the outputs at c5b2ec8 (stage A) and f348ffb (stage B).
The subagent could not write files in reviews/ (harness policy), so the PI transcribed this report from its
hand-back and copied its scripts and aggregate outputs here unchanged:
- predict_folds.py: per-fold predictions from the persisted fold models, with no refit.
- reproduce_g6.py: all statistics, the bootstrap, the hypotheses, set A, and the comparator recomputation.
- compare.py: comparison against the project outputs. Its outputs are cmp_boot.txt and repro_out_boot.json
  (aggregates only).
- The fold-prediction parquet stayed outside the repo.

The reviewer read the observed data directly, after G6, without calling the guard:
- the harmonised set-B profiles (cache/harmonised/<arm>.parquet, column `final`; for Bowden, `affine`);
- set A from data/CRISPRGeneEffect.csv.

These reads are logged in results/test_access_log.csv as 19 rows with access_type
independent_reproduction_read, timestamped 2026-10-28T03:43:44. They change no model or decision.

## Summary verdict
- **Reproduction:** every G6 primary number reproduces with the reviewer's own code. 288 values were compared:
  per-entry skills, background points, 95% CIs and SEs, ensemble points, absolute APs, swap nulls, ceiling
  points, consensus QC, H1/H2/H3 and set-A per-line values.
  - Skills and CIs match to 3e-16.
  - The largest difference is 4.9e-5, and it comes only from the 4-decimal rounding of the frozen ceiling and QC
    values.
  - No discrepancy exceeds 1e-4.
- **Blinding and leakage:** no leakage, blinding or protocol violation that affects a verdict.
- **Verdicts:** all four background × question verdicts follow correctly from the frozen rules:
  - Q3a: transfers in RPE1 and in hPSC.
  - Q3b: inconclusive in RPE1 and in hPSC.
  - H1 pass, H2 pass, H3 fail, H4 inconclusive by declaration.
- **The catch:** the "transfers" verdicts are weaker than the word suggests (findings M1–M3). They need careful
  wording in REPORT.

## 1. Reproduction

**Method.**
- `essml.transfer` and `evaluate_transfer.py` were not imported.
- Metrics are the reviewer's own: step-wise weighted AP and recall at 5% FPR. Its AP equals sklearn's
  `average_precision_score` to 1e-17, also with duplicated genes as weights.
- The bootstrap is the reviewer's own implementation of the declared scheme (`SeedSequence(20261026).spawn(3)`,
  study draws, then paralog-block draws).
- Per-fold predictions:
  - The saved table holds only fold means (see m1).
  - So the reviewer loaded the persisted fold pickles (no refit) and predicted the declared inputs.
  - The fold means equal both saved prediction tables to 8.9e-16. The table hashes match the summaries.

| Quantity | Reviewer | Project | Diff |
|---|---|---|---|
| RPE1 Q3a skill (T, D7) | 0.0320714 [0.0223049, 0.0400690], SE 0.0046476 | same | < 2e-16 |
| RPE1 Q3b skill | 0.0290901 [−0.0027455, 0.0605502] | same | 0 |
| hPSC Q3a skill | 0.0089362 [0.0039731, 0.0140254] | same | < 2e-16 |
| hPSC Q3b skill | 0.0079594 [−0.0288845, 0.0375274] | same | 0 |
| RPE1 / hPSC Q3a, T ensemble | 0.0327212 / 0.0107971 | same | < 1e-16 |
| Per-entry skills (5 RPE1 + 1 hPSC entries × 5 statistics) | all | all | < 1e-16 |
| RPE1 ceiling Q3a / Q3b (LOSO point) | 0.030724 / 0.001916 | 0.0307 / 0.0019 | 2.4e-5 / 1.6e-5 (rounding) |
| hPSC ceiling Q3a / Q3b | −0.098386 / −0.064904 | −0.0984 / −0.0649 | 1.4e-5 / 3.9e-6 (rounding) |
| Swap RPE1 Q3a, MEF / laminin | −0.01089 [−0.0192, −0.0037] / −0.01615 [−0.0251, −0.0083] | same | < 1e-16 |
| Swap hPSC Q3a, RPE1 expression | −0.00032 [−0.0050, 0.0041] | same | < 1e-16 |
| Consensus SSMD and medians (6 entries) | e.g. Brown WT −2.33477 | −2.3348 | ≤ 4.9e-5 (rounding) |
| Plan genes / blocks | RPE1 16839 / 10645; hPSC 16270 / 10447 | same | 0 |
| Set A Q3a per line | 0.03534, 0.04255, 0.01810, 0.04003, 0.03680 (CIs match); median pct 93.33 | same | < 4e-16 |
| Set A Q3b per line | 0.06283, 0.08314, 0.04241, 0.05835, 0.05921; median pct 95.29 | same | < 2e-17 |
| H1 model Δ MDM2 / PPM1D | −0.62298 / −0.28118 | same | 0 |
| H1 observed Δ MDM2 / PPM1D | −1.61215 / −1.70376 | same | 0 |
| H2 pooled | n 15, ρ 0.95714, sign 0.80 | same | 0 |
| H2 Brown / Bowden | ρ 0.92527 / 0.94545 | same | 0 |
| H3 | 4 observed depleted; T 0/4; no-EC 0/4 | same | 0 (no-EC per-gene ≤ 4e-5, 4-dp rounding) |
| B0 / B1 comparators recomputed from the raw CSV (1022 dev lines; B1 κ = 10) | max abs diff vs comparators.parquet | | 2.7e-8 |
| `frac_below` recomputed | 0 genes cross the 0.5 threshold | | ≤ 5e-6 |

**Verdicts re-derived from the reported CIs, frozen ceilings and MDEs, harmonisation pass and anchor gate:**
- **RPE1 Q3a: transfers.**
  - Lower CI 0.0223 > 0, and 0.0321 ≥ 0.5 × 0.0307.
  - Swap lower CIs are below 0, so no cap.
- **hPSC Q3a: transfers.** Lower CI 0.0040 > 0; the ratio condition is vacuous (ceiling −0.098).
- **RPE1 Q3b and hPSC Q3b: inconclusive.** Lower CI < 0, and "fails" is blocked by the anchor gate.
- **Overall:** Q3a transfers, Q3b inconclusive. All match the project.

## 2. Verified clean
- **Freeze:**
  - Tag `freeze-transfer` is commit 987bdec (2026-10-26 20:21:45 UTC).
  - `git diff freeze-transfer HEAD` touches none of the frozen configs, `transfer.py` or `harmonise.py`. The only
    change is `scripts/transfer_criteria.py`, the accepted D21c change.
  - Every `inputs_sha256` entry matches, except that script, which matches its pinned D21c content.
- **Order of events:**
  - Set-B harmonised files were last written at 19:00:47 on 10-26, before the freeze.
  - `finalists.json` was committed at 84efb5c (00:00 on 10-28).
  - T fold models were persisted at 00:36. The no-EC, m5_avg and B2 models were persisted at 01:55–01:59.
  - Set-A access was at 02:37–02:39.
  - Commit c5b2ec8 (02:41:28) holds the sha256 of a prediction table. That table contains T's WT predictions for
    every set-B expression input (RPE1 state, both H1 states), and these are identical in stage B.
  - Set-B access was at 03:22–03:34.
- **Predictions fixed before access:** the WT-flag set-B predictions were committed by hash before any set-B
  access. The null-flag predictions come deterministically from the persisted models, and the reviewer
  regenerated them exactly.
- **One access per finalist × screen:** 55 (finalist, scope) target pairs, each granted once, with no refused
  rows. After the freeze there are no set-B reads other than the G6 accesses and the dev-line anchor profiles
  (D22f).
- **No leakage into training:**
  - All 5 T and 5 no-EC fold training sets are dev-only.
  - No-EC excludes the 4 EC lines (ACH-001578, -002291, -001374, -001668); T includes them.
  - The RPE1 patient appears only as non-cancerous lines; no dev line shares it.
  - There are no hESC lines in DepMap (the DepMap line named "H9" is a lymphoid cancer line).
  - Anchors are dev lines scored by their single out-of-fold model: HAP1 and DLD1 in fold 2, HCT116 in fold 3,
    A375 in fold 4.
- **Expression inputs are the declared ones:**
  - `state:DepMap_RPE1` is the TPM-mean of the 6 DepMap RPE1-ss rows (max diff 9.5e-7).
  - hPSC uses the per-arm H1 states, MEF and laminin.
  - Each entry gets its own TP53 flag.
- **Harmonisation followed the pre-registration:**
  - Gene scores are median guide LFC; Chronos was rejected by the pre-set D11 rule.
  - The affine step was applied. `final` equals `affine` for all 14 verdict arms.
  - The gene-wise correction was evaluated by leave-one-anchor-out (LOAO) and not applied for TKOv1/v2/v3: LOAO gain
    0.0026 / 0.0007 / −0.0029, all below the 0.01 margin. It is not identifiable for Brunello.
  - The consensus QC reproduces the frozen values.
- **Skill is not trivial:**
  - The comparator is the gene mean B0. The scaling set (DepMap common essentials) is already excluded from U.
  - Excluding all DepMap common essentials, T's per-entry Q3a skill is RPE1 0.023–0.056 and hPSC 0.019.
  - Restricted to expressed genes (log2 TPM+1 ≥ 1), it is 0.024–0.049 and 0.011.
  - A trivial rule (B0, but 0 for unexpressed genes) scores −0.003 to +0.006 in RPE1 and −0.019 in hPSC.
  - Swapping in another background's expression removes the skill: matched minus swap is +0.043 / +0.048 in RPE1
    and +0.009 in hPSC.
  - For context, kNN on expression (B2) gets 0.023 in RPE1, so part of the signal is generic expression
    similarity.

## 3. Findings

### Critical
None.

### Major

**M1. The verdict pattern is decided by the uninformative anchor gate.** No anchor has positive skill, which the
pre-access dry run forecast (D21e).
- In all four cells, the "fails" statistical condition holds (upper CI < MDE), and harmonisation passed:
  - RPE1 Q3a: 0.0401 < 0.101
  - hPSC Q3a: 0.0140 < 0.082
  - RPE1 Q3b: 0.0606 < 0.1617
  - hPSC Q3b: 0.0375 < 0.0816
- Only the gate stops "fails": 0 of 4 library-matched anchors have lower CI > 0. HAP1 is −0.0022 [−0.0067, 0.0021]
  for Q3a and +0.0285 [−0.041, 0.058] for Q3b.
- If any anchor had been positive:
  - Q3a would be "inconclusive (transfers and fails conditions both met)" in both backgrounds;
  - Q3b would be "fails" in both.
- So "transfers" here means: skill is above 0, but the whole CI lies below the pre-registered MDE.
- The design was effectively one-sided: "fails" was unreachable, and this was known before access (D21e/f).
- Nuance: with the fixed-pool MDE sensitivity, RPE1 Q3a would no longer meet the fails condition (0.040 > 0.020);
  hPSC Q3a would (0.014 < 0.022). Evidence: results/transfer/setB/{RPE1,hPSC}.json → anchor_gate;
  transfer.py:442-457.

**M2. "Transfers" effectively means "lower CI > 0 against the gene mean", and the effect is small.**
- The ratio condition is near-vacuous in RPE1 (0.5 × 0.0307 = 0.0154) and vacuous in hPSC.
- The PI declined the reviewer's absolute floor (transfer_freeze_response.md, F5 and "Points for the PI" #1) after
  the ceiling values were known, though still blind to T.
- Effect size:
  - RPE1 per-entry AP of T is 0.43–0.80, against 0.42–0.76 for B0 (about +0.03, roughly 5% relative).
  - hPSC is 0.414 against 0.405.
- Robustness:
  - Under DESIGN §7.6's original rule (≥ 50% of the same-pipeline anchor skill), T still clears, because anchor
    Q3a skills are ≤ 0.0033.
  - RPE1 skill is 1.04× the replicate ceiling. Mostly this says the cross-study ceiling is low.
  - Set A agrees in magnitude: 0.018–0.043 on DepMap's own RPE1-ss screens, against a median of 0.011 for cancer
    test lines.

**M3. hPSC "transfers" rests on one study, with a CI that understates the uncertainty.**
- Mair 2019 is one lab with 2 arms. Its CI carries gene-block variance only; resampling one study is the identity.
- In RPE1, study-level variance dominates: ceiling SE 0.051 against 0.010 with fixed pools.
- The hPSC ceiling is negative (−0.098): the independent hPSC screens (Ihry, Yilmaz) predict Mair worse than the
  gene mean, a sign of large study or library effects in this background.
- The swap null supports expression specificity: RPE1 expression gives −0.0003.
- But the claim must be limited to "the one hPSC screen evaluated".

### Minor
- **m1. Primary statistic not reproducible from saved artefacts.** `predictions_with_confidence_*.parquet` stores
  only fold means (evaluate_transfer.py:243-248), while the primary D7 statistic averages over folds. The reviewer
  had to regenerate per-fold predictions from the cache/models pickles. Recommendation: save per-fold predictions,
  or pin the pickle hashes.
- **m2. Non-control set-B reads before the freeze.**
  - 18 `replicate_ceiling_aggregate` log rows (10-26, 16:20–20:16) read full set-B profiles. That is 5 ceiling
    iterations before D13 and 3 D13 ceiling runs.
  - These are required by DESIGN §7.2.5/§7.6, logged, aggregate-only, and involve no model. But they are not in the
    §9 exception list, and the D13 criteria were finalised after seeing them.
  - Add to the same REPORT disclosure the Chronos joint fit (D11) and the unlogged column-name read (D21d).
- **m3. Harmonised set-B files are not hash-pinned.** cache/harmonised/*.parquet is not in `inputs_sha256`.
  Integrity currently rests on file mtimes (before the freeze) and on the reviewer's independent match of the
  frozen consensus QC and ceiling points (≤ 5e-5). Recommendation: pin them.
- **m4. Named tags do not exist.** The criteria `_status` names tags `freeze-transfer-criteria` and
  `freeze-transfer-models`, which were never created. The requirement is met in substance by `freeze-transfer`
  (987bdec) and `freeze-finalists` (84efb5c), which predate access. Note this in the record.
- **m5. Undocumented change to the fails rule.** D13 dropped DESIGN §7.6's "fails requires skill ≤ B0/B1" and added
  "both conditions met → inconclusive". This change is in the criteria file but not in the RESEARCH_TREE D13
  bullets. Under DESIGN's original wording, a gate-positive counterfactual would still give Q3a "transfers".
  Record it as a deviation.
- **m6. POU5F1 is missing from Mair.** It is absent from both Mair harmonised arms (the `raw` column), so H3's
  canonical gene is unobservable. H3 was therefore tested on 4 genes (NANOG, SOX2, LIN28A, PRDM14), with T 0/4
  and no-EC 0/4. The "fail" stands, and was predicted before access (D20e). Report the missing gene; the cause was
  not determined.
- **m7. Q3b positive sets are inflated by noisy screens.**
  - Mair has 3744 of 16270 (23%) positives and Haapaniemi 3571 of 16767 (21%), against 266 for Olivieri.
  - Their NEG evaluation-half medians are −0.076 and −0.105.
  - Q3b in these entries largely scores noise below −0.5. This is part of why Q3b is inconclusive.
- **m8. Set A and set-B RPE1 are not independent on the prediction side.**
  - Set-B RPE1's expression input is the mean of the set-A lines plus ss111.
  - Set-A targets come from DepMap's joint Chronos fit with dev lines (the D26-5 limitation also applies to set A).
- **m9. H1/H2 test the TP53 flag.** Both are pure flag flips at fixed expression. They pass, but they show the TP53
  feature reproduces isogenic p53 effects (H2 ρ 0.96). They are not evidence of expression-based transfer. The
  project already labels them this way; keep that wording.

## 4. Can the G6 results be relied on, and how to word them
- **The numbers:** yes. They are correct, and the protocol was followed: one blinded access per finalist × screen,
  no leakage, and the declared inputs.
- **Q3a:** the verdicts are correct under the frozen rules. Suggested REPORT wording:
  > "In RPE1, T ranks the 10% most depleted genes better than the DepMap gene mean in each of 4 studies
  > (ΔAP 0.032, 95% CI 0.022–0.040; about 5% relative). This matches cross-study RPE1 agreement (ceiling 0.031)
  > and is specific to the RPE1 expression input. In the single hPSC screen (Mair 2019) the gain is 0.009
  > (0.004–0.014; gene-block CI only, no study replication; independent hPSC screens agree with Mair less well
  > than the gene mean does). By the pre-registered rule this is 'transfers'. But the effect is small: its whole
  > CI lies below the pre-registered MDE. The verdict is 'transfers' rather than 'inconclusive' only because no
  > transformed anchor screen showed positive skill, so 'fails' could not be reached (forecast before access,
  > D21e)."
- **Q3b:** report as
  > "Not detectable with these screens. Context-specific gains are below the MDE, and the anchors cannot certify
  > the pipeline's sensitivity."

  Do not call this "fails", and do not call it a "trend toward transfer".
- **H3:** report "fail", with POU5F1 unobservable.
- **H1/H2:** report as TP53-feature checks.
