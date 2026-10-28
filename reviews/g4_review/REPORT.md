# Independent review of G4 (one-time cancer test-set evaluation)

Reviewer: an independent agent, 2026-10-28. It reviewed the run at code commit 84efb5c (the phase-A commit) and
outputs at 07234b8/14bfe4c. The code in src/, scripts/ and configs/ is identical between 84efb5c and HEAD.
The subagent could not write files itself (harness policy), so the PI transcribed this file from its report.
The scripts, outputs and access-log row are the agent's own.

The reviewer wrote its own metric code (pandas, numpy and scipy only; nothing imported from essml). It read the
test targets directly from data/CRISPRGeneEffect.csv. This is a post-G4 reproduction read and changes no model
or decision. Before the read, one row was appended to results/test_access_log.csv
(access_type independent_reproduction_read, 2026-10-28T01:37:06). The reviewer refit no model, made no guard call,
and ran neither evaluate_test.py nor evaluate_transfer.py.

Files (aggregate numbers only):
- reproduce.py writes reproduce_output.json.
- compare.py writes comparison.csv (298 comparisons) and prints every difference above 1e-4.

## Verdict

The G4 numbers can be relied on.
- Every headline, per-class, per-tier and Spearman number reproduces within 2.1e-8 (tolerance 1e-4).
- No leakage, blinding or protocol violation was found. All findings are minor.

| | m5_avg (headline) | m1_E_resB0 |
|---|---:|---:|
| P1 variable (test) | 0.3194 [0.3033, 0.3341] | 0.3101 [0.2938, 0.3250] |
| dP1 vs B1 | +0.1532 [0.1376, 0.1700] | +0.1439 [0.1284, 0.1607] |
| dP1 vs B2 | +0.0703 [0.0623, 0.0783] | +0.0610 [0.0534, 0.0688] |

CIs come from the reviewer's own patient bootstrap (lineage strata, seed 777, 2000 draws). The project's CI for
m5_avg is [0.3042, 0.3334].

Paired diagnostic: m5_avg − m1_E_resB0 on variable-gene P1 = +0.0093 [+0.0085, +0.0100], SE 0.00037.
- 0 of 2000 resamples are ≤ 0.
- m5_avg is better on 69.7% of the 2118 variable genes.
- Dev CV gave +0.0094 (D23), so the gain replicates, and it remains small (< 0.01).

## Task 1: reproduction

Checks before computing metrics:
- All 26 prediction files match the committed manifest (bf3faa3), and the manifest is unchanged since that commit.
- The raw CSV sha256 equals configs/depmap_local_sha256.txt. Its shape is 1361 × 18434, and every test line and
  eval gene is present.
- All 17,915 eval genes have ≥ 50 observed test lines.
- n_obs is identical to the project's per_gene files.
- Casting the truth to float32 (the guard serves float32) flips 0 labels at −0.5.
- Per-gene r, nAP and AUROC agree within 5e-6, the 5-decimal rounding of per_gene.csv.gz.

Definitions (DESIGN §4; evaluate.py consulted only for min_obs, the label threshold and the gene sets):
- Lines: the 255 test lines. Genes: eval_genes, with ≥ 50 observed lines.
- Gene classes: gene_classes.csv.
- Pearson r: an exact-constant prediction or truth gives r = 0.
- P2a: step-wise AP per line over the scored "all" genes (label < −0.5, score −prediction), averaged over lines.
- P2b: selective genes with ≥ 3 positives and ≥ 3 negatives; nAP = (AP − π)/(1 − π); AUROC by Mann–Whitney;
  genes with prevalence ≤ 0.5.
- Tiers: predictable_genes.csv (primary) and predictable_genes_m1.csv (secondary, m1 only).

Reproduction table (no discrepancy above 1e-4; the maximum over all 298 comparisons is 2.1e-8;
full list in comparison.csv):

| finalist | metric | reproduced | project | difference |
|---|---|---:|---:|---:|
| m5_avg | P1_all | 0.162246 | 0.162246 | 5.6e-12 |
| m5_avg | P1_variable | 0.319369 | 0.319369 | 5.4e-11 |
| m5_avg | P1_selective | 0.296465 | 0.296465 | -1.5e-11 |
| m5_avg | P1_rare_selective | 0.175513 | 0.175513 | 2.2e-11 |
| m5_avg | dP1_vs_B1_variable | 0.153198 | 0.153198 | 1.1e-11 |
| m5_avg | dP1_vs_B2_variable | 0.070269 | 0.070269 | 2.0e-10 |
| m5_avg | P2a | 0.914959 | 0.914959 | 0.0e+00 |
| m5_avg | dP2a_vs_B0 | 0.010175 | 0.010175 | 3.5e-18 |
| m5_avg | P2b_nap_prev05 | 0.184214 | 0.184214 | 0.0e+00 |
| m5_avg | P2b_auroc_prev05 | 0.662558 | 0.662558 | 0.0e+00 |
| m5_avg | P2b_n_prev05 | 1753 | 1753 | 0.0e+00 |
| m5_avg | primary:high:P1 | 0.383337 | 0.383337 | -7.4e-11 |
| m5_avg | primary:high:P2b_nap_prev05 | 0.285158 | 0.285158 | -5.6e-17 |
| m5_avg | primary:medium:P1 | 0.195733 | 0.195733 | 6.8e-12 |
| m5_avg | primary:medium:P2b_nap_prev05 | 0.133364 | 0.133364 | 0.0e+00 |
| m5_avg | primary:low:P1 | 0.105888 | 0.105888 | 1.8e-11 |
| m5_avg | primary:low:P2b_nap_prev05 | 0.078643 | 0.078643 | 0.0e+00 |
| m5_avg | primary:spearman_evaluable | 0.647530 | 0.647530 | 0.0e+00 |
| m1_E_resB0 | P1_all | 0.161150 | 0.161150 | 5.1e-12 |
| m1_E_resB0 | P1_variable | 0.310098 | 0.310098 | 5.0e-11 |
| m1_E_resB0 | P1_selective | 0.290641 | 0.290641 | -1.4e-11 |
| m1_E_resB0 | P1_rare_selective | 0.175317 | 0.175317 | 2.0e-11 |
| m1_E_resB0 | dP1_vs_B1_variable | 0.143927 | 0.143927 | 6.8e-12 |
| m1_E_resB0 | dP1_vs_B2_variable | 0.060998 | 0.060998 | 2.0e-10 |
| m1_E_resB0 | P2a | 0.914266 | 0.914266 | -1.1e-16 |
| m1_E_resB0 | dP2a_vs_B0 | 0.009482 | 0.009482 | -6.9e-18 |
| m1_E_resB0 | P2b_nap_prev05 | 0.179476 | 0.179476 | 0.0e+00 |
| m1_E_resB0 | P2b_auroc_prev05 | 0.660178 | 0.660178 | 0.0e+00 |
| m1_E_resB0 | P2b_n_prev05 | 1753 | 1753 | 0.0e+00 |
| m1_E_resB0 | primary:high:P1 | 0.374612 | 0.374612 | -7.0e-11 |
| m1_E_resB0 | primary:high:P2b_nap_prev05 | 0.274023 | 0.274023 | -5.6e-17 |
| m1_E_resB0 | primary:medium:P1 | 0.194963 | 0.194963 | 3.1e-12 |
| m1_E_resB0 | primary:medium:P2b_nap_prev05 | 0.132438 | 0.132438 | 0.0e+00 |
| m1_E_resB0 | primary:low:P1 | 0.105826 | 0.105826 | 1.9e-11 |
| m1_E_resB0 | primary:low:P2b_nap_prev05 | 0.078329 | 0.078329 | 0.0e+00 |
| m1_E_resB0 | primary:spearman_evaluable | 0.645289 | 0.645289 | 1.9e-08 |
| m1_E_resB0 | m1_own:high:P1 | 0.376313 | 0.376313 | -4.1e-11 |
| m1_E_resB0 | m1_own:medium:P1 | 0.195474 | 0.195474 | -9.0e-12 |
| m1_E_resB0 | m1_own:low:P1 | 0.105953 | 0.105953 | 2.1e-11 |
| m1_E_resB0 | m1_own:spearman_evaluable | 0.644824 | 0.644824 | 2.1e-08 |
| t_E_status | P1_all | 0.160374 | 0.160374 | 9.0e-12 |
| t_E_status | P1_variable | 0.310773 | 0.310773 | 6.3e-11 |
| t_E_status | P1_selective | 0.291235 | 0.291235 | -2.1e-12 |
| t_E_status | P1_rare_selective | 0.173893 | 0.173893 | 2.5e-11 |
| t_E_status | dP1_vs_B1_variable | 0.144602 | 0.144602 | 1.9e-11 |
| t_E_status | dP1_vs_B2_variable | 0.061673 | 0.061673 | 2.1e-10 |
| t_E_status | P2a | 0.914181 | 0.914181 | 1.1e-16 |
| t_E_status | dP2a_vs_B0 | 0.009397 | 0.009397 | 1.6e-17 |
| t_E_status | P2b_nap_prev05 | 0.180411 | 0.180411 | 2.8e-17 |
| t_E_status | P2b_auroc_prev05 | 0.660502 | 0.660502 | 0.0e+00 |
| t_E_status | P2b_n_prev05 | 1753 | 1753 | 0.0e+00 |
| t_E_status | primary:high:P1 | 0.374898 | 0.374898 | -6.1e-11 |
| t_E_status | primary:high:P2b_nap_prev05 | 0.274968 | 0.274968 | 0.0e+00 |
| t_E_status | primary:medium:P1 | 0.194120 | 0.194120 | 3.9e-12 |
| t_E_status | primary:medium:P2b_nap_prev05 | 0.133367 | 0.133367 | 0.0e+00 |
| t_E_status | primary:low:P1 | 0.104919 | 0.104919 | 2.3e-11 |
| t_E_status | primary:low:P2b_nap_prev05 | 0.079251 | 0.079251 | 0.0e+00 |
| t_E_status | primary:spearman_evaluable | 0.645590 | 0.645590 | 3.0e-09 |
| t_E_status_noEC | P1_all | 0.160664 | 0.160664 | 9.8e-12 |
| t_E_status_noEC | P1_variable | 0.311015 | 0.311015 | 6.6e-11 |
| t_E_status_noEC | P1_selective | 0.291573 | 0.291573 | -3.5e-12 |
| t_E_status_noEC | P1_rare_selective | 0.174104 | 0.174104 | 2.6e-11 |
| t_E_status_noEC | dP1_vs_B1_variable | 0.144844 | 0.144844 | 2.2e-11 |
| t_E_status_noEC | dP1_vs_B2_variable | 0.061915 | 0.061915 | 2.1e-10 |
| t_E_status_noEC | P2a | 0.914206 | 0.914206 | 0.0e+00 |
| t_E_status_noEC | dP2a_vs_B0 | 0.009422 | 0.009422 | 3.5e-18 |
| t_E_status_noEC | P2b_nap_prev05 | 0.180381 | 0.180381 | 0.0e+00 |
| t_E_status_noEC | P2b_auroc_prev05 | 0.660642 | 0.660642 | 0.0e+00 |
| t_E_status_noEC | P2b_n_prev05 | 1753 | 1753 | 0.0e+00 |
| t_E_status_noEC | primary:high:P1 | 0.375112 | 0.375112 | -6.0e-11 |
| t_E_status_noEC | primary:high:P2b_nap_prev05 | 0.274841 | 0.274841 | 0.0e+00 |
| t_E_status_noEC | primary:medium:P1 | 0.194357 | 0.194357 | 6.7e-12 |
| t_E_status_noEC | primary:medium:P2b_nap_prev05 | 0.133404 | 0.133404 | 0.0e+00 |
| t_E_status_noEC | primary:low:P1 | 0.105254 | 0.105254 | 2.3e-11 |
| t_E_status_noEC | primary:low:P2b_nap_prev05 | 0.079257 | 0.079257 | 1.4e-17 |
| t_E_status_noEC | primary:spearman_evaluable | 0.646030 | 0.646030 | 1.3e-08 |
| b0 | P1_all | 0 | 0 | 0.0e+00 |
| b0 | P1_variable | 0 | 0 | 0.0e+00 |
| b0 | P1_selective | 0 | 0 | 0.0e+00 |
| b0 | P1_rare_selective | 0 | 0 | 0.0e+00 |
| b0 | dP1_vs_B1_variable | -0.166171 | -0.166171 | -4.4e-11 |
| b0 | dP1_vs_B2_variable | -0.249100 | -0.249100 | 1.5e-10 |
| b0 | P2a | 0.904784 | 0.904784 | 0.0e+00 |
| b0 | P2b_nap_prev05 | 0 | 0 | 0.0e+00 |
| b0 | P2b_auroc_prev05 | 0.500000 | 0.500000 | 0.0e+00 |
| b0 | P2b_n_prev05 | 1753 | 1753 | 0.0e+00 |
| b0 | primary:high:P1 | 0 | 0 | 0.0e+00 |
| b0 | primary:high:P2b_nap_prev05 | 0 | 0 | 0.0e+00 |
| b0 | primary:medium:P1 | 0 | 0 | 0.0e+00 |
| b0 | primary:medium:P2b_nap_prev05 | 0 | 0 | 0.0e+00 |
| b0 | primary:low:P1 | 0 | 0 | 0.0e+00 |
| b0 | primary:low:P2b_nap_prev05 | 0 | 0 | 0.0e+00 |
| b1 | P1_all | 0.087622 | 0.087622 | 1.2e-11 |
| b1 | P1_variable | 0.166171 | 0.166171 | 4.4e-11 |
| b1 | P1_selective | 0.155980 | 0.155980 | 1.5e-11 |
| b1 | P1_rare_selective | 0.095515 | 0.095515 | -4.8e-11 |
| b1 | dP1_vs_B2_variable | -0.082929 | -0.082929 | 1.9e-10 |
| b1 | P2a | 0.907638 | 0.907638 | 0.0e+00 |
| b1 | dP2a_vs_B0 | 0.002855 | 0.002855 | 8.2e-18 |
| b1 | P2b_nap_prev05 | 0.075128 | 0.075128 | 0.0e+00 |
| b1 | P2b_auroc_prev05 | 0.585677 | 0.585677 | 0.0e+00 |
| b1 | P2b_n_prev05 | 1753 | 1753 | 0.0e+00 |
| b1 | primary:high:P1 | 0.210921 | 0.210921 | -7.7e-11 |
| b1 | primary:high:P2b_nap_prev05 | 0.120128 | 0.120128 | 0.0e+00 |
| b1 | primary:medium:P1 | 0.101181 | 0.101181 | -2.1e-12 |
| b1 | primary:medium:P2b_nap_prev05 | 0.051040 | 0.051040 | 0.0e+00 |
| b1 | primary:low:P1 | 0.059336 | 0.059336 | 3.5e-11 |
| b1 | primary:low:P2b_nap_prev05 | 0.033518 | 0.033518 | -6.9e-18 |
| b1 | primary:spearman_evaluable | 0.438318 | 0.438318 | -9.1e-09 |
| b2 | P1_all | 0.121805 | 0.121805 | -2.1e-11 |
| b2 | P1_variable | 0.249100 | 0.249100 | -1.5e-10 |
| b2 | P1_selective | 0.228211 | 0.228211 | 2.1e-11 |
| b2 | P1_rare_selective | 0.135272 | 0.135272 | -1.2e-12 |
| b2 | dP1_vs_B1_variable | 0.082929 | 0.082929 | -1.9e-10 |
| b2 | P2a | 0.909839 | 0.909839 | 0.0e+00 |
| b2 | dP2a_vs_B0 | 0.005055 | 0.005055 | -1.7e-18 |
| b2 | P2b_nap_prev05 | 0.139645 | 0.139645 | 0.0e+00 |
| b2 | P2b_auroc_prev05 | 0.625136 | 0.625136 | 0.0e+00 |
| b2 | P2b_n_prev05 | 1753 | 1753 | 0.0e+00 |
| b2 | primary:high:P1 | 0.312713 | 0.312713 | -1.8e-10 |
| b2 | primary:high:P2b_nap_prev05 | 0.216421 | 0.216421 | 0.0e+00 |
| b2 | primary:medium:P1 | 0.148720 | 0.148720 | -3.6e-11 |
| b2 | primary:medium:P2b_nap_prev05 | 0.100703 | 0.100703 | -1.4e-17 |
| b2 | primary:low:P1 | 0.074369 | 0.074369 | 1.5e-11 |
| b2 | primary:low:P2b_nap_prev05 | 0.060372 | 0.060372 | 0.0e+00 |
| b2 | primary:spearman_evaluable | 0.590013 | 0.590013 | 4.7e-09 |

all 298

## Task 2: red team

Verified clean:
- **Unblinding order** (git reflog shows no rewrite):
  - 84efb5c, finalists.json, 00:00:03;
  - preflight, 00:00:35 (all checks pass);
  - predictions written 00:29–00:43;
  - bf3faa3, manifest commit, 00:43:08;
  - first access 00:54:53.
- **One access per finalist:**
  - Access-log rows 919–925 are exactly 7 granted test target accesses, one per config_hash, each preceded by a log
    commit (beeb32d … 6448cda); 1cb4981 closes the run.
  - No access was refused.
  - The only earlier test access is the mask_only count read on 2026-10-26.
- **Configs:**
  - Every config hashes to its committed value.
  - finalists.json equals finalists.example.json apart from _comment.
  - No change in configs/, src/ or scripts/ between 84efb5c and HEAD.
- **Frozen files and split:**
  - The frozen files match frozen_hashes.json and are unchanged since freeze-eval.
  - splits.csv has a single commit (3834bc9).
  - No PatientID crosses roles, and none is missing.
- **Gene sets:** eval_genes and gene_classes were committed at G0 (5181e51); predictable_genes in 6c5e004, before
  finalists.json.
- **Training paths:**
  - refit_predict loads dev Y only (heldout_common.py:406).
  - guarded_fit passes training rows only (run.py:145–156).
  - FoldTransformer and check_fit_ids (features.py:391–439) cover the SD filter and scaling.
  - The kernel is built and centred on training rows only (kernel_ridge.py:680, 726).
  - B2's top-variance genes and scaling come from training rows (baselines.py:133–146).
  - M3 screening and early stopping use training rows (lgbm.py:197–199).
  - The M5 inputs are refit from the training Y with fresh caches (evaluate_test.py:293–299).
  - The train_ids sha256 values equal the dev set (and dev minus the 4 EC lines for no-EC).
  - All phase-A load_targets calls request dev IDs only. A test request would have been logged before the read
    (guard.py:213–223), and none was logged before 00:54:53.
- **Features:** no CRISPR or RNAi data. Paralogs come from Ensembl Compara; ORCS is rejected (priors.py:187–207);
  partners are unused.
- **Population:**
  - Lineages: 27 of 28 appear in test, none is test-only; χ² p 0.18; maximum share difference 3.0 points.
  - Library/NaN groups: test 167/32/36/9/11 vs dev 621/154/143/47/57; χ² p 0.60.
- **Near-duplicates:** the maximum expression correlation from a test line to a dev line (median 0.919, maximum 0.951)
  is no higher than between dev lines of different patients (median 0.918, maximum 0.969). The top pairs are
  distinct lines.

## Findings

Critical: none. Major: none.

Minor:
1. **Spearman(dev mean_r, test r) CI.**
   - The normal CI [0.622, 0.673] is centred on 0.6475.
   - The bootstrap mean is 0.084 lower, and the percentile interval is [0.537, 0.587]
     (results/test/m5_avg/summary.json per_tier).
   - Report the CI as approximate.
2. **M5 assertion coverage.** input_fit_log_ covers only the outer input refits (m5_stack.py:260–267); the nested
   fits are not checked. This is a gap in the check, not a leak.
3. **Eval genes.** The eval-gene set used test NaN counts (mask-only, pre-registered).
4. **DEEP_DEL_THR.** It was set from the features of all lines (scripts/build_features.py:32–51). No targets were
   involved.
5. **Chronos.** The joint Chronos fit is disclosed; keep it in REPORT limitations.
6. **dP1 vs B1, test vs CV.**
   - Test +0.153 vs CV +0.134; the difference comes mostly from B1 (0.166 on test vs 0.184 in CV, inside the CV
     fold range of 0.163–0.197).
   - m5_avg barely moves: 0.319 vs 0.3175.
   - The test refit trained on 1022 lines vs about 818 in CV.
   - Quote absolute P1 next to dP1.
7. **Variable genes only.** m5_avg equals m1 exactly on the 16,275 non-variable genes, so P1_all differs by only
   0.0011.
8. **Static test scope.** test_no_direct_crispr_reads.py does not scan reviews/; this read was logged by hand.

## Bottom line

- The protocol was followed, and no test line entered any fit.
- All numbers reproduce exactly.
- Report m5_avg P1 0.3194 (dP1 vs B1 +0.153, vs B2 +0.070) and the paired gain over m1 of +0.0093
  [0.0085, 0.0100], with the caveats in findings 1, 6 and 7.
