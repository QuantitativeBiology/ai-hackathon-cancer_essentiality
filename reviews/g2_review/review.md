# Independent review of G2 (kernel ridge M1/M1b/R1/T, ablation, M2)

Reviewer: independent agent, 2026-10-26. Code reviewed at commit e44a4c8 (plus HEAD `run.py`, whose
later change only affects model saving). I wrote my own metric code, as LESSONS 10 asks
(`essml.evaluate` / `essml.metrics` are not imported). CRISPR values were read only through
`essml.guard.load_targets` with dev IDs. No test or noncancerous line was read. No frozen file,
config, registry row or result was changed.

Scripts (all in `reviews/g2_review/`; outputs are numbers only, no data):

| script | what it does | output |
|---|---|---|
| `recompute_p1.py` | Recomputes P1 from the stored OOF predictions, plus Δ/SE, per-line centring, class-wise centring, offset-only models, within-library and within-lineage centring, and training-fold variable sets | `recompute_p1.json` |
| `neg_genes.py` | Hart NEGv1 puzzle: r distribution; own-gene expression/CN; partial correlations | `neg_genes.json` |
| `neg_lowrank.py` | NEGv1 puzzle: projects out training-fold PCs of the NEG block | `neg_lowrank_nonessential.json` |
| `check_solver.py` | Brute-force check of the pattern ridge and closed-form leave-group-out with p > n, merged and fallback patterns, patient groups straddling NaN boundaries, and a B1 residual; also the B1 closed form vs a full-pipeline LOGO | stdout (PASS) |
| `check_b2.py` | B2 vectorised predict vs the old loop | stdout (PASS) |
| `reproduce_fold.py` | Refits fold 2 with the committed code, with non-training feature rows removed from the raw cache during fit, and compares with the stored predictions | `reproduce_fold2.json` |

The stored OOF predictions were copied from `cache/runs/` to a scratch directory before a
concurrent cache cleanup. Runs used: b0 150901, b1 163644, b2 164759, m1_E 170015, r1_E 171800,
t_E_status 193048, m1_E_resB0 202213, m1b_E_C 173206, m1b_full 191141.

---

## Verdict

| # | finding | severity | evidence | concrete fix |
|---|---|---|---|---|
| 1 | No validation/test CRISPR reaches any fit, λ/κ choice, offset or M2 prior. Features are fitted on training rows; patient groups are respected in the outer folds, κ CV, ridge LOGO and M2 inner CV. | none | §1; fold 2 refitted with non-training feature rows deleted reproduces stored predictions (max \|Δ\| 2.4e-7) | — |
| 2 | B1-residual LOO criterion is not an exact leave-group-out of the pipeline. The offset is leave-line-out, and other lines' residuals contain y_i. | minor | Toy n=300: criterion biased by −0.02 to −0.03 r and picks λ 1e5 vs 1e4 (`check_solver.py`). B0 path is exact. | Accept for the B0 path (D14). If the B1 arm is ever revived, use a leave-patient-out offset and recompute the residuals per held-out group, or tune λ by inner grouped CV. |
| 3 | Pattern-wise fit and closed-form LOO/LOGO are correct. | none | Brute force with p>n, merged/fallback patterns and groups across NaN boundaries: max \|Δr\| 5e-13, predictions 9e-12 | — |
| 4 | Unit tests are too weak for the production regime: p<n only, no R≠Y LOO, launcher tests run in frame mode only, held-out *features* never perturbed. | minor | §2d | Add `check_solver.py` part 1 and a real-path held-out-feature perturbation test to `tests/`. |
| 5 | B2 vectorisation is equivalent to the old loop. | none | max \|Δ\| 4e-16 incl. edge cases | — |
| 6 | Reported P1 and Δ/SE are reproduced exactly by independent code. | none | §3a–b (all to 4 decimals) | — |
| 7 | Provenance: all G2 runs and the G2 gate came from dirty trees on commits that predate the model code. The re-runs excluded m1_E, R1 and resB0. | major → mitigated | Registry `dirty_tree=1`; commits 347cde3/083578a/987bdec. Committed code reproduces fold 2 of m1_E, m1_E_resB0, R1 and T to 2.4e-7. | Re-run the finalists (m1_E_resB0/T, R1, b1, b2) on a clean tagged tree after `freeze-eval`, and cite those run IDs. |
| 8 | D16 hypothesis (NEG skill = predicted per-line offsets) is refuted. NEG skill is broad, partly library and partly lineage, not low-rank, and proportional to every model's skill, B1 included. | minor (for the headline) | §4: per-line centring keeps 99–100%; offset-only models give ≤0.024 on NEG; within-library 0.097, partial-lineage 0.081, top-20 PCs removed 0.050 | Replace the D16 per-line-centring diagnostic with within-library-group and within-lineage P1, reported beside P1. Do not treat a null line-offset result as clearance. |
| 9 | Library/batch-predictable structure inflates absolute P1. | minor | Within-library: M1 0.272 (91%); ΔM1−B1 +0.101 vs +0.116 | Report within-library P1/ΔP1 as the pre-registered "by NaN/library group" column in the headline table. |
| 10 | The headline mixes runs: "0.308" is m1_E_resB0/T (post-hoc D14), while the doc and gate select m1_E (0.2997). | minor | docs/g2_results.md vs RESEARCH_TREE status | Name one model and one run ID everywhere. Add resB0/T to the selection table as a post-D14 candidate with the D14 note. |
| 11 | The permuted control is easy. A global shuffle removes lineage and library information too; it was run only on the B1 path, with λ at the grid edge. | minor | §5c | Add expression permuted within lineage (and lineage × library); expect ≤ B1. Add a permuted B0-path run; expect ≈ B0. |
| 12 | Fold SE is optimistic (overlapping training sets, df=4). | minor | §5f | Show Nadeau–Bengio-corrected SE (×≈1.5) beside the SE. No current conclusion changes. |
| 13 | Runtimes are inflated about 1.4–1.5× by contention, and tuning-inclusive times are not comparable across models. | minor | Fold 2 refit 58–65 s at 6 threads vs registry 80–100 s | Re-time finalists alone (LESSONS 14). Report fit and tuning separately. |
| 14 | Dev-wide vs training-fold variable-gene set makes no difference. | none | Δ ≤ 0.0003 for all models | — |
| 15 | M2: local-λ CV reuses the OOF-generating splits, and the base is still the B1-residual M1. | minor | m2_local.py:186–195, 122 | Nest the OOF for λ, or accept as training-only. Re-run M2 on the D14 base. |

**Can the G2 numbers be relied on?** Yes, for dev-CV model comparison:
- No leakage was found.
- The solver is correct.
- Every headline P1 and Δ/SE reproduces exactly with independent code.
- The committed code reproduces the stored predictions of the four E-only models.

The caveats are about framing, not validity:
- Quote one model (m1_E_resB0 or T, 0.308) with its run ID.
- Report the within-library value next to it (0.279; ΔP1 vs B1 +0.108).
- Re-run the finalists on a clean tree before freezing.
- Treat runtimes as indicative only.

---

## 1. Leakage

**Verdict: no leakage of validation- or test-line CRISPR into any fit or selection step was found.**
One approximation in the B1-residual LOO λ criterion uses training data only (see 1e).

**a. Feature transforms.**
- `FoldTransformer.fit` (features.py:428–439) selects training rows before computing any
  statistic. This covers SD filter, min-count filter, means/SDs, imputation means, categorical
  levels, and Hallmark/p53 score means and SDs (`_fit_score`, :470).
- The raw `build_*` blocks are row-wise. The one exception is the fusion pre-filter (≥5 lines
  overall, features.py:137). It counts feature presence including non-dev lines, but it is a
  superset of the fit-time ≥5-training-lines filter, so it changes nothing.
- `check_fit_ids` (:391) refuses test and noncancerous IDs. It does not know the outer fold. For
  M1, though, `train_ids = Y_tr.index`, and Y_tr has passed the run.py sentinel.
- Empirical check (`reproduce_fold.py`): fold 2 was refitted with the raw E (and status) rows of
  every non-training line *deleted* from the feature cache during `fit`, then restored for
  `predict`.
  - The predictions equal the stored G2 predictions to float32 precision (max |Δ| 2.4e-7 for
    m1_E and m1_E_resB0; R1/T see §3).
  - So no held-out or test feature row reaches any fitted quantity.

**b. B0/B1 offsets.**
- `_offset_fit` (kernel_ridge.py:397) fits `GeneMean` / `LineageMean` on `Y_tr` only.
- κ comes from `tune_kappa` (:294): StratifiedGroupKFold on training lines, grouped by
  PatientID, with variable genes from the inner-training part.
- The leave-line-out residual (baselines.py:94) uses training rows only.

**c. λ / κ selection.**
- λ comes from the closed-form LOO on the training-fold kernel and training-fold targets
  (kernel_ridge.py:178–214, 217–281). Gene groups (variable/other) come from the training fold
  (:409).
- The outer held-out fold never enters, so the closed form cannot leak the validation fold. It
  sees only `Y_tr`, and run.py:228 slices `Y.loc[tr]` before any model code runs.

**d. M2 priors.**
- Paralog and partner tables are external and line-independent. `LocalFeatureBuilder`
  (priors.py:234) is row-wise and has no CRISPR access.
- Per-gene standardisation uses observed training rows (`normal_equations`).
- The inner OOF M1 predictions and the local-λ inner folds are both inside the training fold.

**e. Patient groups.**
- Outer folds: no PatientID spans two folds or two roles (checked on `configs/splits.csv`; 72
  dev lines share a patient).
- `tune_kappa`, run.py inner CV and M2 inner splits all use StratifiedGroupKFold by PatientID.
- The ridge LOO is leave-patient-out: block formula over patient groups inside each cluster,
  verified by brute force.
- *Approximation (minor):* the B1 offset is leave-**line**-out, not leave-patient-out. Also,
  other training lines' residuals R_j still contain y_i through the lineage means. So for the
  B1-residual model (pre-registered M1) the closed-form criterion is not an exact
  leave-group-out of the whole pipeline.
  - On a realistic toy (n=300, 6 lineages of 20–90 lines, κ=10) the criterion is biased
    downwards by 0.02–0.03 in r. It selects λ=1e5 where a true full-pipeline LOGO selects 1e4
    (`check_solver.py`).
  - This is training-only, so it is not leakage, but it handicaps the B1-residual arm (see 5d).
  - An extreme toy (lineage of 6 lines, row shifts confounded with lineage) gives a gap of up to
    0.10. How large it is on the real data is unknown.
  - The B0-residual models (R1, T, m1_E_resB0, the D14 choice) are exact: the unpenalised
    intercept absorbs the gene mean.

**f. CRISPR from validation/test lines as model input.** None.
- Dev-wide `gene_classes.csv` (which includes held-out lines) is used only for reporting, as
  pre-registered.
- `eval_genes.csv` uses mask-only test counts.

---

## 2. Correctness

**a. Pattern-wise exact fit.**
- Each cluster's kernel is the submatrix of one fold-level K, double-centred on the cluster
  rows. Prediction is `rbar + (k − kbar) α`. This equals ridge with an unpenalised intercept on
  the cluster rows.
- `check_solver.py` confirms it against primal ridge at new points (max |Δ| 9e-12). The test
  covers p > n (150 features, 90 rows), rows inside NaN regions shifted by +25 to catch
  centring bugs, a merged pattern (intersection rows), an unmerged near-duplicate (the merge
  guard at :129 works as documented) and a fallback pattern (prediction exactly 0).

**b. LOO / LOGO formula.**
- The hat matrix is 11ᵀ/n + V diag(d/(d+λ)) Vᵀ. LOO residuals are e/(1−h_ii); the group blocks
  use (I − H_gg)⁻¹ e_g.
- The criterion is r(y − e, y) from per-gene sums.
- Against brute-force leave-group-out with fixed targets, the max |Δr| is 5e-13 (B0 residual)
  and 3e-13 (B1 residual). Groups straddle NaN-pattern boundaries.
- With the residual offset, p = (y − R) + ridge_LOO is the intended quantity, and it is computed
  correctly. Its limitation is 1e above.

**c. B2 vectorisation.**
- It is equivalent to the old loop: max |Δ| 4e-16, including all-NaN neighbourhoods, k=1,
  k=n_train and k>n_train (`check_b2.py`).
- Neighbour indices from argsort are distinct, so the 0/1 weight matrix is exact.
- The G1 and G2 B2 rows agree to 5 decimals on P1, P2a and P2b.

**d. Are the unit tests strong enough? Not fully. Gaps:**
- All solver tests have p < n (p=8, n=60). The real regime is p ≈ 18k ≫ n ≈ 800.
- LOO is never tested with R ≠ Y (B1 residual).
- Group LOO is tested only on fully observed rows with `merge_frac=0`.
- Merging is tested for prediction but not for the LOO criterion.
- The launcher leakage test runs in *frame mode* (X passed as a DataFrame). The production path
  (X=None, raw blocks from cache/features, `FoldTransformer` with `check_roles=True`) is never
  exercised.
- Held-out invariance perturbs held-out *targets* only, never held-out *features*.
- `tune_kappa`, `permute_seed` and the T modality guard are untested.

My scripts cover the solver gaps (PASS) and the held-out-feature gap (fold-2 reproduction).
Recommendation: move `check_solver.py`'s first part and a small real-path feature-perturbation
test into `tests/`.

---

## 3. Results integrity

**a. Independent P1 recomputation.**
- Method: per-gene masked Pearson within each outer fold, genes ∈ eval_genes ∩ class with ≥40
  observed held-out lines, constant → 0, mean over genes, then over folds.
- It reproduces every number I checked, to 4 decimals:

| run | P1 variable (mine) | doc | P1 all | P1 NEGv1 |
|---|---|---|---|---|
| b1 | 0.1837 ± 0.0131 | 0.1837 ± 0.0131 | 0.0928 | 0.0676 |
| b2 | 0.2527 ± 0.0087 | 0.2527 ± 0.0087 | 0.1205 | 0.0827 |
| m1_E | 0.2997 ± 0.0067 | 0.2997 ± 0.0067 | 0.1498 | 0.1104 |
| r1_E | 0.3073 ± 0.0050 | 0.3073 ± 0.0050 | 0.1546 | 0.1140 |
| t_E_status | 0.3082 ± 0.0040 | 0.3082 ± 0.0040 | 0.1548 | 0.1143 |
| m1_E_resB0 | 0.3081 ± 0.0047 | 0.3081 ± 0.0047 | 0.1547 | 0.1141 |

- Fold-mean per-gene r also matches the committed `per_gene.csv.gz` within its 4-decimal
  rounding (max |Δ| 5e-5).

**b. Δ / SE** (paired fold differences, SE = SD/√5). All match docs/g2_results.md:
- m1_E − b1: +0.1159 (0.0036), 5/5 folds.
- m1_E − b2: +0.0470 (0.0021), 5/5.
- resB0 − m1_E: +0.0084 (0.0017), 5/5.
- r1 − m1_E: +0.0076 (0.0019).
- T − m1_E: +0.0085 (0.0017).
- T − resB0: +0.0001 (0.0004).
- resB0 − r1: +0.0008 (0.0008).
- m1b_E_C − m1_E: −0.0135 (0.0005).
- Permuted control (from `per_fold.csv`; predictions not cached): ΔP1 vs B1 = +0.00002.

**c. Provenance (major, now largely mitigated).**
- Every G2 registry row, and the G2 gate itself, has `dirty_tree=1`. Their recorded commits
  (347cde3, 083578a, 987bdec) predate `kernel_ridge.py` and `m2_local.py`, which first appear in
  e44a4c8. So no stored number was tied to committed code.
- The "19 bit-identical re-runs" were also dirty, sometimes on the same commit. They did not
  include m1_E, r1_E or m1_E_resB0, i.e. the selected model and the D14 comparison.
- I refitted fold 2 of m1_E, m1_E_resB0, r1_E and t_E_status with the committed code. They
  match the stored predictions to float32 precision (table in §6). Committed configs hash to the
  run IDs for all 29 latest runs.
- Residual risk: M1b/M2 rows are unverified against committed code; their same-commit
  duplicates agree with each other.

---

## 4. The non-essential-gene puzzle (D16)

**Per-line offsets do not explain it. The D16 hypothesis is not supported.**

**Per-line centring** (line mean over observed eval genes, subtracted from both observed and
predicted values, within each fold):

| run | P1 var raw → centred (survives) | P1 NEG raw → centred (survives) |
|---|---|---|
| b1 | 0.1837 → 0.1838 (100%) | 0.0676 → 0.0674 (100%) |
| b2 | 0.2527 → 0.2530 (100%) | 0.0827 → 0.0824 (100%) |
| m1_E | 0.2997 → 0.3000 (100%) | 0.1104 → 0.1098 (99.5%) |
| r1_E | 0.3073 → 0.3076 (100%) | 0.1140 → 0.1133 (99%) |
| t_E_status | 0.3082 → 0.3085 (100%) | 0.1143 → 0.1136 (99%) |
| m1_E_resB0 | 0.3081 → 0.3084 (100%) | 0.1141 → 0.1134 (99%) |

- **Class-wise centring** (line mean over the evaluated class's own genes): variable m1_E
  0.3006, NEG 0.1094 (99%). No change either.
- **Offset-only models:**
  - B0 + M1's predicted per-line shift over all genes: P1 var 0.006, NEG 0.012.
  - The same shift per class: NEG 0.024.
  - Oracle observed line mean (upper bound): var 0.054, NEG 0.034.
- **Line-mean scale:** the predicted line means vary very little (SD 0.002–0.003 vs observed
  0.008). Their correlation with observed line means is 0.20 over all genes and 0.40 over NEG
  genes. B0 alone already gets 0.34 on NEG, via NaN-pattern composition.

**What the NEG skill is:**
- **Broad:** 89% of the 804 NEG genes have r > 0; median r 0.106.
- **Not own-gene driven:** r(own expression, effect) = −0.012 and r(own CN, effect) = −0.041.
  Partialling both out of observed and predicted values leaves 0.105.
- **Partly lineage:** additionally partialling lineage leaves 0.081. Within-lineage centring
  gives 0.086.
- **Partly library:** within NaN/library-group centring (4 largest row-NaN patterns + other)
  gives 0.097.
- **Diffuse, not low-rank:**
  - The training-fold NEG block has no dominant factor (PC1 3.3% of variance).
  - Projecting out the top 1/3/10/20 training PCs leaves M1 0.106/0.086/0.058/0.050 (B1
    0.065/0.052/0.032/0.025).
  - PC2 is 35% library and 26% lineage (R² of held-out scores). M1 predicts the PC2 scores with
    r 0.55.
- **Not M1-specific:** the NEG/variable P1 ratio is the same for B1 (0.37), B2 (0.33), M1-E
  (0.37) and m1_E_resB0 (0.37). Every model, including the pure lineage mean, gets NEG skill in
  proportion to its general skill.

**Interpretation:**
- Chronos scores for NEGv1 genes still carry line-to-line variation that is shared across many
  genes. It is structured by library/screen group, by lineage, and by a diffuse multi-factor
  component that expression predicts.
- Some of it is technical: the library group, and PC2's library loading. The rest cannot be
  separated from biology with dev data alone.
- It is not a per-line scalar, so a "B0 + predicted line scalar" model and per-line-centred P1
  (the D16 plan) will come out null and should not be read as clearance.
- Rough scale of the effect on the headline:
  - If the NEG-type structure had the same *absolute* size in variable genes (NEG SD ≈ 0.13 vs
    variable ≈ 0.33), it could account for at most about 0.11 × 0.13/0.33 ≈ 0.04 of variable P1.
  - Within-library centring removes 0.028 of M1's P1 and 0.015 of ΔP1 vs B1.

---

## 5. Other threats to the headline

**a. The headline mixes runs.**
- "M1-E/R1/T ≈ 0.308":
  - m1_E (the pre-registered M1, which docs/g2_results.md names as the selected model) is
    **0.2997**;
  - R1 is 0.3073;
  - 0.308 is m1_E_resB0 and T, chosen post hoc in D14.
- RESEARCH_TREE quotes 0.308 while the G2 doc and the gate say "selected: m1_E". The selection
  table lists neither resB0 nor T as candidates.
- The post-hoc optimism is negligible: the top three are within 0.0009. Still, the docs should
  name one number and one run ID, per LESSONS 5.

**b. Library/batch-predictable structure.**
- Within-library centring keeps 91% of M1's variable P1 (0.2717). ΔP1 vs B1 drops to +0.101
  (SE 0.003, 5/5 folds) and vs B2 to +0.046.
- Within-lineage P1 is 0.237 for M1, 0.177 for B2 and 0 for B1 by construction. The skill
  beyond lineage is real.
- The headline should carry the within-library number as a secondary column, as DESIGN §4
  already pre-registers "results by NaN/library group".

**c. The permuted control is easy.**
- It shuffles ModelIDs globally, so features lose lineage and library information as well. It
  was run only for the B1-residual M1, and λ hits the grid's upper edge (1e7).
- It shows the LOO λ does not fit noise. It does not test the B0-residual headline path, or
  whether skill exceeds lineage/library proxies.
- Better controls:
  - expression permuted *within lineage*: should not beat B1;
  - expression permuted within lineage × library group;
  - a permuted run of the B0 path, which should give P1 ≈ 0.
- The gate threshold (ΔP1 ≤ 0.01) is lenient but was easily met.

**d. D14 interpretation.**
- The decision (switch to the B0 residual) is safe: it is evaluated on outer folds, and the B0
  LOO criterion is exact.
- The stated reason ("a fixed shrunk lineage mean is a worse offset") is not established. The
  B1-residual arm also has an inexact, downward-biased λ criterion (1e) that can pick a
  different λ.

**e. Variable-gene set definition.**
- Using the dev-wide classes (which include held-out lines) instead of training-fold classes
  changes P1 by ≤0.0003 for every model (m1_E 0.2994 vs 0.2997). Not a concern.

**f. Fold SD vs SE.**
- The table correctly shows P1 ± fold SD and Δ ± SE with t₄ intervals.
- Fold differences share about 75% of the training data, so SD/√5 understates uncertainty. A
  Nadeau–Bengio-style correction inflates the SE about 1.5×.
- No headline or D14 conclusion changes (D14: t ≈ 4.9 → ≈3.3). The "Δ > 1 SE" rule is weak with
  df = 4.

**g. Runtime.**
- On a lightly loaded machine, the committed code fits fold 2 in 58–65 s at 6 threads (table in
  §6). The registry reports 80–100 s for the same fold: m1_E 100 s, r1 85 s, T 80 s at 12
  threads; resB0 80 s at 6 threads.
- So the G2 runtimes are inflated about 1.4–1.5× by contention, as LESSONS 14 says. More
  threads did not help.
- fit_s includes inner tuning (B1/B2 grid CV; M1-B1 κ tuning), so per-model costs are not
  comparable. Re-time finalists alone before quoting any runtime.

**h. M2 details (minor).**
- The local-λ inner CV reuses the same inner splits that produced the OOF M1 features. So
  training-row features in split a were built from models that saw split b. This is training-only
  optimism in λ, not leakage.
- M2 still uses the B1-residual M1 base. Re-run on the D14 base.

---

## 6. Reproduction of fold 2 with committed code

Non-training feature rows were removed during fit. Fit times are at 6 threads.

| run | max abs diff vs stored | λ (other/var) | κ | fit s |
|---|---|---|---|---|
| m1_E | 2.4e-7 | 21544 / 21544 | 20 | 64.8 |
| m1_E_resB0 | 2.4e-7 | 46416 / 21544 | — | 57.9 |
| r1_E | 2.4e-7 | 46416 (single λ) | — | 58.0 |
| t_E_status | 2.4e-7 | 46416 / 21544 | — | 58.1 |
