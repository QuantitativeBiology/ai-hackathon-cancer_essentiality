# G3 nonlinear track — M4 pre-check / M4, M3 LightGBM, R2/R4, published-comparable block

Status: **pre-registration** (written and committed before any run of this track). Results are
appended below the pre-registration; the pre-registered text is not edited afterwards (changes are
logged in "Deviations").

Current best (D17): `m1_E_resB0`, run 20261026-222117, dev-CV P1 (variable genes) 0.3081.
Selection rule (DESIGN §4): a more complex model replaces it only if the paired fold-wise ΔP1
(variable genes, dev-wide `configs/gene_classes.csv`) > δ = 0.005 **and** > 1 SE (SE = SD of the 5
fold differences / √5). All runs go through `scripts/run.py` (registry), dev lines only, 6 cores.
Runtimes are indicative (machine shared with the kernel track, D18 / lesson 14).

## Pre-registration

### 1. M4 pre-check (DESIGN §5; time box for all of M4: 1 day)
- Inputs (both arms identical): expression (E) → `FoldTransformer` fitted on the training lines →
  PCA fitted on the training lines (256 PCs; held-out lines only projected). Target: y standardised
  per gene with the training mean/SD on observed rows (= B0 residual / SD). NaN targets are masked,
  never filled.
- Inner split: ONE patient-grouped, lineage-stratified 80/20 split of the outer training lines
  (`split_seed` 20261026), used for early stopping and for the small grids of both arms. During
  the inner phase the feature transform, PCA and target standardisation are refitted on the
  inner-training lines only. Criterion: mean per-gene r on the inner-validation lines over the
  variable genes of the inner-training part.
- Arm A, `M4pre_RidgePC`: multi-output ridge on the PCs (exact per NaN pattern), grid λ/n ∈
  10^[-4..2] (13 points) × PC whitening {no, yes}. Deterministic.
- Arm B, `M4pre_MLP`: 256 → 512 (ReLU, dropout) → G, masked MSE, AdamW (lr 1e-3, weight decay 1e-4),
  batch 64, ≤ 200 epochs, patience 20 on the inner criterion; grid whitening {no, yes} × dropout
  {0.2, 0.5}. Refit on all training lines for the chosen number of epochs. One run per model seed
  0–4 (`configs/seeds.json`); the seed sets initialisation and minibatch order only.
- Both arms are run on outer fold 0 only (`--folds 0`).
- **Decision: GO for M4 iff the mean over the 5 seeds of the MLP's fold-0 P1 (variable genes) is
  greater than RidgePC's fold-0 P1. Otherwise M4 is dropped (NO-GO) and not built.**
  The fold-0 P1 of `m1_E_resB0` is reported as context only.

### 2. M4 (only if the pre-check says GO)
- y[c,g] = b_g + ⟨u(x_c), v_g⟩: encoder 256 PCs → 512 → 64; v_g initialised from the top
  right singular vectors of the training-fold standardised target matrix; masked MSE; early
  stopping on the same inner grouped split; CPU. No local term and no prior-knowledge
  initialisation (time box).
- 5 outer folds × 5 model seeds. Selection uses the seed-averaged fold-wise P1; the seed spread
  (SD over seeds of the 5-fold mean P1) is reported. Permuted-feature control (rows of E permuted
  across lines) if put forward.

### 3. M3 LightGBM on m1_E_resB0 residuals (D18; time box 1 day; 6 processes × 1 thread)
- Per outer fold: out-of-fold `m1_E_resB0` predictions for the training lines from 5 inner
  patient-grouped folds of the training lines; residual r = y − m1_oof. One LightGBM per
  training-fold variable gene, fitted on r; prediction = m1_E_resB0 (refit on all training lines)
  + LightGBM correction. Non-variable genes pass m1_E_resB0 through unchanged.
- Features per gene: inner-fold screening (top-k |Pearson r| with the residual, computed only on
  the lines the model is fitted on) from a pool of E, C, M, H and lineage, plus the gene's own
  M1 prediction. Early stopping on an inner patient-grouped split; screening is redone inside
  that split. Exact settings are in `configs/runs/g3_m3_*.json`, committed before the run.
- Coverage rule: all training-fold variable genes if a timing probe on fold 0 (100 genes of the
  subset below) projects ≤ 4 h wall for 5 folds; otherwise the pre-declared random subset
  `configs/g3_variable_subset500.csv` (500 of the 2,118 variable eval genes, seed 20261026),
  compared with `m1_E_resB0` on the same 500 genes.
- Negative control: the same M3 with the tree features permuted across lines (`permute_seed`);
  expected ΔP1 vs m1_E_resB0 ≤ 0.

### 4. Published-method re-runs
- **R2** — Shi-style elastic net, α (l1_ratio) = 0.5, per variable gene, on standardised E; the
  penalty is chosen per gene by inner patient-grouped CV (MSE). Raw y with intercept. Non-modelled
  genes get B0. Time box 2 h: a timing probe on fold 0 decides; if the full 5-fold run projects
  > 2 h wall, R2 runs on `configs/g3_variable_subset500.csv`.
- **R4** — TrPLet-like: per gene, top-k E features by |r| with y on the fitting lines, k ∈ {5, 10,
  20, 50, 100, 200} × learner ∈ {ridge (LOO-tuned), RBF-SVR} = 12 candidates. **Nested**: the
  candidate is chosen per gene by inner 3-fold patient-grouped CV inside the training fold (the
  top-k screen is recomputed inside each inner split). **Non-nested** (leaky=1, diagnostic only):
  the same 12 candidates fitted on the outer training folds, and the candidate per gene chosen by
  its mean outer-CV r — TrPLet's best-of-N selection inside the reported CV. Optimism =
  P1(non-nested) − P1(nested). Same coverage rule as R2.

### 5. Published-comparable block ("loosely comparable"; never used for selection)
For B1, B2, R1, m1_E_resB0 and the G3 runs, from the stored per-fold per-gene results:
- mean per-gene Spearman on the DeepDEP 1,298-gene and LEAP 1,539-gene sets, **only if** the gene
  lists can be obtained from published supplements or code repositories without bot-check
  circumvention (otherwise stated as not obtainable);
- TrPLet-style count of genes with r ≥ 0.2 (all eval genes);
- Shi-style count: genes with r > 0.2 among genes with ≥ 5 held-out lines in each label class
  (effect < −0.5 vs ≥ −0.5).
Counts are per fold and then averaged (D3); pooled-OOF versions are reported beside them.

### Reporting
For every run: P1 variable ± fold SD, ΔP1 vs m1_E_resB0 ± SE (paired, 5 folds), P1 all, P2b nAP and
AUROC (selective genes, prevalence ≤ 0.5), fit time. Subset runs are compared with m1_E_resB0 on
the same genes from its stored per-gene results.

## Deviations
- **D-G3-1 (M3 early stopping, after the independent leakage review, before any full-M3 result).**
  The reviewer found that M3's early-stopping training rows got M1 features and residual targets
  from out-of-fold M1 models that had seen the early-stopping validation rows' targets (the split
  coincided with OOF fold 0). This is inside the training fold, so it is not outer leakage, and toy
  simulations showed it is pessimistic (it suppresses corrections). Fix (`nested_es=True`, configs
  `g3_m3n*`): the early-stopping split is OOF split 0; its training rows get M1 features from a
  nested OOF pass inside them; its validation rows keep the M1 fitted on the training rows only.
  The first full run `g3_m3` (20261026-235648) was killed during fold 0 and marked failed; it
  produced no numbers. Seen before this decision: only the 100-gene fold-0 probe.
- **D-G3-2 (negative control).** `permute_seed` originally permuted rows over every ID in the
  feature cache (including test/set-A/unscreened lines' *features*; no targets). Fixed before any
  control was run: rows are permuted only among the outer-training lines present in each block;
  held-out lines keep their own features (DESIGN §8 "features permuted across training lines").
- **D-G3-3 (R2 coverage).** The R2 timing probe (`g3_r2_probe`, 30 genes, fold 0, 40 penalties,
  5 inner folds) had not finished after 21 min at 6 processes (≈ 220 s per gene per process) and
  was killed (registry row marked failed). Projected: full ≈ 108 h, 500-gene subset ≈ 25 h, both far
  beyond the 2 h box. A fold-0 micro-benchmark on 6 subset genes (scratch script, training lines
  only; not a result) gave ≈ 45 s per gene with 20 penalties down to 0.05 × the largest and 3
  inner folds (chosen penalties were 0.03–0.07, never at the path edge); restricting to the 5,000
  most variable genes saved only ~1.4×. To stay near the box with the exact method, R2 runs on the
  **first 150 genes of the pre-declared random subset** (≈ 1.6 h): config `g3_r2_sub150`. This
  replaces the pre-declared 500-gene fallback; the genes are still a random sample fixed before
  any run. The R4 timing probe was also interrupted (queue re-planning; no result) and is re-run.

## Results

### R1. M4 pre-check → **NO-GO (M4 dropped)**
Fold 0 only (`--folds 0`, registry status `partial_folds`), 6 threads, torch 2.6.0+cpu.

| arm | run id | fold-0 P1 variable | P1 all | P2b nAP | P2b AUROC | inner-val P1 | chosen | fit s |
|---|---|---|---|---|---|---|---|---|
| RidgePC (256 PCs) | 20261026-231523 | **0.3068** | 0.1531 | 0.1889 | 0.6634 | 0.2912 | no whitening, λ/n = 0.032 | 40 |
| MLP seed 0 | 20261026-231657 | 0.3025 | 0.1581 | 0.1832 | 0.6581 | 0.2868 | dropout 0.5, no whitening, 34 ep | 274 |
| MLP seed 1 | 20261026-232221 | 0.3028 | 0.1616 | 0.1836 | 0.6590 | 0.2872 | dropout 0.2, 17 ep | 268 |
| MLP seed 2 | 20261026-232818 | 0.3033 | 0.1617 | 0.1823 | 0.6585 | 0.2889 | dropout 0.5, 23 ep | 348 |
| MLP seed 3 | 20261026-233541 | 0.3017 | 0.1597 | 0.1809 | 0.6574 | 0.2888 | dropout 0.2, 17 ep | 316 |
| MLP seed 4 | 20261026-234238 | 0.3038 | 0.1611 | 0.1836 | 0.6592 | 0.2900 | dropout 0.5, 28 ep | 335 |
| **MLP mean (5 seeds)** | | **0.3028** (seed SD 0.0008) | 0.1604 | 0.1827 | 0.6584 | | | |
| context: m1_E_resB0 fold 0 | 20261026-222117 | 0.3129 | | | | | | |

- Pre-registered rule: GO iff mean MLP fold-0 P1 > RidgePC fold-0 P1. 0.3028 < 0.3068
  (Δ −0.0040; every seed loses by 0.0030–0.0051). **M4 is dropped and was not built.**
- The MLP also lost on the inner-validation criterion for every seed (0.287–0.290 vs 0.291), so
  the outcome does not hinge on the outer fold. Whitened PCs were clearly worse for the MLP
  (inner 0.256–0.264) and were never chosen; for ridge, whitening also lost.
- The MLP is better on P1 *all* genes (0.160 vs 0.153), i.e. on low-variance genes, but worse on
  the selection metric and on P2b. Not a selection criterion; noted only.
- Ridge on 256 training-fold PCs is itself 0.006 below the full-kernel m1_E_resB0 on fold 0, so the
  PC bottleneck costs accuracy before any nonlinearity is added.

### R2. M3 coverage decision (pre-registered rule)
- Timing probe `g3_m3_probe` (run 20261026-234919, fold 0, first 100 genes of the subset; 99
  were training-fold variable): base (5 inner M1 fits + full M1) 254 s, screening 4 s, trees 97 s
  for 99 genes at 6 processes. Projection for all ~2,120 training-fold variable genes × 5 folds:
  ≈ 5 × (254 + 2,070 + 90) s ≈ 3.3 h ≤ 4 h → **full coverage** (`g3_m3`).
- Probe numbers (fold 0, 100-gene subset; not a selection result): P1 0.3417 vs m1_E_resB0 0.3344
  on the same genes (+0.0073); 69% of genes got no correction (early stopping found no
  improvement over a zero correction).

### R3. R2/R4 coverage decisions
- R2: see deviation D-G3-3 → `g3_r2_sub150`.
- R4 timing probe `g3_r4_probe` (run 20261027-004320, fold 0, 29 modelled genes): 55 s at 6
  processes ≈ 1.9 s per gene → full ≈ 5.6 h > 2 h → **pre-declared 500-gene subset**
  (`g3_r4_sub500`, nested) and the leaky non-nested run on the same genes (`g3_r4_oracle_sub500`).

### R4. M3 LightGBM on m1_E_resB0 residuals → **NO-GO (does not replace m1_E_resB0)**
Run `g3_m3n` 20261027-014824 (nested early stopping, D-G3-1), all training-fold variable genes
(2,143–2,179 modelled per fold), 5 folds, 6 processes × 1 thread for the trees, 6 BLAS threads
for the M1 base. LightGBM 4.7.0.

| class / metric | M3 (g3_m3n) | ΔM3 − m1_E_resB0 (SE; folds better) |
|---|---|---|
| P1 variable (selection) | 0.3099 ± 0.0046 | **+0.0018 (0.0002; 5/5)** |
| P1 selective | 0.2907 | +0.0026 (0.0001; 5/5) |
| P1 all | 0.1549 | +0.0002 (0.0000) |
| P1 rare-selective | 0.1725 | +0.0001 |
| P1 non-essential (NEGv1) | 0.1142 | +0.0001 |
| P2b nAP (selective, prev ≤ 0.5) | 0.1931 | +0.0023 (0.0002) |
| P2b AUROC | 0.6614 | +0.0007 (0.0001) |

- Selection rule: ΔP1 = +0.0018 > 1 SE (0.0002) but **< δ = 0.005** → M3 does not replace
  m1_E_resB0. The gain is real (5/5 folds) but small: per gene-fold the change is ±0.04 (5–95%),
  improved in 37% of gene-folds, unchanged (no trees) in 28% of genes; median 4 boosting rounds.
  The largest per-gene gains (r +0.2 to +0.47) are in a handful of genes (e.g. Entrez 4893 NRAS,
  9130, 131118, 474170, 1964).
- Runtime: 2,087 s fit+predict for 5 folds (wall ≈ 39 min, including 3 folds whose inner M1 OOF
  was not cached: ≈ 500 s each for 10 M1 fits); the trees themselves take ≈ 47 s per fold for
  ~2,170 genes. Peak RSS 4.2 GB (above the DESIGN §10 per-stage guide of 4 GB).
- Interpretation: once expression is used linearly, the residual nonlinearity that boosted
  trees on screened E/C/M/H/lineage features can find is worth ~0.002 P1 — consistent with the
  literature (RF ≈ EN). Not worth the extra rung under the pre-registered δ.
