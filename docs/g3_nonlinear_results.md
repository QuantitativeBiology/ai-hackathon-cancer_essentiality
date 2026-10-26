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
(none yet)
