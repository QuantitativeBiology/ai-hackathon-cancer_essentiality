# G3 nonlinear track — M4 pre-check / M4, M3 LightGBM, R2/R4, published-comparable block

Status: **complete** (2026-10-27). The pre-registration below was committed before any run of
this track and is not edited; changes are logged in "Deviations"; results follow.

**Outcome:** M4 NO-GO (MLP loses to ridge on the same PCs); M3 NO-GO (ΔP1 +0.0018, SE 0.0002,
5/5 folds, < δ; un-nested sensitivity +0.0043); R2 elastic net −0.030 and R4 top-k −0.038 vs
m1_E_resB0 on the same random genes (R4 non-nested optimism +0.018). `m1_E_resB0` remains the
headline model.

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
  The largest per-gene gains (mean r +0.19 to +0.47) are a handful of genes with known
  threshold-type dependencies: NRAS (hotspot-mutation addiction), and X/Y or paralog-buffered
  genes (FAM50A, EIF1AX, RPS4X, NAA10, RPP25L, DNAJC19, LRRC37A2) whose dependency switches on
  when the paralog is lost or unexpressed — exactly where trees beat a linear kernel.
- Runtime: 2,087 s fit+predict for 5 folds (wall ≈ 39 min, including 3 folds whose inner M1 OOF
  was not cached: ≈ 500 s each for 10 M1 fits); the trees themselves take ≈ 47 s per fold for
  ~2,170 genes. Peak RSS 4.2 GB (above the DESIGN §10 per-stage guide of 4 GB).
- Interpretation: once expression is used linearly, the residual nonlinearity that boosted
  trees on screened E/C/M/H/lineage features can find is worth ~0.002 P1 — consistent with the
  literature (RF ≈ EN). Not worth the extra rung under the pre-registered δ.

### R5. R2 Shi-style elastic net (published-method re-run; not a selection candidate)
Run `g3_r2_sub150` 20261027-022710: first 150 genes of the pre-declared random variable subset
(D-G3-3), all 5 folds, exact coordinate descent (sklearn 1.5.2 `ElasticNetCV`, l1_ratio 0.5,
20 penalties to 0.05 × max, 3 patient-grouped inner folds), ~18.1k standardised E features.
Wall 2 h 08 min (fit+predict 7,485 s) at 6 processes — the 2 h box was met only by cutting to
150 genes.

Same 150 genes (dev-wide variable, evaluable in every fold), per-fold means then averaged:

| model | P1 | ΔP1 vs m1_E_resB0 (SE; folds) | mean Spearman | genes r ≥ 0.2 (per fold) |
|---|---|---|---|---|
| m1_E_resB0 | 0.3270 | — | 0.2944 | 122.4 |
| R1 (ridge) | 0.3263 | −0.0007 (0.0007) | 0.2931 | 123.2 |
| **R2 elastic net** | **0.2967** | **−0.0303 (0.0014; 0/5)** | 0.2600 | 105.4 |
| B2 (kNN) | 0.2726 | −0.0545 | 0.2411 | 101.2 |
| B1 (lineage mean) | 0.2040 | −0.1230 | 0.1818 | 72.6 |

- Elastic net keeps a median of ~50 features per gene and loses 0.030 P1 to dense multi-output
  ridge (better on only 30% of gene-folds). This matches the literature's "ridge ≥ sparse linear"
  on held-out lines, and it is ~1,000× more expensive per gene than the kernel ridge.
- P1 all / P2b of this run are not meaningful (all other genes are predicted by the training mean).

### R6. R4 TrPLet-like top-k: nested vs non-nested (published-method re-run; not a selection candidate)
Runs `g3_r4_sub500` 20261027-043531 (nested) and `g3_r4_oracle_sub500` 20261027-052843 (leaky=1),
the pre-declared 500-gene random variable subset, all 5 folds, 6 processes. 12 candidates per gene:
top-k E features by |r| (k ∈ 5…200) × {ridge (LOO α), RBF-SVR}. Wall: nested 53 min, oracle 22 min.

Same 500 genes, per-fold means then averaged:

| model | P1 ± fold SD | ΔP1 vs m1_E_resB0 (SE; folds) | mean Spearman | genes r ≥ 0.2 (per fold) |
|---|---|---|---|---|
| m1_E_resB0 | 0.3182 ± 0.0054 | — | 0.2915 | 410.0 |
| R1 (ridge) | 0.3172 ± 0.0062 | −0.0010 (0.0008) | 0.2897 | 409.8 |
| **R4 nested** | **0.2800 ± 0.0065** | **−0.0382 (0.0016; 0/5)** | 0.2487 | 353.8 |
| R4 non-nested (leaky, best-of-12 by outer-CV r) | 0.2981 ± 0.0065 | −0.0201 (0.0014; 0/5) | 0.2651 | 381.4 |
| B2 (kNN) | 0.2615 | −0.0567 | 0.2345 | 330.8 |
| B1 (lineage mean) | 0.1915 | −0.1267 | 0.1729 | 218.4 |

- **Optimism of non-nested selection: +0.0181 P1 (SE 0.0003, 5/5 folds)** — about half of R4's gap to
  ridge, from choosing 1 of only 12 candidates per gene on the reported folds (TrPLet chose 1 of 60).
  Even with this optimism R4 stays 0.020 below multi-output ridge.
- Nested choices favour the largest k (ridge k=200 33%, SVR k=200 23%; 2,500 gene-folds); the
  oracle also picks k=200 most often. k=200 is the grid edge, so R4 with larger k might close part
  of the gap (TrPLet used larger sets); not run (time box), noted as a limitation.
- P1 all / P2b of these runs are not meaningful (non-modelled genes get the training mean).

### R7. M3 sensitivity: the originally pre-registered (un-nested) early stopping, and controls
After D-G3-1 the un-nested config `g3_m3` was re-run unchanged (20261027-060743; the killed row
20261026-235648 has no numbers) as a **sensitivity analysis**, decided after the nested result was
known. It is reported, not selected on.

| run | ΔP1 vs m1_E_resB0 (SE; folds) | P1 variable | genes without correction | median rounds |
|---|---|---|---|---|
| `g3_m3n` (nested ES; the M3 of record) | +0.0018 (0.0002; 5/5) | 0.3099 | 28% | 4 |
| `g3_m3` (un-nested ES; sensitivity) | **+0.0043 (0.0001; 5/5)** | 0.3124 | 69% | 0 |
| `g3_m3n_perm` (control: pool rows permuted among training lines, no M1 feature) | −0.0042 (0.0001; 0/5) | 0.3039 | 45% | 1 |
| `g3_m3_perm` (same control, un-nested ES) | −0.0037 (0.0001; 0/5) | 0.3044 | — | — |

- Both versions stay **below δ = 0.005**, so the M3 decision (NO-GO) does not depend on the
  early-stopping variant.
- The un-nested variant scores higher because its per-gene gate is stricter: 69% of genes get no
  trees at all, versus 28% with nested ES. The permuted controls show what weak corrections cost
  when there is no signal (−0.004), so the nested gate admits too many small, noisy corrections.
  The reviewer's toy prediction ("un-nested is pessimistic") holds for the inner criterion but not
  for outer P1. A stricter gate (e.g. a minimum validation-MSE improvement) is the obvious next
  idea, but it would be chosen after seeing these results, and the ceiling suggested here
  (about +0.004 to +0.005) is at δ at best.
- Controls pass: with no line information the trees cannot help (ΔP1 < 0 in every fold).

### R8. Published-comparable block (**loosely comparable**; not used for selection)
DeepDEP 1,298 genes: 1,286 in our evaluable set; LEAP 1,539: 1,538. Per-fold = metric within each
held-out fold (~204 lines), averaged over folds (D3); pooled = all 1,022 out-of-fold dev
predictions at once. Full-coverage runs only (subset runs predict the training mean elsewhere).

| model | Spearman DeepDEP (per-fold / pooled) | Spearman LEAP (per-fold / pooled) | genes r ≥ 0.2 of 17,915 (per-fold / pooled) | Shi count r > 0.2 (per-fold of 4,028 eligible / pooled of 9,028) |
|---|---|---|---|---|
| B1 lineage mean | 0.173 / 0.160 | 0.178 / 0.165 | 2,203 / 1,207 | 1,235 / 1,112 |
| B2 kNN | 0.233 / 0.229 | 0.245 / 0.240 | 3,858 / 2,828 | 2,010 / 2,467 |
| R1 ridge | 0.282 / 0.278 | 0.301 / 0.297 | 5,511 / 4,367 | 2,630 / 3,583 |
| **m1_E_resB0** | 0.284 / 0.280 | 0.303 / 0.300 | 5,507 / 4,362 | 2,631 / 3,578 |
| M3 nested (`g3_m3n`) | 0.285 / 0.281 | 0.304 / 0.301 | 5,482 / 4,333 | 2,619 / 3,561 |
| M3 un-nested (`g3_m3`) | 0.286 / 0.283 | 0.306 / 0.302 | 5,505 / 4,366 | 2,630 / 3,584 |

Subset-only re-runs (compared with m1_E_resB0 on the same genes; see R5/R6): R2 mean Spearman
0.260 vs 0.294 (150 genes); R4 nested 0.249 and non-nested 0.265 vs 0.292 (500 genes).

How loosely comparable:
- Published numbers: DeepDEP mean per-gene ρ 0.16–0.18 (ridge 0.276 in Chang & Zhang's
  re-analysis) on 2018 CERES with 278 lines; LEAP PS-LASSO 0.321 (23Q4, 1,019 lines, 10 × 80/20
  splits; its README asks for CRISPRGeneDependency, so its target is probably the dependency
  probability, not gene effect); TrPLet 6,283 of 16,845 genes at r ≥ 0.2 (23Q2, non-nested
  best-of-60); Shi 1,966 of 7,260 (1Q21 CERES, genes with ≥ 5 dependent lines over all data).
- Our per-fold counts at n ≈ 204 lines are inflated by sampling noise (a gene with true r ≈ 0.1
  often crosses 0.2); the pooled counts (n = 1,022) are the closer analogue of the published
  pooled-OOF counts. Ours are patient-grouped, lineage-stratified, Chronos 26Q3, without
  non-nested selection.
- Reading: the linear expression model sits near LEAP's best on LEAP's genes (0.30 vs 0.32) and well
  above DeepDEP's published deep model on DeepDEP's genes, consistent with Chang & Zhang. Our
  TrPLet-style pooled count (4,362) is below TrPLet's 6,283; R6 shows that non-nested best-of-N
  selection alone adds about +0.018 r per gene.

### R9. Summary table (all G3 runs; generated by `scripts/g3_report.py`)
P1 = mean per-gene Pearson r on dev-wide variable genes, per fold then averaged (± fold SD); Δ is
paired by fold against m1_E_resB0 on the same genes (SE = SD/√n folds); fit+pred = sum over folds
of fit and predict seconds on a shared machine (indicative). `results/g3_nonlinear_summary.json`
holds the same numbers.

| run | run id | genes | folds | P1 var ± SD | ΔP1 vs m1_E_resB0 ± SE (folds better) | P1 all | P2b nAP | P2b AUROC | fit+pred s | replaces |
|---|---|---|---|---|---|---|---|---|---|---|
| b1 | 20261026-220921 | all variable | 5 | 0.1837 ± 0.0131 | -0.1244 ± 0.0048 (0/5) | 0.0928 | 0.0867 | 0.5941 | 105 | no |
| b2 | 20261026-221456 | all variable | 5 | 0.2527 ± 0.0087 | -0.0554 ± 0.0023 (0/5) | 0.1205 | 0.1524 | 0.6279 | 129 | no |
| r1_E | 20261026-223037 | all variable | 5 | 0.3073 ± 0.0050 | -0.0008 ± 0.0008 (1/5) | 0.1546 | 0.1898 | 0.6603 | 360 | no |
| m1_E_resB0 | 20261026-222117 | all variable | 5 | 0.3081 ± 0.0047 | +0.0000 ± 0.0000 (0/5) | 0.1547 | 0.1908 | 0.6607 | 331 | no |
| g3_m4pre_ridgepc | 20261026-231523 | all variable | 1 | 0.3068 ± — | -0.0061 ± — (0/1) | 0.1531 | 0.1889 | 0.6634 | 40 | no |
| g3_m4pre_mlp_s0 | 20261026-231657 | all variable | 1 | 0.3025 ± — | -0.0105 ± — (0/1) | 0.1581 | 0.1832 | 0.6581 | 275 | no |
| g3_m4pre_mlp_s1 | 20261026-232221 | all variable | 1 | 0.3028 ± — | -0.0102 ± — (0/1) | 0.1616 | 0.1836 | 0.6590 | 269 | no |
| g3_m4pre_mlp_s2 | 20261026-232818 | all variable | 1 | 0.3033 ± — | -0.0097 ± — (0/1) | 0.1617 | 0.1823 | 0.6585 | 349 | no |
| g3_m4pre_mlp_s3 | 20261026-233541 | all variable | 1 | 0.3017 ± — | -0.0113 ± — (0/1) | 0.1597 | 0.1809 | 0.6574 | 317 | no |
| g3_m4pre_mlp_s4 | 20261026-234238 | all variable | 1 | 0.3038 ± — | -0.0092 ± — (0/1) | 0.1611 | 0.1836 | 0.6592 | 335 | no |
| g3_m3_probe | 20261026-234919 | subset100 | 1 | 0.3417 ± — | +0.0073 ± — (1/1) | 0.1550 † | 0.1939 † | 0.6668 † | 372 | no |
| g3_m3n | 20261027-014824 | all variable | 5 | 0.3099 ± 0.0046 | +0.0018 ± 0.0002 (5/5) | 0.1549 | 0.1931 | 0.6614 | 2087 | no |
| g3_m3n_perm | 20261027-055023 | all variable | 5 | 0.3039 ± 0.0045 | -0.0042 ± 0.0001 (0/5) | 0.1542 | 0.1895 | 0.6599 | 793 | no |
| g3_m3 | 20261027-060743 | all variable | 5 | 0.3124 ± 0.0046 | +0.0043 ± 0.0001 (5/5) | 0.1552 | 0.1936 | 0.6618 | 772 | no |
| g3_m3_perm | 20261027-062459 | all variable | 5 | 0.3044 ± 0.0045 | -0.0037 ± 0.0001 (0/5) | 0.1543 | 0.1898 | 0.6601 | 793 | no |
| g3_r2_sub150 | 20261027-022710 | subset150 | 5 | 0.2967 ± 0.0095 | -0.0303 ± 0.0014 (0/5) | 0.0025 † | 0.0082 † | 0.5062 † | 7485 | no |
| g3_r4_probe | 20261027-004320 | subset30 | 1 | 0.3239 ± — | -0.0328 ± — (0/1) | 0.0005 † | 0.0021 † | 0.5014 † | 60 | no |
| g3_r4_sub500 | 20261027-043531 | subset500 | 5 | 0.2800 ± 0.0065 | -0.0382 ± 0.0016 (0/5) | 0.0078 † | 0.0216 † | 0.5152 † | 2973 | no |
| g3_r4_oracle_sub500 (leaky) | 20261027-052843 | subset500 | 5 | 0.2981 ± 0.0065 | -0.0201 ± 0.0014 (0/5) | 0.0083 † | 0.0228 † | 0.5159 † | 1101 | no |

† subset run: P1 all and P2b cover every gene, most of them passed through (M3: M1; R2/R4: training mean); ΔP1 and 'replaces' refer to the subset genes only.

### R10. Runtimes (indicative; shared machine, 6 of 12 cores) and versions
| run | wall-clock | notes |
|---|---|---|
| M4 pre-check, ridge on PCs | 1.5 min (fold 0) | 6 BLAS threads |
| M4 pre-check, MLP | 4.5–5.8 min per seed (fold 0) | torch 6 threads; 4 grid points with early stopping + refit |
| M3 `g3_m3n` | 39 min (5 folds) | about 8 min per fold for 10 uncached inner M1 fits; trees about 47 s per fold for ~2,170 genes |
| M3 `g3_m3`, both controls | 13–17 min each | inner M1 OOF cached |
| R2 `g3_r2_sub150` | 2 h 08 min | 150 genes; about 45 s per gene per process |
| R4 `g3_r4_sub500` / oracle | 53 / 22 min | 500 genes |
| (reference) m1_E_resB0 | about 6 min | registry, 6 threads |

Peak RSS: M3 4.2–4.3 GB (above the DESIGN §10 guide of 4 GB per stage); others ≤ 2.6 GB.
Versions: Python 3.13.1, numpy 2.1.3, scipy 1.15.2, pandas 2.2.3, scikit-learn 1.5.2, torch
2.6.0+cpu, lightgbm 4.7.0, threadpoolctl 3.6.0. py-spy 0.4.2 was installed (user site) to
diagnose a stalled run.

Provenance: every G3 row has `tree_clean = 0` (the other track's untracked/modified files;
accepted for G3). `dirty_tree = 1` for MLP seeds 0–2 (an uncommitted refactor of a
`g3_common.py` helper the MLP does not use); all other G3 rows have `dirty_tree = 0`.
Killed/failed rows, none with numbers: `g3_m3` 20261026-235648 (D-G3-1), `g3_m3n`
20261027-004513 (main-process bottleneck in building tree inputs; fixed in 3910d77),
`g3_r2_probe` 20261027-000958 (time box), `g3_r4_probe` 20261027-003100 (queue re-plan).

## Interpretation and go/no-go
- **M4 (neural): NO-GO** at the pre-registered pre-check. An MLP on 256 training-fold PCs loses to
  ridge on the same PCs (fold-0 P1 0.3028 vs 0.3068; all 5 seeds lose), and ridge on PCs is itself
  0.006 below the full kernel. M4 was not built.
- **M3 (LightGBM on m1_E_resB0 residuals): NO-GO.** The gain is real but small: +0.0018 (nested,
  the M3 of record) or +0.0043 (un-nested sensitivity), 5/5 folds, both below δ = 0.005. It sits in
  a few threshold-type dependencies (NRAS hotspot; X/Y- or paralog-buffered genes such as RPS4X,
  EIF1AX, FAM50A, NAA10). Permuted controls are negative.
- **Published methods:** elastic net (R2) −0.030 and TrPLet-like top-k (R4) −0.038 vs m1_E_resB0
  on the same random genes. R4's non-nested selection inflates P1 by +0.018 and still trails ridge
  by 0.020. Dense multi-output ridge on expression remains the best model; nothing in the nonlinear
  track replaces `m1_E_resB0`.
- This agrees with the literature (linear models match or beat deep and tree models on held-out
  lines). The residual nonlinear signal is concentrated in a few genes. If it is wanted, a targeted
  per-gene stage for paralog/sex-chromosome and hotspot genes (cf. the M2 priors) is a cheaper route
  than a general rung.

## Open issues for the PI
1. M3 early-stopping variant: the nested fix (D-G3-1, decided before results) gave the lower
   number; the un-nested original gives +0.0043. Neither changes the decision, but if M5 stacking
   is pursued, the un-nested M3 is the better ingredient, and choosing it now would be post hoc.
2. R2 ran on 150 genes, not the pre-declared 500 (D-G3-3), because exact elastic net on ~18k
   features costs about 45 s per gene. R4's k grid hit its upper edge (k = 200).
3. LEAP's published target is probably the dependency probability, not gene effect; the LEAP
   comparison is loose.
4. M3 peak memory (4.2 GB) exceeds the 4 GB per-stage guide.
