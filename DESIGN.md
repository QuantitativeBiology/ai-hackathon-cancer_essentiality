# DESIGN — predicting CRISPR gene essentiality for unscreened human cells

Version: v0 (draft for review round 1), 2026-10-26. Goal and parts are fixed by PROMPT.md. Evidence base: LITERATURE.md.

## 1. Questions
- Q1. Given a never-screened cancer cell line's molecular profile, how well can we predict its Chronos gene effect
  profile, and which feature types matter (ablation)?
- Q2. How much better than per-gene mean and per-gene lineage-mean baselines, at what runtime?
- Q3. Does the model transfer to non-transformed human cells (hTERT, ESC/iPSC and derived, primary haematopoietic)?
  Separately: (a) core-essential genes; (b) context-specific dependencies. Answered only after screen batch effects
  are harmonised and the harmonisation itself is validated.

## 2. Data
- **Targets:** DepMap 26Q3 `CRISPRGeneEffect` (Chronos), 1361 models × 18434 genes; 3.2% NaN.
  - NaNs are never imputed in evaluation; metrics use observed entries only.
  - For training, NaN targets are filled with the training-fold gene mean (ridge needs a complete Y). A sensitivity
    check drops genes with >10% NaN.
- **Inputs (DepMap, local):**
  - E: expression, log2(TPM+1), protein-coding.
  - C: gene-level relative CN; log2(CN+0.01) transform; missing → 0 + missing indicator.
  - M: damaging mutations (0/1/2 → binary).
  - H: hotspot mutations (binary).
  - F: fusions (gene-partner indicators; partners present in >20% of models dropped as artefact-like).
  - Meta: lineage, primary disease, sex, age, PrimaryOrMetastasis, ModelType (cell line/organoid), GrowthPattern.
  - Excluded as screen-status proxies (audit flag): media, FormulationID, CatalogNumber, PublicComments,
    DataSource, AnchorDrug, PassageNumber, EngineeredModel, CulturedResistanceDrug, ModelTreatment.
- **Not available locally:** signatures, methylation, proteomics, metabolomics. They need a manual portal download.
  If the user supplies them, they enter the ablation on their line subsets. Open alternative to try: Sanger Cell
  Model Passports proteomics (matched through SangerModelID). Their absence is reported as a limitation, not
  silently dropped.
- **Prior knowledge (allowed, non-text-mined, non-CRISPR):**
  - Ensembl 116 paralogs, with % identity.
  - BioPlex 3.0, HuRI, BioGRID physical (low-throughput + high-throughput, excluding ORCS).
  - Reactome, MSigDB Hallmarks, HGNC gene groups.
  - CORUM if obtainable.
  - Gene embeddings built from training-fold CRISPR profiles only (no test lines).
- **Derived line features:** mutation burden; CN-based genome-instability score (fraction of genes with |log2CN|>0.3);
  TP53 functional score (TP53 damaging/hotspot status + expression of a fixed p53 target set: CDKN1A, MDM2, RRM2B,
  ZMAT3, TNFRSF10B, FDXR, SESN1, DDB2); Hallmark expression scores.

## 3. Populations and splits (locked before modelling; `configs/splits.csv`, seed 20261026)
- **Cancer population:** screened models with expression, excluding Non-Cancerous (1277).
  - Locked **test** set: 255 lines. Grouped by PatientID, stratified by lineage. Touched once, at the end, for the
    pre-declared finalists.
  - **Dev** set: 1022 lines, 5 patient-grouped, lineage-stratified CV folds. All model selection, tuning and
    ablation happen here, with nested inner splits for hyper-parameters.
- **Non-transformed set A:** 15 DepMap Non-Cancerous models with RNA-seq (same library and scoring; no screen batch).
  Patients overlapping them are excluded from the cancer sets.
- **Non-transformed set B:** external published screens (§7).
- **Secondary split:** leave-one-lineage-out on dev, measuring transfer to an unseen lineage (closest analogue to
  non-transformed transfer).

## 4. Evaluation (pre-registered)
### Gene classes (defined on dev-training data only)
- *all*: all genes.
- *variable*: SD of gene effect > 0.25.
- *selective*: 5%–95% of dev lines with effect < −0.5.
- *common essential*: effect < −0.5 in >95% of dev lines.
- *non-essential*: Hart NEGv1.

### Metrics
- **Primary metric P1 (model selection):** mean per-gene Pearson r across held-out lines on *variable* genes.
  - Constant predictions get r = 0, so the gene-mean baseline scores 0 by construction. This is correct: it carries
    no information across lines.
  - Also reported: median; mean on *all* and *selective*; number of genes with r > 0.3; fraction of genes beating
    the lineage baseline.
- **Primary metric P2:** AUPRC for essentiality, label = observed effect < −0.5. Two views:
  - P2a, per line: AUPRC across all genes, averaged over held-out lines. This mostly measures core essentials.
  - P2b, per gene across held-out lines, for *selective* genes with ≥3 positives and ≥3 negatives in the evaluated
    lines. Reported as mean AUPRC and mean (AUPRC − prevalence).
  - Sensitivity check: the −1.0 threshold.
- **Secondary metrics:** per-line Pearson/Spearman of the predicted profile, RMSE, per-gene Spearman, calibration
  of predicted effect.
- **Uncertainty:**
  - CV: mean ± SD over folds and over the 5 model seeds (stochastic models).
  - Test: 95% bootstrap CIs resampling patients (1000 resamples).
  - Model-vs-baseline comparisons: paired gene-level differences, with paired bootstrap over patients.
- **Runtime:** wall-clock fit and predict on 12 threads, peak RSS, logged by a common harness.

### Selection rule
The finalist is the configuration with the best dev-CV P1. A simpler model is preferred when within 1 SE (fold SD/√5).
P2b is reported alongside but does not drive selection. At most 3 finalists, plus baselines, run on test, once.

## 5. Models (each must beat the previous rung in dev CV to justify its cost)
- **B0** per-gene mean.
- **B1** per-gene lineage mean, shrunk toward the gene mean (pseudo-count κ tuned in CV); unseen lineage → B0.
- **B2** kNN in expression space (correlation on top-variance genes; k tuned).
- **M1** multi-output kernel ridge (linear kernel) on standardised features, predicting the residual from B1. λ per
  gene chosen by inner-CV or GCV. Seconds per fold.
- **M2** M1 + per-gene local features in a per-gene second-stage ridge on M1 residuals:
  - own expression, CN and mutation;
  - expression of top-3 paralogs (by identity) and of the paralog with minimum expression;
  - mean expression of complex/interaction partners.
- **M3** LightGBM, per gene, for *variable* genes only (inner-fold feature screening: local features + top
  correlated expression features + PCs). For other genes, M2 predictions are used.
- **M4** shared-parameter neural model (CPU). Predicts y[c,g] from a cell embedding (from features), a gene
  embedding (from training-fold Y and prior knowledge) and per-(c,g) local features. Trained across all genes, so
  sparse genes borrow strength.
- **M5** ensemble (stacked or averaged) of the best of M1–M4, with weights fit on inner folds.
- Ideas considered and not prioritised: foundation-model gene/cell embeddings (beaten by raw expression in 2025
  benchmarks; GPU-heavy); graph networks over STRING (text-mined, excluded).

## 6. Ablation
- On the best linear model (M1/M2) and the final model:
  - add-one modality to E;
  - leave-one-out from the full set {E, C, M, H, F, Meta, derived, prior knowledge};
  - E-free variants (genomics only), to quantify how far the model rests on expression.
- Prior-knowledge ablation in M2/M4: none / paralogs / +complexes / +PPI.
- Line-subset analysis: models with all of E, C, M (the 1214 intersection). Modality comparisons use the same lines
  so that availability does not bias them.
- Same dev folds, same metrics; report P1, P2b and runtime per configuration.

## 7. Non-transformed generalisation
### 7.1 Assembly
- Download the least processed data per screen (counts > LFC > BF).
- Matched expression: same lab or same cell line (ENCODE, HipSci, GEO; DepMap for RPE1).
- CRISPRi screens (iPSC neurons and microglia) form a separate stratum.
- Inventory: `docs/external_screens_inventory.md`.

### 7.2 Screen harmonisation (validated before use)
1. **Uniform re-scoring.** When counts exist, compute guide LFC vs plasmid/T0, then fit Chronos fresh per library
   when feasible (otherwise median guide LFC per gene).
2. **Affine scaling** per screen: median of a *scaling* control set → −1, median of non-essentials → 0.
   - Scaling set: DepMap common essentials defined on dev lines, minus a held-out half of CEGv2.
   - Non-essentials: half of NEGv1.
   - The other halves of CEGv2/NEGv1 are kept for evaluation only.
3. **Gene-wise location/scale correction** (empirical-Bayes, ComBat-like) per lab/library, fitted on anchor lines
   screened both externally and in DepMap (HAP1, HCT116, DLD-1, A-375, K-562, HT-29, AML lines, RPE1). Applied
   only if leave-one-anchor-out validation improves agreement with DepMap.
4. **Duration.** Covered by per-screen affine scaling. A per-gene depletion-kinetics term is tested only if anchors
   show residual time-dependent bias.
5. **Validation checks:**
   - held-out CEG/NEG separation (AUROC);
   - NEG median ≈ 0;
   - anchor same-line correlation with DepMap, before and after each step;
   - cross-study agreement for the same cell type (e.g. three RPE1 screens; DepMap RPE1-SS vs external RPE1).
   - This agreement is the empirical **ceiling** for any model.

### 7.3 Expression harmonisation
- Recompute log2(TPM+1) on DepMap's protein-coding gene set. Quantile-normalise to the DepMap expression
  distribution.
- Validate on ENCODE/GEO profiles of DepMap lines (K-562, HepG2, MCF-7, A-549, HCT116): same-line nearest
  neighbour in DepMap, and correlation with the DepMap profile.

### 7.4 Questions and tests
- **Core essentials transfer:** per screen, AUPRC/AUROC of predicted effects for held-out CEGv2 vs NEGv1. Also
  AUPRC for harmonised effect < −0.5 across all genes; model vs gene-mean baseline vs observed ceiling.
- **Context-specific dependencies:**
  - Per screen, correlation across *selective* genes between predicted and observed deviations from the
    pan-cancer gene mean (pred − mean vs obs − mean).
  - Where ≥4 non-transformed samples share a pipeline (DepMap set A), per-gene correlation across them.
  - Pre-specified hypotheses:
    - TP53-WT backgrounds: positive TP53/CDKN1A/USP28 KO effects and MDM2/PPM1D dependency, contrasted with the
      TP53-null RPE1 pairs.
    - Lineage-like dependencies: pluripotency and haematopoietic regulators where measured.
    - Loss of cancer-specific dependencies (oncogene addictions such as KRAS/BRAF/CTNNB1 predicted near 0 in
      untransformed cells).
- **Failure attribution:** a transfer gap counts as model failure only if it remains after harmonisation, exceeds
  the cross-study ceiling, and is absent for anchor lines processed the same way.

## 8. Verification and controls
- Unit tests with toy data and known answers:
  - metric functions (perfect → r = 1, AUPRC = 1; random → AUPRC ≈ prevalence; constant → r = 0);
  - split integrity (no patient across roles);
  - default-entry loading.
- **Leakage tests:** feature scaling, feature selection, gene classes and lineage means are computed on training
  folds only. An automated assertion checks that no test/non-cancerous ModelIDs appear in any training matrix.
- **Positive control:** synthetic targets planted as a known function of features (own-gene expression + noise).
  The pipeline must recover the planted r.
- **Negative controls:**
  - permuted line labels (features shuffled across lines within the training fold) → P1 ≈ 0;
  - random-feature models.
- **Harmonisation controls:** §7.2 step 5; a simulated batch (known affine + gene-wise distortion of DepMap anchor
  profiles) must be undone.
- **Independent review:** every headline number is re-checked by an agent that did not produce it, from the
  committed scripts and outputs.

## 9. Reproducibility and record keeping
- Scripts only (`scripts/`, `src/essml/`), configs in `configs/`. Seeds in `configs/seeds.json`.
- An experiment registry (`results/registry.csv`) lists config, commit, metrics and runtime for every run,
  including failed or negative ones.
- `DATA.md` lists source, version, checksum and licence. No raw or restricted data in git.
- `RESEARCH_TREE.md` holds decisions; `LESSONS.md` holds process rules.

## 10. Compute plan (12 cores, 16 GB, CPU)
- M1: seconds per fold.
- M2: minutes (≈18k small per-gene ridges).
- M3: ≈2k genes × 5 folds of LightGBM; estimated ≤2 h with early stopping, parallel over genes.
- M4: ≤2 h per CV run on CPU (minibatch over (cell, gene) pairs).
- External harmonisation: Chronos fits ≤ a few hours each.
- Memory: Y and X in float32 ≈ 80 MB each; per-fold copies are fine.

## 11. Risks
- External data availability (iPSC-derived and primary categories are thin).
- Expression-pipeline mismatch.
- CRISPRi vs KO incomparability.
- TP53-WT DNA-damage response in untransformed cells depresses all cutting screens.
- Small n in set A.
- Ceiling effects from short primary-cell screens.
