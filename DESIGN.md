# DESIGN — predicting CRISPR gene essentiality for unscreened human cells

Version: **v3 — FROZEN** (after review rounds 1–3; responses in `reviews/round{1,2,3}_response.md`), 2026-10-26.
Changes after freeze are logged in RESEARCH_TREE.md with the reason.
Goal and parts are fixed by PROMPT.md. Evidence base: LITERATURE.md.

## 1. Questions
- **Q1.** Given a never-screened cancer cell line's molecular profile, how well can we predict its Chronos gene effect
  profile? Which feature types matter (ablation)?
- **Q2.** How much better is that than the per-gene mean (B0) and the per-gene lineage mean (B1), at what runtime?
- **Q3.** Does the model transfer to non-transformed human cells? Asked separately for (a) core-essential genes and
  (b) context-specific dependencies. Answered only after screen batch effects are harmonised, the harmonisation is
  validated on controls and anchors, and transfer criteria are committed in advance (§7.6).

## 2. Data and features
### Targets
DepMap 26Q3 `CRISPRGeneEffect` (Chronos), 1361 × 18434.
- **NaNs are library-structured:** 49 gene-level NaN patterns and 4 major line groups (827/221/118/195 lines;
  R1 data review). NaNs are never filled.
  - Ridge models are fitted exactly on observed rows, per NaN pattern (one kernel solve per pattern).
  - Neural and stacking models use masked losses.
  - A gene is evaluated in a held-out set only if it has ≥50 observed held-out lines (CV folds: ≥40).
  - The NaN line group is kept as a *library proxy* for stratified reporting. CRISPRScreenMap is not local.
- **Chronos joint fit (disclosure):** test-line gene effects come from a Chronos run that also saw training lines.
  The shared parameters are gene-level (guide efficacy, gene means), so per-gene correlation across lines is barely
  affected, while absolute-scale metrics may be slightly favoured. This is stated in the report.

### Inputs (all transforms fitted on training folds only)
- **E** — expression, log2(TPM+1) on DepMap's protein-coding genes. Standardised; genes with training SD < 0.1 dropped.
- **C** — gene-level relative CN as log2((CN + 0.1)/1.1), clipped to [−3, 3], so that 0 = diploid (relative CN 1).
  No per-gene missing indicator: CN NaN is constant per gene, so all-NaN genes are dropped and other NaNs filled
  with 0 (= diploid).
  Derived: genome-instability score on **autosomes only**, computed on raw relative CN (fraction of genes with
  |CN − 1| > 0.3).
- **M** — damaging mutations as 0/1/2 dosage. Kept for genes mutated in ≥5 training lines.
  Derived: log mutation burden.
- **H** — hotspot mutations (binary).
- **F** — fusions. Self-fusions and fusions with non-coding partners are dropped. Gene-partner indicators kept if
  present in ≥5 training lines.
- **Meta** — lineage, primary disease, sex, age.
  - GrowthPattern, ModelType and PrimaryOrMetastasis are **ablation-only**: they proxy the library/NaN group
    (e.g. "Dome" ≡ organoid).
  - Excluded entirely as screen-status proxies: media, FormulationID, CatalogNumber, PublicComments, DataSource,
    AnchorDrug, PassageNumber, EngineeredModel, CulturedResistanceDrug, ModelTreatment.
- **P** — proteomics, from the open Nusinow 2020 CCLE TMT dataset (~375 lines). A subset ablation arm.
- **Other modalities:** signatures, methylation, metabolomics and Sanger MS proteomics need a manual portal
  download. The exact file list goes to the user (`docs/USER_DOWNLOAD_REQUEST.md`); they enter the ablation only if
  supplied.
- **Genotype/state composites (derived):**
  - p53 pathway integrity: TP53 damaging/hotspot status, TP53 CN loss, MDM2/MDM4 amplification,
    PPM1D/ATM/CHEK2 damaging status, plus a p53-target expression score.
  - Rb integrity: RB1 status, plus CDKN2A deep deletion (from CN, not expression).
  - Viral-oncoprotein flag where known (SV40/HPV/EBV; metadata only, for set A tiers).
  - Hallmark expression scores.
- **Convention:** zero-filled genomic features mean "diploid, wild-type". For external cells, engineered genotypes
  (e.g. TP53−/−) are encoded explicitly; everything else uses the convention.

### Gene-level features (for the local and gene-embedding models)
- Own gene: expression, not-expressed indicator (TPM < 1) and hinge, CN, deep-deletion indicator, mutation.
- Paralogs: top-3 paralogs by Ensembl 116 identity — their expression, minimum expression, deep deletion and
  mutation; paralog family size and maximum identity.
- Complex/interaction partners: mean expression of partners from BioPlex 3.0, HuRI and BioGRID physical (excluding
  ORCS); CORUM if obtainable.
- Gene embeddings: built from **training-fold** CRISPR profiles only, plus non-CRISPR prior knowledge.

### Excluded
STRING; OmniPath text-mined layers; BioGRID ORCS as a feature; Project Score; co-essentiality networks not rebuilt
inside folds; DepMap CRISPRInferredCommonEssentials as a feature.

### IDs
Entrez is the master ID, mapped via HGNC. Coverage is reported for every external resource.

## 3. Populations and splits (locked; `configs/splits.csv`, seed 20261026)
- **Cancer:** 1277 screened models with expression, excluding Non-Cancerous models and their patients.
  - **Test** (255 lines): grouped by PatientID, stratified by lineage. Guarded in code (§9) and opened once.
  - **Dev** (1022 lines): 5 grouped, stratified CV folds. All selection happens here, with inner grouped splits for
    hyper-parameters.
- **Set A** — DepMap Non-Cancerous models with RNA-seq. Tiered by transformation mechanism (R1 biology/data):
  - **A1** — hTERT-immortalised, p53/Rb intact. The RPE1-SS sublines form *one genetic background* (5 aneuploid
    sublines), so A1 is reported as a **single-background case study**.
  - **A2** — immortalised with viral or SV40 oncoproteins (HA1E SV40-LT; prostate SV40/HPV lines). Secondary.
  - **A3** — virally transformed lymphoid/myeloid lines (C8166, HG-3, ARH-77, IM-9) and the fibroblast Hs 729
    (origin to be verified). Descriptive.
  - Tiers are assigned from metadata (immortalisation mechanism), never from the lines' own CRISPR outcomes.
    A functional p53 check (TP53 KO > +0.5, MDM2 < −1) is reported only as a descriptive diagnostic.
  - Unit of analysis = line. Each line is ranked within the cancer-test distribution of per-line metrics. No
    per-gene correlation across RPE1 sublines (near-isogenic).
- **Set B** — external published screens (§7).
- **Secondary splits on dev:**
  - Leave-one-lineage-out.
  - Leave-out of the 7 germ-cell/embryonal-carcinoma lines (EC-capability test for hPSC transfer).
    *Erratum (RESEARCH_TREE D20d): only 4 of these are dev lines (NCCIT, NTERA-2 cl.D1, PA-1, SUSA); the
    leave-out uses those 4. Criteria unchanged.*
- **Gene classes (§4):** for *evaluation reporting* everywhere (CV, test, set A, set B), one list computed from all
  1022 dev lines is committed as `configs/gene_classes.csv` before any model is compared.
  - It is identical for every model, so it favours none, and fold differences never reflect class membership.
  - Wherever classes feed into *training* decisions (e.g. grouped λ, M3 gene choice), they are recomputed on each
    training fold.

## 4. Evaluation (pre-registered)
### Evaluable genes
- Cancer dev/test: one fixed evaluable gene set, `configs/eval_genes.csv`, with genes having ≥40 observed lines in
  every dev fold and ≥50 in test.
- Set A and set B: `eval_genes_external`, i.e. genes observed in the screen (and, for set A, in ≥3 of its lines). Test counts come from a logged, mask-only guard call that never returns values. All metrics use this set, so
fold differences never reflect which genes were scored.

### Gene classes
Computed on training data, with prevalences over observed (non-NaN) lines only:
- **all**
- **variable** — SD > 0.25.
- **selective** — 5–95% of lines < −0.5.
- **rare-selective** — ≥3 lines < −0.5 and <5%. Needed because POU5F1-type dependencies fall here.
- **common essential** — < −0.5 in >95% of lines.
- **non-essential** — Hart NEGv1.

### Primary metric P1 (selection)
- **Definition:** mean over genes of per-gene Pearson r across held-out lines, on *variable* genes. Computed
  **within each held-out fold** and then averaged over folds. Pooling out-of-fold predictions gives B0 a spurious
  negative r, so pooling is not used.
- Constant predictions get r = 0.
- **Headline form:** ΔP1 vs B1 (skill over the lineage mean).
- **Also reported:**
  - median r;
  - P1 on *all*, *selective* and *rare-selective* genes;
  - number of genes with r > 0.3;
  - fraction of genes beating B1;
  - within-lineage per-gene r (lineages with ≥20 held-out lines);
  - results by NaN/library group.

### Primary metric P2 (essentiality, label = observed effect < −0.5; sensitivity check at −1.0)
- **P2a, per line:** AUPRC across genes, averaged over lines. Always reported as ΔAP vs B0, per line.
- **P2b, per gene across lines:** for *selective* genes with ≥3 positives and ≥3 negatives.
  - Normalised AP = (AP − π)/(1 − π), plus AUROC.
  - The headline uses genes with held-out prevalence ≤ 0.5.
  - Ties are handled by sklearn's step-wise AP. A unit test asserts that a constant predictor scores AP = π.

### Published-comparable block (not used for selection; labelled "loosely comparable")
- Mean per-gene Spearman on the DeepDEP 1,298-gene and LEAP 1,539-gene sets, where obtainable.
- TrPLet-style count of genes with r ≥ 0.2.
- Shi-style count of genes with r > 0.2 among genes with ≥5 lines of each label.

### Uncertainty
- **CV:** fold-wise values; paired fold-wise differences between models.
- **Test:** 95% bootstrap CIs, stratified by lineage, resampling patients (2000 resamples). Paired patient
  bootstrap for model vs baseline. No tests that treat genes as independent.
- **Ablation claims:** Holm-adjusted paired CIs across the pre-declared contrasts.
- **Power statement for the test set:** the per-gene r standard error at n = 255 is ≈0.06 per gene, but averages
  over ≥2,000 genes. The binding uncertainty is line resampling.

### Selection rule
- Candidates are ordered by complexity: B1 < M1 < M2 < M3/M4 < M5.
- A more complex candidate replaces a simpler one only if its paired fold-wise ΔP1 > δ = 0.005 **and** the mean
  paired difference exceeds 1 SE.
- Exactly one **headline model** is named in `configs/finalists.json` and committed before the test is opened.
  At most 2 additional finalists are allowed, plus baselines.

### Per-gene predictability (diagnostic; added pre-G4, RESEARCH_TREE D22)
This section does not change the headline metric (P1), the selection rule or δ. It is fixed before any test-set
access and before the per-gene distribution is examined for this purpose.

**Source:** dev-CV out-of-fold predictions of the G4 headline model, which is `m1_E_resB0` unless G3 selects
another. Per gene g, compute r_g,k (Pearson r across held-out lines in fold k) for the 5 outer folds.
- mean_r = the mean over folds;
- min_r = the minimum over folds;
- sd_r = the fold SD.

A gene is scored in a fold only where evaluate's min_obs rule makes it evaluable. A gene must be evaluable in all
5 folds to receive a tier; otherwise it is "not evaluable".

**Tiers (frozen to `configs/predictable_genes.csv`):**
- **predictable (high):** mean_r ≥ 0.30 and min_r ≥ 0.20.
- **medium:** not high, mean_r ≥ 0.15 and min_r > 0.
- **low:** all other evaluable genes.

**Confidence score** (shipped with every prediction): mean_r − sd_r, clipped to [−1, 1], plus the tier.

**Null calibration:** the same tiers are computed on the permuted-expression B0-path control
(`ctl_resB0_permuted`). The expected number of false "high" genes is the number of control genes passing. If that
exceeds 5% of the observed high set, the 0.30 threshold is raised to the smallest value in {0.35, 0.40, 0.50} at
which it is ≤ 5%. This is the only allowed adjustment, and it is recorded.

**Characterisation (dev only):** both groups are compared on:
- gene-effect SD and mean, and fraction of lines < −0.5;
- gene class (variable / selective / common essential / non-essential);
- screen reliability: cross-library agreement inside DepMap, with the method declared in config before it is
  computed;
- known biomarker-dependency pairs, from a literature list committed before the comparison;
- own-expression association with the gene's own effect;
- paralog status;
- pathway enrichment (MSigDB Hallmark and Reactome, hypergeometric over the evaluable background, BH).

**Test reporting (G4):**
- P1 and P2b are reported separately for high / medium / low genes;
- also reported: the across-gene Spearman correlation between dev mean_r and test r_g, which checks whether the
  confidence score itself transfers.
- The headline P1 is still the variable-gene mean.

### Mutation-driven dependencies (diagnostic; added pre-G4, RESEARCH_TREE D24)
This section does not change P1, the selection rule or δ. It uses dev CV only.

**Pair set (fixed here):**
- Oncogene addiction, "mutant" = hotspot unless stated otherwise:
  - BRAF (BRAF hotspot), KRAS (KRAS hotspot), NRAS (NRAS hotspot), PIK3CA (PIK3CA hotspot);
  - EGFR (EGFR hotspot or amplification), ERBB2 (ERBB2 amplification), CTNNB1 (CTNNB1 hotspot);
  - CCND1 (CCND1 amplification), MDM2 (TP53 wild-type), MYC (MYC amplification).
- Paralog / collateral loss, "mutant" = damaging mutation or deep deletion of the first gene, dependency on the
  second:
  - ARID1A→ARID1B, SMARCA4→SMARCA2, STAG2→STAG1, ENO1→ENO2, VPS4B→VPS4A;
  - MTAP→PRMT5, MTAP→MAT2A, EP300→CREBBP, CREBBP→EP300.
- The amplification and deep-deletion thresholds are taken from the existing feature definitions
  (src/essml/features.py). Their values are written to the pair file before any metric is computed. A pair
  with < 5 mutant dev lines is reported as "too few" and excluded from pooled statistics.

**Models compared** (existing out-of-fold predictions):
- m1_E_resB0, expression only;
- m1bt_full, expression plus tuned genomics;
- g3_m3n, LightGBM, which sees genomics;
- m5_avg, the headline;
- B1;
- a "biomarker-only" reference: the training-fold mean effect by mutant status, per gene.

**Metrics per pair:**
- per-gene r over all lines;
- r within mutant and within WT lines;
- AUPRC and normalised AP for effect < −0.5, over all lines and within each stratum;
- "biomarker discrimination": AUROC of the prediction for separating mutant from WT among dependent lines.

All metrics are fold-wise and paired.

**"Genomics clearly helps"** holds if both conditions are met:
1. over eligible pairs, the mean paired ΔAUPRC (m1bt_full or g3_m3n minus m1_E_resB0) has a fold-paired 95% CI
   excluding 0;
2. at least 3 pairs individually gain ΔAUPRC ≥ 0.10 in ≥ 4 of 5 folds.

If it holds, one targeted model is specified and pre-registered before it runs:
- M1 plus per-gene own-genomic features (the gene's own hotspot/damaging/CN status and, for paralog pairs, the
  partner's loss status) in the M2 local stage;
- it is evaluated under the unchanged selection rule against the current headline;
- its pair-set gains are reported regardless of the outcome.

### Runtime
Wall-clock fit and predict, thread count and peak RSS, logged by the common launcher.

## 5. Models (each rung must beat the previous one in dev CV; time-boxed, §12)
### Baselines
- **B0** — gene mean.
- **B1** — lineage mean shrunk toward the gene mean, with pseudo-count κ tuned in grouped inner CV. Unseen
  lineage → B0.
- **B2** — kNN in expression space (correlation on top-variance genes; k tuned).

### Published-method re-runs (rung R)
Run on the same folds, for fair prior-work comparison:
- **R1** — Chang & Zhang-style default multi-output ridge, predicting raw y.
- **R2** — Shi-style elastic net (α = 0.5) on variable genes. Time-boxed to 2 h, else a 500-gene subsample.
- **R3** — DepMap-style random forest: top-1000 features chosen inside the fold, on a pre-declared random 500
  variable genes (compute bound).
- **R4** — TrPLet-like top-k correlated-feature selection, nested. Also one non-nested run, to quantify its optimism.

### Our models
- **M1** — multi-output kernel ridge.
  - Linear kernel on standardised features; one `eigh` per fold and NaN pattern; closed-form LOO over λ ∈ [1e2, 1e7].
  - λ is shared within gene groups (variable / other), not free per gene. Per-gene LOO λ gave no gain in the R1
    timing check.
  - NaN patterns whose row sets differ by ≤1% of training lines are merged. Patterns with <50 observed training lines
    get B1 only. Pattern sub-kernels are submatrices of one fold-level K, centred on observed rows. Inner fits are
    cached by (fold, pattern).
  - Target = residual from B1, computed **leave-line-out** in training.
  - **M1b** — a sum of per-modality kernels with block weights.
    - Weights are tuned once per outer fold by coarse coordinate descent on the training LOO criterion (≤25 eigh).
    - Default 1/√p_b weights are also reported.
    - Ablations re-tune λ only.
- **M2** — M1 plus a per-gene local stage: a ridge on out-of-fold M1 predictions plus gene-level features
  (own gene, paralogs, partners). Batched normal equations.
- **M3** (go/no-go) — LightGBM on M2 residuals for variable genes. Inner-fold feature screening; 12 processes ×
  1 thread.
- **M4** (go/no-go) — factorised neural model, y[c,g] = b_g + ⟨u(x_c), v_g⟩ + w·local[c,g].
  - Inputs: cell encoder on 256 training-fold PCs; gene embedding initialised from training-fold Y and prior
    knowledge.
  - Masked MSE, early stopping on an inner split. CPU, 5–15 min per fold.
  - Pre-check first: an MLP on PCs vs ridge on one fold. M4 is dropped if the MLP loses.
- **M5** — stacking. Nested out-of-fold predictions grouped by patient, non-negative weights per gene class.
  Plain averaging is the baseline it must beat.
- **Cut for executability (round 3):** reduced-rank regression and low-rank denoising of Y.
- **Deferred unless time allows:** M1b weight tuning (default 1/√p_b first); R3 random forest; Celligner-style
  expression alignment.
- **Not prioritised:** foundation-model embeddings (beaten by raw expression in 2025 benchmarks), STRING graph networks.

### Transfer model T (pre-registered)
- Uses only features available for external cells: E and engineered genotype (TP53/RB1 status flags).
- **No lineage** is primary. A pre-committed lineage mapping is a secondary variant; mapping hPSC to Germ Cell would
  route B1 through the EC lines and confound H3.
- CN and mutation are set to the "diploid WT" convention.
- **T is trained and tuned on dev like any other model, and reported on the cancer test set beside the full-feature
  model.**

### Frozen models for evaluation
- Cancer test set: the model refit on all 1022 dev lines.
- Set A and set B: the average of the 5 dev fold-models.
- Anchors: ~~the average of the 4 fold-models not trained on that anchor~~ (impossible with 5-fold CV; superseded by
  RESEARCH_TREE D7: set A/B skill is averaged over per-fold models, and each anchor uses its single out-of-fold model). Ensembling is then matched between
  anchors and external screens, so ensembling cannot bias the "≥50% of anchor skill" criterion.
- `configs/finalists.json` freezes the headline model, model T and the H3 no-EC variant at gate G4.

## 6. Ablation (dev CV; λ re-tuned per configuration)
- **Modality add-one** to E: +C, +M, +H, +F, +Meta, +derived, +P (on the P subset).
- **Modality leave-one-out** from the full set.
- **E-free variants** (genomics + Meta only).
- **Prior-knowledge ablation** in M2/M4: none / paralogs / +partners / +gene embeddings.
- **Common line subset:** comparisons between modalities use the same lines (the C∩M∩E intersection; 1206 cancer
  lines overall).
- Each configuration reports P1, P2b and runtime.
- Holm correction over the declared contrasts. Claims of "modality X helps" need an adjusted CI excluding 0.

## 7. Non-transformed generalisation
### 7.1 Assembly
- Least-processed data per screen: counts > LFC > BF. Inventory in `docs/external_screens_inventory.md`.
- **Fitness-readout rule:** only proliferation/viability dropout arms are used; FACS-sorted arms are excluded.
- **Strata:**
  - KO screens are the primary stratum.
  - CRISPRi screens (iPSC neurons and microglia) are descriptive only. They need their own KO↔CRISPRi anchor
    (K562 CRISPRi growth screens only; RPE1 is excluded as an anchor, per §7.2) and serve as a no-DSB control.
  - Post-mitotic neurons: divergence is expected and declared in advance (no proliferation component). Results are
    split into proliferation/cell-cycle genes vs maintenance (housekeeping, mitochondrial, proteostasis) genes.
- **QC gate:** NNMD (or the CEG/NEG control gap) must pass; short screens (e.g. 6-day primary T cells) are gated
  individually.
- **BF/MAGeCK-only screens** get rank metrics and AUPRC only.
- **Screen plan:** before any model is applied, a per-screen table is committed as `configs/screen_plan.csv`.
  Columns: anchors available, readout, matched expression source, data level, and which verdicts are reachable.
  Expected "inconclusive" outcomes are therefore declared in advance.

### 7.2 Screen harmonisation (choices fixed on controls and anchors only, committed before any model predictions on set B)
1. **Uniform re-scoring.** When counts exist, compute guide LFC vs plasmid/T0, then a gene score.
   - Default: median guide LFC.
   - Chronos is refit (fresh, per library) only where anchors show that median LFC fails the §7.2 validation.
2. **Per-screen affine scaling:**
   - scaling set = DepMap common essentials (dev) minus a held-out half of CEGv2 → −1;
   - half of NEGv1 → 0;
   - the other halves are kept for evaluation.
   - Validated by held-out CEG/NEG medians and SSMD, not AUROC (AUROC is invariant to affine maps).
3. **Gene-wise location/scale correction** per lab/library: empirical-Bayes with strong shrinkage, fitted on
   **transformed** anchor lines only.
   - Anchors are DepMap dev lines also screened externally: HAP1, HCT116, DLD-1, A-375, K-562, HT-29, AML lines.
   - RPE1 and all non-transformed lines are **excluded** (circular).
   - HepG2/MCF-7 are test lines, so they are used for expression validation only.
   - Applied only if leave-one-anchor-out validation improves anchor agreement measured on **deviations from the
     gene mean**.
   - With fewer than 8 anchors per library, the correction is location-only.
   - Anchor skill (§7.6) is always computed from profiles corrected **without that anchor**.
   - The p53/DNA-damage-response gene set (TP53, CDKN1A, TP53BP1, USP28, CHEK2, ATM, MDM2, MDM4, PPM1D and the
     p53-target score genes) is excluded from fitting and application. All anchors are cancer lines without
     competent p53, so a correction could otherwise absorb p53-competent cutting toxicity as "batch". H1/H1′ are
     reported before and after correction.
4. **Duration** is not treated as a global scale. Residual time-dependent bias is checked on anchors. Short screens
   are gated by QC.
5. **Validation outputs:**
   - held-out CEG/NEG medians and SSMD;
   - NEG median ≈ 0;
   - anchor same-line agreement with DepMap (deviation correlation), before and after each step;
   - same-cell-type cross-study agreement (Hart 2015 / Brown 2019 / Olivieri 2020 RPE1; DepMap RPE1-SS vs external
     RPE1). This agreement is the **replicate ceiling**.
   - Simulated-batch control: a known affine plus gene-wise distortion applied to DepMap anchor profiles must be
     undone.

### 7.3 Expression harmonisation
- Recompute TPM on **all genes**, then subset to DepMap's protein-coding set, then log2(TPM+1). This replicates
  DepMap's convention (per-line TPM sums ≈0.94M).
- Candidates, chosen on matched ENCODE/GEO profiles of DepMap *dev* lines (K-562, HCT116, A-549 …) only:
  - (i) none;
  - (ii) quantile normalisation to DepMap;
  - (iii) quantile + gene-wise location correction fitted on matched lines;
  - (iv) rank-expression features;
  - (v) Celligner-style alignment (contrastive PCA + mutual nearest neighbours), if it can be implemented without
    tumour data.
  - Candidate (iii) is the ComBat-type recipe of Shi et al.'s GTEXDEPMAP (quantile + ComBat), the nearest prior
    art for transfer to normal tissue.
- Criterion: same-line nearest-neighbour recall and profile correlation.
- Profiles must be state-matched to the screened cells (same cell type and state). QuantSeq/3′ data are flagged as
  non-TPM-compatible.
- An **applicability-domain distance** (to the nearest dev lines in harmonised expression space) is reported for
  each external profile.

### 7.4 Core-essential transfer
Labels: harmonised effect < −0.5, plus a **scale-free rank label** (the 10% most depleted genes per screen) that does
not depend on CEG-based scaling.
Per screen, using observed harmonised labels (effect < −0.5), with the evaluation halves of CEGv2/NEGv1 also
reported:
- **Primary evidence (Q3a):** the absolute AUPRC of model T, per screen and per label, relative to the replicate
  ceiling (the AUPRC of an independent same-cell-type screen used as the predictor).
- **Secondary:** ΔAUPRC vs B0. It is structurally small for core essentials, which B0 already captures.
- Also reported: AUPRC/AUROC of T and the full model vs B0, B1 and B2 (kNN in harmonised expression space), plus the
  observed replicate ceiling.
- **Lost/gained essentials:** cancer common essentials that are not essential in the screen, and vice versa. Did
  the model predict each shift (direction and rank)?
- RPE1 results are reported with and without CEGv2 genes, because CEGv2 partly derives from RPE1/HCT116 screens.

### 7.5 Context-specific dependencies
- **Primary statistic (Q3b):** Δrecall vs B1 (and vs B0) of screen-specific dependencies, at a fixed 5% FPR on the
  **evaluation half** of NEGv1.
  - Screen-specific = observed < −0.5 **and** dependency fraction ≤ 0.5 in cancer dev lines (headline).
  - Secondary definition: dependency fraction ≤ 5% (rare).
  - Clopper–Pearson CIs for recall. The ceiling is a same-cell-type replicate screen. B2 is also a comparator.
- **Deviation correlation:** across *selective + rare-selective* genes, corr(pred − gene mean, obs − gene mean).
  Note B0 has deviation ≡ 0.
- **Nulls:**
  - **profile-swap** — predict each screen from other screens' expression profiles. Matched must beat mismatched.
  - **nearest-lineage comparator** — the B1 prediction for the closest cancer lineage.
- **Pre-specified hypotheses** (thresholds in §7.6):
  - **H1, TP53 pathway, learnable part:** MDM2/PPM1D dependency is predicted stronger in TP53-WT backgrounds.
    Scored by rank/slope. USP7 is excluded: it is broadly essential (<−0.5 in 72% of lines).
  - **H1′, suppressor tail as a diagnostic:** positive effects for TP53/CDKN1A/TP53BP1/USP28/CHEK2/ATM. This is a
    screen diagnostic, not a model test. The observed effects are 4–8× beyond the cancer range.
  - **H2, isogenic TP53 counterfactual:** the RPE1 TP53-WT vs TP53-null delta (Brown 2019; Bowden 2020 if obtained)
    against the model's delta with the TP53 feature flipped. Batch-cancelling.
  - **H3, pluripotency:** POU5F1/LIN28A-type dependencies in hPSC.
    - Tested with the model trained with and without the 7 germ-cell/EC lines (erratum D20d: the 4 EC dev lines).
    - If they are predicted only when the EC lines are present, the claim is "narrow EC transfer".
  - **H4, haematopoietic positives:** GATA3/MYB/RUNX1/IKZF1 dependency in the primary CD8+ T-cell screen
    (Shifrut 2018).
    - Each gene must be expressed in the matched profile.
    - TCR-stimulation coupling is declared.
    - The unexpressed-gene false-positive rate is the control.
    - CD34+ (Ting 2018) is used only if its deposit resolves; otherwise H4 for CD34+ is inconclusive by declaration.
    - SPI1 belongs to the descriptive microglia CRISPRi stratum.
- **Negative controls / declared expected failures:**
  - oncogene addictions (KRAS/BRAF/CTNNB1) predicted ≈0 (a negative control, near-tautological);
  - TCR-signalling genes in T-cell screens;
  - the suppressor tail;
  - neuron-specific biology.
- **CN null expectation:** RPE1-SS sublines differ by chromosome *gains*, while collateral lethality needs *loss*,
  and Chronos already corrects CN. Collateral/paralog genes (SMARCA2, ARID1B, STAG1, ME2/ME3, ENO1) are therefore
  expected near 0 in both prediction and observation. This is a null check.

### 7.6 Decision rule (thresholds committed in `configs/transfer_criteria.json` before set-B predictions)
- Only Q3a and Q3b get primary verdicts. Each is one of **transfers / fails / inconclusive**.
- Verdicts are given per **biological background**: the RPE1 datasets are one background, hPSC lines another,
  primary T cells another. The overall verdict is "transfers" if a majority of evaluable backgrounds transfer.
- H1–H4 are secondary. Criteria and thresholds are committed in `configs/transfer_criteria.json` at G5, before any
  set-B prediction. A hypothesis whose criterion is not committed by then is reported descriptively only.
- **Thresholds also fixed in that file at G5, from anchors and controls only:**
  - QC: NNMD < −1.0, or control-gap SSMD for non-Avana libraries;
  - shrinkage prior;
  - anchor-agreement improvement margin;
  - MDE: 2 × the bootstrap SE of the replicate-ceiling skill.
- **Transfers:** the lower CI bound of skill over B0 (Q3a) or B1 (Q3b) is > 0, **and** skill is ≥ 50% of the
  same-pipeline anchor skill.
- Comparison: skill over B0/B1 in non-transformed screens vs the same skill in transformed anchors processed through
  the same lab/library pipeline.
- The minimum detectable effect comes from the replicate ceiling and the anchor noise.
- Uncertainty: patient bootstrap where several lines exist; otherwise a gene block bootstrap (blocks = gene
  families/complexes).
- "Fails" requires all four:
  - skill ≤ B0/B1 bound;
  - the upper CI bound lies below the MDE;
  - the anchors in the same pipeline show positive skill;
  - harmonisation validation passed.
- **Expression-swap negative control:** must show no skill. If swapped profiles still show skill, the verdict is
  capped at "inconclusive".
- Otherwise: inconclusive.

## 8. Verification and controls
- **Unit tests (toy data, known answers):**
  - metrics: perfect → r = 1, AP = 1; constant → r = 0, AP = π; random → AP ≈ π; B0 per-fold r = 0;
  - split integrity;
  - loader default entries;
  - pattern-wise masked ridge equals per-gene ridge on observed rows;
  - leave-line-out lineage mean.
- **Leakage assertions:** fit-time functions receive only training ModelIDs. An assertion checks that no
  test/set-A/set-B IDs appear in any fit call. Scaling, selection, gene classes and lineage means are computed
  inside folds.
- **Positive control:** planted synthetic targets (y = f(own expression) + noise, with known r). The pipeline
  recovers the planted r.
- **Negative controls:**
  - features permuted across training lines → P1 ≈ 0;
  - random features;
  - near-duplicate check: the max feature correlation between test and dev lines is reported; dev-only reference
    sets.
- **Harmonisation controls:** §7.2 step 5.
- **Probe bugs become tests:** validation-fold λ selection, and imputation that uses validation lines.

## 9. Process safeguards (agent-failure reviews, rounds 1–2)
- **Frozen evaluation code:**
  - `src/essml/{metrics,evaluate,genesets,guard}.py` are hashed into `configs/frozen_hashes.json` at freeze, and
    git tag `freeze-eval` marks the commit.
  - Transfer and harmonisation code are frozen the same way under tag `freeze-transfer`, together with
    `configs/harmonisation_frozen.json`.
  - Tests fail if a hash changes. Any post-freeze change is logged in RESEARCH_TREE with reason, and baselines are
    re-run.
- **Test guard:**
  - Targets for test, set A or set B are released only to a run whose config hash is listed in a
    `configs/finalists.json` committed earlier, on a clean tree, with `ESSML_UNBLIND=1`.
  - Exactly one logged access per finalist (test) and per finalist × screen (set A/B, staged per screen); every
    access goes to `results/test_access_log.csv`.
  - Logged exceptions:
    - mask-only counts (`observed_counts`);
    - control-gene reads (CEGv2/NEGv1/scaling-set genes only) of external screens for harmonisation and QC before
      the transfer freeze.
  - A static test fails if any module other than `io.py`/`guard.py` reads the CRISPR file directly.
  - Set-A and set-B prediction code refuses to run until `configs/harmonisation_frozen.json` and
    `configs/transfer_criteria.json` are committed.
- **Gate evidence:** each gate G0–G7 writes `results/gates/G<k>.json` (checks run, pass/fail, commit). The next gate's
  scripts refuse to start without it.
- **Run registry:** a single launcher (`scripts/run.py`) writes `results/registry.csv`.
  - A row is written at run start and updated at the end.
  - Columns: run id, config hash, git commit, dirty-tree flag, splits.csv hash, roles used, n training IDs,
    metrics, runtime, status.
  - Failed runs are included. Leaky diagnostic runs (e.g. non-nested R4) are flagged `leaky=1`. Every number in REPORT.md cites a run id,
  and `scripts/check_report_numbers.py` verifies them.
- **Independent reproduction:** headline numbers are re-computed by a separate agent writing its own metric code.
  Review prompts never state the expected answer. A red-team agent is assigned to each key finding.
- **Literature:** `[unverified]` claims are never used for design decisions without checking.
- **External-data intake:** each new dataset gets a DATA.md row, sha256 and licence before use.
- **Post-freeze bug protocol:** fix, re-run all affected runs including baselines, and record the change.
- **Exclusion log:** every dropped line, gene or screen and the reason goes to `results/exclusions.csv`.
- **Scope ledger:** every planned arm (models, ablations, hypotheses, screens) with status run / skipped (reason)
  is kept in `results/scope_ledger.csv` and reproduced in REPORT.md.
- **Implementation check:** at G0 and G4, an agent verifies that each accepted review item in
  `reviews/round*_response.md` is implemented, and writes `reviews/implementation_check.md`.
- **Stopping:** stop a rung when two consecutive ideas fail to beat the current best by δ in dev CV, or when its time
  box expires. Stop the project when the idea list is exhausted or the time budget nears its end.

## 10. Compute plan (12 cores, 16 GB; memory < 4 GB per stage)
| Stage | Estimated cost |
|---|---|
| M1 | ~10 s per fold per configuration (eigh + LOO grid) |
| M2 | < 1 min |
| M3 | ~15–20 min per CV |
| M4 | 5–15 min per fold |
| R2/R3 | ≤ 1 h each (bounded gene sets) |
| Chronos fits | ≤ a few hours each |

- Threads: env vars are set before numpy is imported; `threadpoolctl` is used. BLAS stages run 12 threads;
  process-parallel stages run 1 thread each. `__main__` guards are used throughout.

## 11. Risks
- Thin iPSC-derived and primary categories.
- Expression-pipeline mismatch.
- CRISPRi vs KO incomparability.
- TP53-WT cutting toxicity.
- Small, tiered set A.
- Short primary screens.
- Missing DepMap modalities (user download).

## 12. Work plan and gates
| Gate | Work | Pass condition |
|---|---|---|
| G0 | Verification suite + frozen metrics + launcher/registry + test guard | All tests pass |
| G1 | B0/B1/B2 and R1 on dev CV | Sanity: B0 r = 0; B1 > 0 |
| G2 | M1, M1b, M2 + ablation | Choose features |
| G3 | M3/M4 go/no-go (time box 1 day each); M5 | Choose headline model |
| G4 | Commit finalists.json, gene_classes.csv | Unblind test once; independent reproduction |
| G5 | In parallel from now: external acquisition → harmonisation validation | Commit criteria/transfer_criteria.json |
| G5a | Set A (batch-free) is evaluated before set B | — |
| G6 | Apply frozen models to set A and set B | Verdicts per §7.6 |
| G7 | REPORT.md | Number check passes; red-team review |
