# DESIGN — predicting CRISPR gene essentiality for unscreened human cells

Version: **v1** (after review round 1; responses in `reviews/round1_response.md`), 2026-10-26.
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
- **C** — gene-level relative CN as log2(CN + 0.1), clipped to [−3, 3]. No per-gene missing indicator: CN NaN is
  constant per gene, so all-NaN genes are dropped and other NaNs filled with 0 (= diploid).
  Derived: genome-instability score on **autosomes only** (fraction of genes with |log2CN| > 0.3).
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
  - p53 pathway integrity: TP53 damaging or hotspot status, plus a p53-target expression score.
  - Rb integrity: RB1 status, plus CDKN2A expression.
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
    sublines). Primary claim set.
  - **A2** — immortalised with viral or SV40 oncoproteins (HA1E SV40-LT; prostate SV40/HPV lines). Secondary.
  - **A3** — virally transformed lymphoid/myeloid lines (C8166, HG-3, ARH-77, IM-9) and the fibroblast Hs 729
    (origin to be verified). Descriptive.
  - Unit of analysis = line. Each line is ranked within the cancer-test distribution of per-line metrics. No
    per-gene correlation across RPE1 sublines (near-isogenic).
- **Set B** — external published screens (§7).
- **Secondary splits on dev:**
  - Leave-one-lineage-out.
  - Leave-out of the 7 germ-cell/embryonal-carcinoma lines (EC-capability test for hPSC transfer).
- **Gene classes (§4)** for test, set A and set B: computed once from all 1022 dev lines and committed as
  `configs/gene_classes.csv` before the test is opened. In CV they are recomputed on each training fold.

## 4. Evaluation (pre-registered)
### Gene classes
Computed on training data:
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
- **R2** — Shi-style elastic net (α = 0.5) on variable genes.
- **R3** — DepMap-style random forest: top-1000 features chosen inside the fold, on a pre-declared random 500
  variable genes (compute bound).
- **R4** — TrPLet-like top-k correlated-feature selection, nested. Also one non-nested run, to quantify its optimism.

### Our models
- **M1** — multi-output kernel ridge.
  - Linear kernel on standardised features; one `eigh` per fold and NaN pattern; closed-form LOO over λ ∈ [1e2, 1e7].
  - λ is shared within gene groups (variable / other), not free per gene. Per-gene LOO λ gave no gain in the R1
    timing check.
  - Target = residual from B1, computed **leave-line-out** in training. Lineage is unpenalised.
  - **M1b** — a sum of per-modality kernels with block weights (tuned), to integrate C/M/F/Meta.
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
- **Diagnostic only:** reduced-rank regression; low-rank denoising of Y.
- **Not prioritised:** foundation-model embeddings (beaten by raw expression in 2025 benchmarks), STRING graph networks.

### Transfer model T (pre-registered)
- Uses only features available for external cells: E, engineered genotype (TP53/RB1 status flags) and a fixed
  lineage mapping (or no lineage).
- CN and mutation are set to the "diploid WT" convention.
- **T is trained and tuned on dev like any other model, and reported on the cancer test set beside the full-feature
  model.**

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
    (K562 Horlbeck 2016; RPE1 Replogle 2022 if obtainable) and serve as a no-DSB control.
  - Post-mitotic neurons: divergence is expected and declared in advance (no proliferation component).
- **QC gate:** NNMD (or the CEG/NEG control gap) must pass; short screens (e.g. 6-day primary T cells) are gated
  individually.
- **BF/MAGeCK-only screens** get rank metrics and AUPRC only.

### 7.2 Screen harmonisation (choices fixed on controls and anchors only, committed before any model predictions on set B)
1. **Uniform re-scoring.** When counts exist, compute guide LFC vs plasmid/T0, then gene score by Chronos (fit fresh
   per library) or median guide LFC. The choice is made on anchors.
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
  - (iv) rank-expression features.
- Criterion: same-line nearest-neighbour recall and profile correlation.
- Profiles must be state-matched to the screened cells (same cell type and state). QuantSeq/3′ data are flagged as
  non-TPM-compatible.
- An **applicability-domain distance** (to the nearest dev lines in harmonised expression space) is reported for
  each external profile.

### 7.4 Core-essential transfer
Per screen, using observed harmonised labels (effect < −0.5), with the evaluation halves of CEGv2/NEGv1 also
reported:
- AUPRC and AUROC of model T and the full model **vs B0 and B1**, plus the observed replicate ceiling.
- **Lost/gained essentials:** cancer common essentials that are not essential in the screen, and vice versa. Did
  the model predict each shift (direction and rank)?
- RPE1 results are reported with and without CEGv2 genes, because CEGv2 partly derives from RPE1/HCT116 screens.

### 7.5 Context-specific dependencies
- **Headline statistic:** recall of screen-specific dependencies at a fixed NEGv1 false-positive rate (5%).
  - Screen-specific = observed < −0.5 and not common-essential in the cancer dev set.
  - Measured against the replicate ceiling and against B0/B1.
- **Deviation correlation:** across *selective + rare-selective* genes, corr(pred − gene mean, obs − gene mean).
  Note B0 has deviation ≡ 0.
- **Nulls:**
  - **profile-swap** — predict each screen from other screens' expression profiles. Matched must beat mismatched.
  - **nearest-lineage comparator** — the B1 prediction for the closest cancer lineage.
- **Pre-specified hypotheses** (thresholds in §7.6):
  - **H1, TP53 pathway, learnable part:** MDM2/PPM1D/USP7 dependency is predicted stronger in TP53-WT backgrounds.
    Scored by rank/slope; the effect size is 4–8× beyond the cancer range.
  - **H1′, suppressor tail as a diagnostic:** positive effects for TP53/CDKN1A/TP53BP1/USP28/CHEK2/ATM. This is a
    screen diagnostic, not a model test.
  - **H2, isogenic TP53 counterfactual:** the RPE1 TP53-WT vs TP53-null delta (Brown 2019; Bowden 2020 if obtained)
    against the model's delta with the TP53 feature flipped. Batch-cancelling.
  - **H3, pluripotency:** POU5F1/LIN28A-type dependencies in hPSC.
    - Tested with the model trained with and without the 7 germ-cell/EC lines.
    - If they are predicted only when the EC lines are present, the claim is "narrow EC transfer".
  - **H4, haematopoietic positives with DepMap support:** e.g. MYB, LMO2/GATA2-type and SPI1-type dependencies in
    CD34+/blood screens. The unexpressed-gene false-positive rate is the control.
- **Negative controls / declared expected failures:**
  - oncogene addictions (KRAS/BRAF/CTNNB1) predicted ≈0 (a negative control, near-tautological);
  - TCR-signalling genes in T-cell screens;
  - the suppressor tail;
  - neuron-specific biology.
- **CN-collateral test:** where aneuploid RPE1-SS lines differ in CN, check whether predicted dosage effects match.

### 7.6 Decision rule (thresholds committed in `configs/transfer_criteria.json` before set-B predictions)
- For each question, the verdict is one of **transfers / fails / inconclusive**.
- Comparison: skill over B0/B1 in non-transformed screens vs the same skill in transformed anchors processed through
  the same lab/library pipeline.
- The minimum detectable effect comes from the replicate ceiling and the anchor noise.
- "Fails" requires all three:
  - skill ≤ B0/B1 bound;
  - the anchors in the same pipeline show positive skill;
  - harmonisation validation passed.
- **Expression-swap negative control:** must show no skill.

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

## 9. Process safeguards (R1 agent-failure review)
- **Frozen evaluation code:** `src/essml/metrics.py` and `src/essml/evaluate.py` are hashed into
  `configs/frozen_hashes.json` at freeze. Tests fail if a hash changes. Any post-freeze change is logged in
  RESEARCH_TREE with reason, and baselines are re-run.
- **Test guard:**
  - the loader refuses test/set-B targets unless `ESSML_UNBLIND=1` *and* `configs/finalists.json` exists in git;
  - every access is appended to `results/test_access_log.csv`;
  - set A and set B predictions use frozen models only.
- **Run registry:** a single launcher (`scripts/run.py`) writes `results/registry.csv`. Columns: run id, config
  hash, git commit, metrics, runtime, status. Failed runs are included. Every number in REPORT.md cites a run id,
  and `scripts/check_report_numbers.py` verifies them.
- **Independent reproduction:** headline numbers are re-computed by a separate agent writing its own metric code.
  Review prompts never state the expected answer. A red-team agent is assigned to each key finding.
- **Literature:** `[unverified]` claims are never used for design decisions without checking.
- **External-data intake:** each new dataset gets a DATA.md row, sha256 and licence before use.
- **Post-freeze bug protocol:** fix, re-run all affected runs including baselines, and record the change.
- **Exclusion log:** every dropped line, gene or screen and the reason goes to `results/exclusions.csv`.
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
| G6 | Apply frozen models to set A and set B | Verdicts per §7.6 |
| G7 | REPORT.md | Number check passes; red-team review |
