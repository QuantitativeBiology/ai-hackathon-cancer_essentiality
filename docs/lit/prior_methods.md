# Prior work: predicting gene dependency from cell-line molecular features

Scope: predicting the CRISPR gene-effect profile of a never-screened line from its omics, plus the evaluation conventions behind the published numbers. Compiled 2026-10-26. Tags: **[opened]** means the paper or full text was read. **[search]** means the bibliographic details and claims come only from search-result snippets or abstracts. **[unverified]** means not confirmed. Reference numbers [n] point to the list at the end.

---

## 1. DepMap's own predictive modelling

- **Tsherniak et al. 2017 (Cell) [1] [search].** RNAi (DEMETER) screens in 501 lines. 769 genes were "differentially required" (6σ outlier criterion). Nonlinear regression on 66,646 features gave predictive models for 426 of them (55%), and in 82% of models the top biomarker was an expression feature. This was the first evidence that expression beats genomics.
- **Dempster et al. 2020 (bioRxiv) [2] [opened].** Setup: one random forest per perturbation (100 trees, depth ≤ 8, ≥ 5 lines per leaf), trained on the top 1,000 features by |Pearson r| inside the training fold. Validation: 10-fold CV over cell lines, with the out-of-fold predictions concatenated and scored by Pearson r for each perturbation; "predictable" means r ≥ 0.4. Data: Achilles, Project Score, DEMETER2 RNAi, GDSC17 and PRISM, with 241–649 lines per dataset.
  - **Features:** expression (genes plus ssGSEA) won almost everywhere. 1,641 perturbations were predictable from expression and 784 from DNA. Adding DNA to the best expression model usually made it slightly *worse*.
  - **Accuracy:** median r for strongly selective perturbations ranged 0.27–0.49 across datasets. Elastic net was only slightly worse than RF (RF won 355 of 500).
  - **Model structure:** the best models are usually sparse, driven by the gene's **own expression** or a **paralog's expression**. Mutation was the best feature only for RAS/BRAF-type oncogenes.
- **DepMap portal "Predictability" [3,4] [opened, forum].** Per-gene RFs, top-1,000 feature filter (f_regression in practice per forum), 5-fold CV, Pearson r. Model variants ("Core Omics", "DNA-based", "Related") and whether "Related" uses STRING-like networks are *[unverified]*. Staff: predictions not stored; low predictability ≠ artefact; correlated features split importance.
- **Genome Biol 2023 [5] [search; authors unverified].** Per-gene models over 15,221 genes (CRISPR and RNAi), "accurate" = r > 0.5; removing confounder-driven models changed the counts.

## 2. Deep learning and "foundation" approaches (2021–2026)

| Method | Architecture | Split | Key point |
|---|---|---|---|
| **DeepDEP**, Chiu et al. 2021, *Sci Adv* [6] [opened] | Autoencoders for mutation, expression, methylation and CN, pre-trained on 8,238 TCGA tumours, plus a gene "fingerprint" (membership in 3,115 MSigDB CGP gene sets) | Held-out lines: 90/10 (28 test lines), also 10-fold | Pooled ρ = 0.87, but **mean per-gene ρ = 0.18** (0.16 under 10-fold) on 1,298 genes, 278 lines, 2018Q2 CERES |
| **Chang & Zhang 2023**, bioRxiv [7] [opened] | Default multi-output **ridge** on concatenated omics, no fingerprint | 10-fold over lines | Pooled ρ 0.88. **Per-gene ρ 0.276 vs 0.137 for DeepDEP**. Ridge on expression alone also beats DeepDEP. |
| **DeepVul**, Jararweh et al. 2024, bioRxiv [9] [opened] | Transformer encoder over expression, multi-task (essentiality → drugs) | **Split not specified** | "46% of genes with Spearman > 0.3", on the top ~1,000 variable genes. The evaluation is too poorly specified to compare. |
| **LEAP**, Bodinier et al. 2026, *Sci Rep* [10] [opened] | Ensemble of 5 masked-autoencoder representations × regressors (LASSO, MLP, LGBM) | 10 × 80/20 **held-out lines**, plus a tissue-held-out test | 23Q4: 1,019 lines, 1,539 genes, **mean per-gene Spearman 0.321** (PS-LASSO). On 22Q4: 0.297 vs 0.256 for OmicsRPZ. Ensembling over representations adds about +0.012. MAE pre-training needs a GPU. |
| **Zinchenko et al. 2025** (Bayer), bioRxiv [13] [opened] | Bi-modal NN benchmark of cell-line encoders | Held-out lines; selective perturbations only | **Raw expression beat Geneformer, CancerFoundation and scVI.** The foundation models did not beat KEGG-pathway scores. Gene encoders matter only for held-out *genes*, and STRING was best there (forbidden for us). |

Also seen but not evaluated (numbers/splits not retrieved): Itzhacky & Sharan 2021 expression NNs [8]; a transformer+PPI-GNN [11] and GATDep graph attention on 24Q4 [12], both network-based (STRING-type sources would be excluded).

Gene-side embeddings (LLM- or PLM-derived, e.g. GenePT and ESM) only help when the task is predicting **unseen genes**. In our setting every gene is seen in training and only the *line* is new, so a per-gene model already holds all the gene-specific information. I found no paper showing that foundation-model cell embeddings beat raw bulk expression for predicting dependency in held-out lines [13]. EGBench [14] [search] probes Geneformer and ESM2 embeddings for gene-level (not line-level) essentiality.

**Takeaway:** at matched held-out-line evaluation, linear models on expression match or beat published deep models [7,10,15]. The deep models' headline numbers are pooled correlations.

## 3. Similarity and transfer methods for unscreened samples

- **Per-gene linear models transposed to tumours.**
  - *Shi et al. 2024 (Nat Cancer), TCGADEPMAP [16] [opened].* Elastic net (α = 0.5) on 1Q21 CERES data, 897 lines, 10-fold CV. Of 7,260 genes with ≥ 5 dependent lines, **1,966 passed r > 0.2** (FDR < 0.001) with expression only and 2,045 with multi-omics. 76% of genes were within Δr = 0.05 between the two feature sets, so multi-omics barely helps. Transfer to TCGA used quantile normalisation plus contrastive PCA, removing 4 stromal components. *They used STRING only downstream, for interpretation.*
  - *TrPLet, Sadagopan et al. 2024/25 (bioRxiv) [15] [opened].* RNA only, 23Q2 Chronos. For each gene, the best of 60 models (EN, ridge, lasso, linear and RBF SVR × 12 top-k correlation feature-set sizes) was chosen by **the same 5-fold CV that is reported, with no nested CV, so the result is optimistically biased**. **6,283 of 16,845 genes reached r ≥ 0.2. Mean per-gene r = 0.20 and pooled r = 0.92.** On DeepDEP's 1,298 genes, mean r was 0.37 vs 0.18. Adding mutations made no difference. Validated prospectively in screened tRCC lines.
- **DeepDR / DeepDEP tumour transfer** (Chiu et al. 2019, *BMC Med Genomics* [17] [search], drug response; and [6]). Encoders pre-trained on TCGA are meant to reduce the cell-line-to-tumour shift. This is unvalidated for dependencies beyond the internal checks.
- **Celligner** (Warren et al. 2021, *Nat Commun* [18] [search]). Unsupervised alignment of 1,249 DepMap lines with 12,236 tumours that removes the tumour-vs-line axes (purity, stroma). It is useful as an **expression alignment step before nearest-neighbour transfer** and for projecting non-cancer samples. It also shows that several hundred lines take on a mesenchymal, undifferentiated state that has distinct dependencies.
- **Nearest neighbours and matrix completion with side information.** I found **no dedicated DepMap paper** that benchmarks kNN in expression space or inductive matrix completion for new lines. LEAP uses PS-KNN (k = 5) only as a weak baseline [10]. The general literature notes that plain low-rank completion fills an entirely missing column with a constant, so side information (expression) is required, which is equivalent to reduced-rank regression. *This is a gap we can fill cheaply.*

## 4. Known biological predictive structure

- **Own-gene expression and oncogene addiction.** The gene's own expression was the top feature in a large share of models [2]. Shi et al. found the target itself among the selected features in 26–31% of models [16]. Pacini et al. 2024 (*Cancer Cell*) [19] [search], working with 930 lines, report that most marker–dependency links are **gain-of-function addictions, not synthetic lethality**.
- **Lack of expression implies a null effect.** DepMap's probability-of-dependency model uses scores of *unexpressed genes* as part of its null [20,21]. Expression = 0 is therefore a near-deterministic predictor of "not dependent", and this partly shapes the label itself (see §5).
- **Paralog buffering.**
  - Paralogs are less often essential than singletons, and 13–17% of variable paralog essentiality is attributable to buffering by the other paralog (De Kegel & Ryan 2019) [22] [search].
  - About half of constitutively expressed genes are never hits, enriched for paralogs, and 24 synthetic-lethal pairs were found by enCas12a (Dede et al. 2020) [23] [search].
  - DUSP4/DUSP6 is a digenic dependency in MAPK-driven lines (Ito et al. 2021) [24] [search].
  - A paralog-SL classifier generalises to unseen lines with AUC around 0.9 (Kebabci et al. 2026) [25] [opened], **but its network features come from STRING v11.5. We may use only the paralog expression, CN and sequence-identity parts** (sequence identity alone gives AUC 0.54–0.57).
  - Paralog annotations themselves (Ensembl Compara) are not text-mined and are allowed.
- **Copy-number CRISPR artefact.** Cas9 cutting in amplified regions reduces fitness independently of the gene (Aguirre et al. 2016) [26] [search]. CERES [27] and Chronos [20] correct for this, but residual CN-driven signal can still make CN look "predictive" for genes inside amplicons. **This should not transfer to largely diploid non-transformed cells.**
- **Collateral lethality.** MTAP is co-deleted with CDKN2A, and its loss sensitises cells to PRMT5 inhibition (Mavrakis et al. 2016) [28] [search]. Kryukov et al. 2016 [unverified] reported the same link. The original passenger-deletion concept (ENO1/ENO2, Muller et al. 2012 [unverified]) is classic. It requires **gene-level deletion calls and expression of the partner**, and it is captured by own-gene and partner expression features.
- **Complex membership / co-essentiality.** Complex members share dependency profiles; paralogs in complexes are more often SL [22]. Useful mainly as a **multi-output / low-rank prior**.
- **Lineage TFs and lineage.** Lineage-restricted dependencies are well captured by expression [2]; lineage was in every Dempster model [2] and in Shi's confounder null [16]. It is the first thing a model learns, hence the lineage-mean baseline.

## 5. Evaluation pitfalls (evidence-based)

1. **Pooled vs per-gene correlation.** Pooled ρ (0.87–0.92) mostly reflects between-gene differences, i.e. the per-gene mean. Per-gene means are about 0.14–0.37 on the same predictions [6,7,15]. LEAP explicitly rejects pooled metrics [10]. *We should report the full per-gene distribution and the mean/median over gene subsets: all genes, strongly selective genes, and high-variance genes.*
2. **Gene subset definition.** Published means depend heavily on the subset: 1,298 variable or Cancer Gene Census genes [6], 1,539 genes [10], selective perturbations [2,13], or all 16,845 genes [15]. **Select "variable/selective" genes using training folds only**, because selecting on the full matrix leaks test-line variance.
3. **Model selection inside the reported CV** (best-of-60 per gene [15], or top-k feature filtering outside the folds) inflates r. Use **nested CV**.
4. **Lineage confounding.** Lineage was a covariate or confounder control in [2,16]. *None of the papers read report a per-gene lineage-mean baseline.* Our design should therefore report (a) Δ over the lineage mean, (b) within-lineage per-gene r, and (c) leave-lineage-out, as in LEAP's tissue hold-out [10].
5. **Multiple models per patient.** DepMap's `Model.csv` `PatientID` links ModelIDs derived from the same individual [29] [opened forum]. **Use GroupKFold on PatientID.** I found no published dependency-prediction paper that does this; it is standard leakage hygiene.
6. **Split type.** Random (line, gene) entry splits leak the line and the gene. DeepVul's split is unspecified [9]. Benchmarks that hold out genes [13] answer a different question. Only held-out-line numbers are comparable to our task.
7. **Essentiality threshold and AUPRC.**
   - The probability of dependency (> 0.5) is a mixture-model posterior whose null uses unexpressed genes [20,21]. Expression features are therefore partly circular with this label. Report results with **gene effect < −0.5 (or −1)** as well.
   - Pooled AUPRC is dominated by common-essential genes, so the per-gene mean already scores high. Report base rates, AUPRC on **selective dependencies only**, and per-gene AUPRC for genes with ≥ k dependent test lines (Shi et al. required ≥ 5 dependent and ≥ 5 non-dependent lines [16]).
8. **Screen and library batch.** Achilles and Project Score libraries disagree for some genes [6]. Sanger-vs-Broad agreement on shared lines is ρ 0.66 (per [6]), which sets a practical ceiling. Treat screen source as a covariate and evaluate within source.
9. **Release and algorithm drift.** CERES (2018–2021) vs Chronos (2021+) and changing line counts make numbers across papers only loosely comparable.

## 6. Reported numbers to beat (held-out cell lines unless noted)

| Method | Data | Split | Metric | Value |
|---|---|---|---|---|
| DeepDEP [6] | 2018Q2 CERES, 278 lines, 1,298 genes | 10-fold over lines | mean per-gene ρ (pooled) | 0.16–0.18 (0.87) |
| Linear/EN/SVM baselines as run by [6] | same | same | mean per-gene ρ | 0.03–0.11 |
| Multi-output ridge [7] | same | 10-fold over lines | mean per-gene ρ (pooled) | **0.276** (0.88) |
| TrPLet [15] | 23Q2 Chronos, ~1,100 lines | 5-fold, best-of-60 (optimistic) | mean per-gene r, 16,845 genes / on DeepDEP's 1,298 genes | 0.20 / **0.37**; 6,283 genes r ≥ 0.2; pooled 0.92 |
| Shi et al. EN [16] | 1Q21 CERES, 897 lines | 10-fold | # genes r > 0.2 (of 7,260) | 1,966 (RNA), 2,045 (multi-omics) |
| Dempster RF [2] | ~2019 Achilles/Score | 10-fold over lines | median per-gene r, strongly selective perturbations | 0.27–0.49; 1,641 predictable perturbations at r ≥ 0.4 (all datasets) |
| LEAP PS-LASSO [10] | 23Q4, 1,019 lines, 1,539 genes | 10 × 80/20 lines | mean per-gene **Spearman** | **0.321** |
| LEAP PS-LASSO / OmicsRPZ [10] | 22Q4, 893 lines, 1,223 genes | same | mean per-gene Spearman | 0.297 / 0.256 |
| PS-KNN (k = 5) [10] | 23Q4 | same | mean per-gene Spearman | < 0.321; exact value only in a figure |
| Tsherniak [1] | RNAi, 501 lines | not retrieved | fraction of differential dependencies predictable | 426/769 (55%) |
| Paralog-SL RF [25] | 22Q4 features, Ito et al. pairs | unseen lines | AUROC (uses STRING) | ~0.90 |

No paper reports **AUPRC for selective dependencies in held-out lines**, a **lineage-mean baseline**, or **PatientID-grouped CV**. All three are open ground for us.

## 7. Recommendations for our design (12 cores, 16 GB RAM, no GPU)

**Data scale:** about 1,100 lines × 19k expression genes and 18k Chronos genes, roughly 80 MB each in float32. All linear models fit comfortably.

**Baselines.**
- B0: per-gene mean.
- B1: per-gene lineage mean, shrunk toward the gene mean for small lineages (empirical Bayes), falling back to B0 for unseen lineages.
- B2: kNN in expression space (correlation distance on top-variance genes, optionally Celligner-aligned), with k tuned in nested CV. This is the transfer baseline that §3 found missing from the literature.

**Strong, cheap models (shortlist, in priority order).**
1. **Multi-output kernel/ridge regression on expression**, predicting residuals from B1 (or B0).
   - Fit: one SVD of the n × p training matrix, then λ per gene by closed-form leave-one-out/GCV. All 18k genes take seconds to minutes.
   - Evidence: this is what beat DeepDEP [7]; LASSO-type linear models top LEAP [10].
   - Variants: expression PCs (100–500), plus CN and mutation blocks with separate penalties (group ridge), plus a lineage one-hot.
2. **Per-gene sparse models with prior-knowledge features**: elastic net or ridge on [own expression, own CN, own damaging mutation, paralog expression and CN (Ensembl Compara), hotspot mutations of ~50 drivers, lineage], plus top-k training-fold-correlated expression features.
   - Evidence: this mirrors the RF findings that a few features dominate [2] and TrPLet's top-k recipe [15], but without best-of-60 selection bias.
   - Compute: about 18k fits × folds, parallelised over 12 cores, should be minutes to tens of minutes (estimate).
3. **Reduced-rank / low-rank regression**: factorise the training-line residual dependency matrix (rank 50–300), predict the factor scores from expression (ridge), then reconstruct. This shares strength across co-essential modules and complexes, and is equivalent to matrix completion with side information.
4. **Gradient-boosted trees (LightGBM, CPU) or RF on the selective-gene subset only** (~1–3k genes, curated features of fewer than 1,000), to catch nonlinearity such as expression-zero thresholds and mutation × lineage interactions. RF matches EN closely in [2], so expect small gains at a higher runtime cost. Budget-check it first.
5. **Stacked ensemble** of models 1–4 with nested CV (LEAP-style diversity gives about +4% [10]).

**Deprioritise:** DeepDEP-style autoencoders (beaten by ridge [7]); transformers and single-cell foundation-model embeddings (no gain over raw bulk expression [13], GPU-bound); LLM/protein gene embeddings (irrelevant when genes are seen); PPI GNNs (usually STRING); methylation/proteomics as primary features (expression dominates [2,16]; fewer lines; ablation only).

**Features that matter, by expected value:**
1. Expression (log TPM, standardised).
2. Lineage.
3. Own-gene and paralog CN and deletion calls (collateral lethality, CN artefact).
4. Driver hotspot mutations (RAS/BRAF/PIK3CA etc.).
5. Fusions, as a small block for lineage-defining drivers (e.g. TFE3 [15]).
6. Signatures, methylation and proteomics, for the ablation only.

**For the non-transformed transfer:** prefer expression-only models, because CN-artefact and mutation features will be out of distribution in diploid cells. Align expression first (Celligner-style or quantile/ComBat). Report core-essential recovery separately from selective dependencies.

---

## References

1. Tsherniak A, et al. Defining a Cancer Dependency Map. *Cell* 2017;170(3):564–576. doi:10.1016/j.cell.2017.06.010 [search]
2. Dempster JM, Krill-Burger JM, McFarland JM, et al. Gene expression has more power for predicting in vitro cancer cell vulnerabilities than genomics. *bioRxiv* 2020. doi:10.1101/2020.02.21.959627 [opened]
3. DepMap forum, "About PREDICTABILITY". https://forum.depmap.org/t/about-predictability/3215 [opened]
4. DepMap forum, "Confounders in the Predictability Tab". https://forum.depmap.org/t/confounders-in-the-predictability-tab/2026 [opened]
5. Partial gene suppression improves identification of cancer vulnerabilities when CRISPR-Cas9 knockout is pan-lethal. *Genome Biol* 2023. doi:10.1186/s13059-023-03020-w [search; authors unverified]
6. Chiu YC, et al. Predicting and characterizing a cancer dependency map of tumors with deep learning. *Sci Adv* 2021;7(34). PMC8378822 [opened]
7. Chang D, Zhang X. Ridge regression baseline model outperforms deep learning method for cancer genetic dependency prediction. *bioRxiv* 2023. doi:10.1101/2023.11.29.569083 [opened]
8. Itzhacky N, Sharan R. Prediction of cancer dependencies from expression data using deep learning. *Mol Omics* 2021;17(1):66–71. doi:10.1039/d0mo00042f [search]
9. Jararweh A, et al. DeepVul: a multi-task transformer model for joint prediction of gene essentiality and drug response. *bioRxiv* 2024. doi:10.1101/2024.10.17.618944 [opened]
10. Bodinier B, et al. Predicting gene essentiality and drug response from preclinical perturbation screens with layered ensemble of autoencoders and predictors (LEAP). *Sci Rep* 2026;16:26261. doi:10.1038/s41598-026-56381-0; arXiv:2502.15646 [opened]
11. Modeling cancer dependency with deep graph models. *bioRxiv* 2024. doi:10.1101/2024.02.26.582022 [search; content not retrieved]
12. Context-aware gene dependency modeling via graph attention networks for precision oncology (GATDep). *J Transl Med* 2025. doi:10.1186/s12967-025-07501-3 [search]
13. Zinchenko V, Schlicker A, Kurilov R, Pouplin A, Villalba S, Horlacher M. Benchmarking chemical, genetic, and cell line encodings for cancer perturbation response prediction. *bioRxiv* 2025. doi:10.64898/2025.12.15.694331 [opened]
14. EGBench. https://github.com/tubic/EGBench [search]
15. Sadagopan A, et al. Cancer target discovery enabled by transcriptome-based virtual CRISPR screening (TrPLet). *bioRxiv* 2024 (v3 2025). doi:10.1101/2024.10.24.620074 [opened]
16. Shi X, Gekas C, Verduzco D, et al. Building a translational cancer dependency map for The Cancer Genome Atlas. *Nat Cancer* 2024;5(8):1176–1194. doi:10.1038/s43018-024-00789-y [opened]
17. Chiu YC, et al. Predicting drug response of tumors from integrated genomic profiles by deep neural networks. *BMC Med Genomics* 2019;12(Suppl 1). doi:10.1186/s12920-018-0460-9 [search]
18. Warren A, et al. Global computational alignment of tumor and cell line transcriptional profiles. *Nat Commun* 2021;12:22. PMC7782593 [search]
19. Pacini C, et al. A comprehensive clinically informed map of dependencies in cancer cells and framework for target prioritization. *Cancer Cell* 2024. doi:10.1016/j.ccell.2023.12.016 [search]
20. Dempster JM, et al. Chronos: a cell population dynamics model of CRISPR experiments that improves inference of gene fitness effects. *Genome Biol* 2021;22:343. doi:10.1186/s13059-021-02540-7 [search]
21. Dempster JM, et al. Extracting biological insights from the Project Achilles genome-scale CRISPR screens in cancer cell lines. *bioRxiv* 2019 (cited by DepMap forum for probability of dependency) [search]
22. De Kegel B, Ryan CJ. Paralog buffering contributes to the variable essentiality of genes in cancer cell lines. *PLoS Genet* 2019. doi:10.1371/journal.pgen.1008466 [search]
23. Dede M, McLaughlin M, Kim E, Hart T. Multiplex enCas12a screens detect functional buffering among paralogs otherwise masked in monogenic Cas9 knockout screens. *Genome Biol* 2020;21:262. doi:10.1186/s13059-020-02173-2 [search]
24. Ito T, et al. Paralog knockout profiling identifies DUSP4 and DUSP6 as a digenic dependence in MAPK pathway-driven cancers. *Nat Genet* 2021;53:1664–1672. doi:10.1038/s41588-021-00967-z [search]
25. Kebabci N, Ajmal H, Adams DJ, Ryan CJ. A predicted cancer dependency map for paralog pairs. *bioRxiv* 2026. doi:10.64898/2026.01.19.700065 [opened]
26. Aguirre AJ, et al. Genomic copy number dictates a gene-independent cell response to CRISPR/Cas9 targeting. *Cancer Discov* 2016;6(8):914–929 [search]
27. Meyers RM, et al. Computational correction of copy number effect improves specificity of CRISPR–Cas9 essentiality screens in cancer cells. *Nat Genet* 2017;49:1779–1784 [search]
28. Mavrakis KJ, et al. *Science* 2016 — MTAP loss sensitises cells to PRMT5 inhibition via MTA accumulation [search; exact title unverified]. Kryukov GV, et al. *Science* 2016 [unverified]. Muller FL, et al. *Nature* 2012 (ENO1 passenger deletion) [unverified]
29. DepMap forum, "DepMap Metadata" / "DepMap Data Structure" (PatientID definition). https://forum.depmap.org/t/depmap-metadata/3694 [search]
