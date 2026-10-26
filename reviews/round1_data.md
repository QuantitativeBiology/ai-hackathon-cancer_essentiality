# Round 1 review: data correctness (DESIGN.md v0)

## 1. Verdict
Sound in outline, but needs revision before modelling. The leakage and circularity rules are mostly right: Meta excludes the obvious provenance columns, and prior knowledge is restricted to non-CRISPR, non-text-mined sources. Five problems need fixing. (a) The NaN handling is not random and mean-fill is wrong. (b) Three retained Meta fields still encode the screening batch. (c) The transfer model is not tied to features that external cells actually have. (d) Several transforms and filters are mis-specified: CN, fusions and the TPM recipe. (e) Set A is not a clean "non-transformed" population. All numbers below were computed from the local cache, `configs/splits.csv` and Model.csv.

## 2. Proposed changes

**1. [MUST] §2 Targets: use masked fitting, not training-mean fill.**
- The NaNs are library-structured, not random. Only 49 distinct gene-level NaN patterns exist. Clustering lines on their NaN signature gives 4 groups: 827, 221, 118 and 195 lines. They are almost certainly screen-library or batch groups.
  - Group 2 (118 lines) is the Sanger-source group: 49 of its lines have Model.csv SourceType Wellcome Sanger Institute, and its ModelCondition DataSource is SANGER or BROAD|SANGER. 1086 genes are NaN in more than 50% of its lines.
  - Group 1 contains all 108 organoids.
- Mean-filling shrinks those lines toward the mean at 1086 genes and makes the NaN pattern a learnable batch signal.
- Change: group genes by NaN pattern (49 groups) and solve the ridge or kernel ridge exactly on the observed rows per group. For M3 and M4, use a masked loss. Drop the "fill with training gene mean" sentence.
- Evaluation: require at least 50 observed held-out lines per gene for P1. Test-set observed counts per gene have median 255, and 519 genes have fewer than 50.
- Evaluation: report P1 split by NaN group. Only 9 test lines are in group 2, so state that group-2 genes are evaluated on the other groups.
- Fully observed variable genes show a small group effect, about 1% of variance explained against 0.17% for permuted labels. This is minor but non-zero.

**2. [MUST] §2 Meta: remove screen-provenance proxies from the default feature set.**
- GrowthPattern "Dome" equals organoid, 108 of 108 lines, all in NaN group 1. "Spheroid" is 40 of 40 in group 1. "Unknown" is 37 of 41 in group 2 (Sanger).
- PrimaryOrMetastasis "Recurrent" is 13 of 13 in group 1.
- Expression alone predicts organoid status with 99.4% cross-validated accuracy (91.6% baseline), so ModelType adds nothing legitimate beyond expression.
- Change: the default Meta set is lineage, OncotreePrimaryDisease, Sex and Age only. Age is 17% missing, so impute with a missing flag. GrowthPattern, ModelType and PrimaryOrMetastasis become ablation-only items, labelled "provenance-risk".
- Report the ablation with and without them, and a model-predicts-NaN-group check.
- Sex is not a leakage problem. CN of Y genes is about 1e-6 in females and about 0.3 in males, so sex is recoverable from CN anyway.

**3. [MUST] §7 and §5: define a transfer-capable model.**
- External cells only have expression, plus known genotype such as TP53 status. They have no CN, mutations, fusions or DepMap Meta. Lineage is "unseen" for iPSC-derived cells, primary haematopoietic cells and RPE1.
- Change: pre-register a transfer model trained on E-only, optionally with a small known-genotype set (TP53 and a few drivers, set by literature). Prior-knowledge features are allowed.
- Report the full-feature and E-only models side by side in dev CV, so the cost of the restriction is visible.
- Set A also lacks CN for 7 of the 15 lines, so only 8 have full omics.

**4. [MUST] §2 C: fix the CN transform, the missing-value rule and the line count.**
- Values are relative to ploidy (median 1.02) and heavy-tailed: 0.09% of values are below 0.1 (homozygous-deletion-like), the maximum is 260, and the minimum is 1e-6.
- log2(CN+0.01) ranges from −6.64 to 8.02. Its 0.1% and 99.9% quantiles are −2.08 and 1.92, so a few values dominate a linear model.
- Change: use log2(CN+0.1), or clip log2(CN+0.01) to [−2, 3], then standardise on the training fold.
- NaN is gene-level and constant across lines: 278 genes are NaN in all 1763 lines, and every line has exactly 1.48% NaN. A per-gene "missing" indicator is therefore meaningless. Drop those 278 genes, and add a single line-level "CN absent" flag only for the 71 lines that have expression but no CN.
- The "1214 intersection" counts 8 non-cancerous lines. For the cancer population it is 1206.
- Derived genome-instability score: use autosomes only, since X and Y reflect sex.

**5. [MUST] §2 F: replace the ">20% of models" fusion filter.**
- It removes only 3 of 63,832 unique fusions. The top fusions are self-fusions or read-throughs (TFAM--TFAM in 520 lines, ACIN1--ACIN1, KRT7--KRT8), and 55,896 fusions are singletons.
- Change: drop fusions where Gene1 equals Gene2. Drop fusions with a non-coding or pseudogene partner (names like RP11-*). Keep only fusions seen in at least 5 training lines (1416 fusions in total). Optionally require FFPM above a threshold, with the threshold tuned in dev CV.
- Better: encode fusions at gene level ("gene g is in a recurrent fusion") instead of as 60k sparse partner pairs.

**6. [MUST] §7.3: match DepMap's TPM definition before normalising.**
- DepMap per-line TPM sums to a median of 0.94M (range 0.63M to 0.97M), not 1M. So TPM was computed over all genes and then subset to protein-coding.
- Change: for external RNA-seq, compute TPM over the full annotation, then subset. Use the same transcript model (GENCODE and RSEM/STAR, as DepMap does if documented; confirm).
- Matching the pipeline is safer than aggressive quantile normalisation. Quantile normalisation forces every sample onto DepMap's distribution, which can erase real global differences in ESCs and neurons. Make it a validated, optional step.
- `data/external_expr/tian2019_wtc11_neuron_quantseq` is QuantSeq (3' counting). It is not TPM-compatible, so treat it as a limitation or use a rank-based or CPM variant.
- Add the same-line validation (K-562, HCT116 and so on) as a go/no-go gate before any external prediction.

**7. [MUST] §3 and §7: Set A is "Non-Cancerous", not necessarily non-transformed.**
- It contains C8166 (HTLV-1) and HG-3 (EBV), which are virally transformed. It contains 5 RPE1-SS lines from one patient, 5 SV40- or hTERT-immortalised prostate lines, HA1E, and an adenoma organoid.
- Hs 729 is annotated Fibroblast, but I recall ATCC listing it as rhabdomyosarcoma-derived; verify before including it as a normal.
- The effective number of independent backgrounds is about 5, not 15.
- Change: report Set A per background (RPE1, prostate, HA1E, B-lymphoblast, adenoma). Resample by PatientID or background in the bootstrap. Exclude virally transformed lines from the headline "non-transformed" claim, and keep ARH-77 and IM-9 out because they have no expression.

**8. [SHOULD] §2 M and H: specify feature construction.**
- Damaging values are 0 / 1 / 2, with counts 46.9M / 122k / 12.6k. At value 2 the median CN is 0.85 and only 7.8% have CN below 0.5, so 2 is not a CN-loss artefact. Keep ordinal 0/1/2 dosage rather than collapsing to binary.
- 694 genes are never mutated and 920 CRISPR genes are absent from the matrix. Treat absence as 0, and keep only genes mutated in at least 5 training lines (or a driver list).
- Mutation burden spans 1 to 1727 per line (99th percentile 671). Use log burden.
- Hotspot covers only 165 genes (TP53 hotspot is in 1067 lines).

**9. [SHOULD] §2 prior knowledge: specify ID mapping and circularity.**
- Use Entrez as the master key, matching `io.py`, and map Ensembl IDs (Compara, HuRI) to Entrez through HGNC.
- Report coverage. 284 CRISPR genes are absent from CN and 69 from expression. 561 CN columns are Ensembl-only and unparseable; only 4 of them map via HGNC, so that loss is negligible.
- Fusion partners carry ENSG IDs and need the same mapping.
- BioPlex (293T and HCT116) is expression-biased, so include it only as an ablation level, as planned.
- CEGv2 and NEGv1 were derived partly from hTERT-RPE1, HCT116 and DLD-1 screens (Hart 2015/2017). Treat them as an independent control only for lines not in that derivation, and report RPE1 core-essential AUPRC with and without CEGv2.
- CORUM is CC BY-NC (acceptable for research, no redistribution) and was unobtainable by script. Use hu.MAP as the fallback, and state which one was used.

**10. [SHOULD] §7.2 anchors: use only dev lines.**
- HAP1, HCT116, DLD-1, A-375, K-562, HT-29, A549, THP-1 and the AML/B-ALL lines are all in the dev set. HepG2 and MCF-7 are test lines.
- Change: fit the gene-wise correction only on dev anchors. Evaluate anchor predictions out-of-fold, since the final model has trained on those lines' targets.
- Use HepG2 and MCF-7 only for expression validation, never for CRISPR correction.
- There are 7 usable anchors, so ComBat-style gene-wise correction is under-powered for 18k genes. Use heavy shrinkage or a gene-class-level correction.

**11. [COULD] §3 and §6: handle organoids explicitly.**
- Organoids (108) differ from cell lines. Within the lineages that have organoids, 479 genes differ by more than 0.3 in mean effect, with a median of 99 observed organoid lines per gene against 277 cell lines. This is confounded with library and medium.
- Change: keep organoids in the pool for model fitting (expression separates them cleanly), but report P1 for cell lines only and for organoids only. Treat the organoid-ablation difference as a sensitivity, not a finding.

**12. [COULD] §2 CRISPR metadata: add the screen-source map.**
- A cluster map built from the NaN pattern is an approximation of the screen source.
- If the user can supply `ScreenSequenceMap.csv` and `CRISPRScreenMap.csv` (already in the "manual download" list), replace the NaN-cluster proxy with the actual library, so batch-stratified reporting is exact.
