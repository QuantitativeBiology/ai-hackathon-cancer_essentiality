# LITERATURE — prior work, datasets and implications for the design

Survey date: 2026-10-26. Four independent agent reviews hold the detail and full reference lists; this file is the
synthesis. References in prior_methods.md and nontransformed_screens.md carry `[opened]`, `[search]` or `[unverified]` tags. harmonisation_labels.md and data_resources.md are largely untagged, so treat their claims as unverified until checked. Treat
`[unverified]` items as leads, not facts.

| Note | Scope |
|---|---|
| [docs/lit/prior_methods.md](docs/lit/prior_methods.md) | Predicting dependency from omics; reported accuracy; evaluation pitfalls |
| [docs/lit/nontransformed_screens.md](docs/lit/nontransformed_screens.md) | CRISPR screens in hTERT lines, ESC/iPSC (and derived), primary blood cells |
| [docs/lit/harmonisation_labels.md](docs/lit/harmonisation_labels.md) | Screen scoring, cross-screen batch correction, essentiality labels |
| [docs/lit/data_resources.md](docs/lit/data_resources.md) | DepMap files, open omics mirrors, prior-knowledge resources, licences, circularity |

## 1. Predicting dependency in unscreened lines

- **Expression carries most of the signal.** DepMap's random-forest study (Dempster et al. 2020, bioRxiv) and elastic
  nets (Shi et al. 2024) found that adding mutations or copy number to expression adds little and sometimes hurts.
  The best per-gene models are sparse, driven by the gene's own expression, its paralogs' expression and lineage.
- **Linear models are hard to beat on held-out lines.** Multi-output ridge reached mean per-gene ρ 0.276 against
  0.137 for DeepDEP on the same data (Chang & Zhang 2023 re-analysis). LASSO-type models top the LEAP benchmark
  (mean per-gene Spearman 0.321, DepMap 23Q4, 1,539 genes). A 2025 benchmark found raw expression beat
  single-cell foundation-model embeddings (Geneformer, CancerFoundation).
- **Reported numbers vary with gene set and split.** Means over all ~17k genes are ≈0.2 (TrPLet, 23Q2). Means over
  selected variable genes are 0.3–0.4. Pooled correlation over all (line, gene) pairs is ≈0.9, but that comes from
  differences between gene means and says nothing about predicting across lines. TrPLet's numbers are optimistic
  (best of 60 models per gene, chosen inside the cross-validation it reports).
- **Gaps, corrected after review round 1:**
  - Shi et al. 2024 already used a confounder null that includes lineage.
  - The EFFECT benchmark (Szalai et al. 2023) already scores with partial correlations controlling for line and
    gene means.
  - What remains uncommon is an explicit shrunk lineage-mean *baseline* reported beside the model, CV grouped by
    PatientID, and AUPRC for *selective* dependencies in held-out lines.
  - Similarity transfer (kNN, matrix completion) is barely benchmarked on DepMap.
  - Our contribution is the combination, not any one of these.
- **Known biological structure:** oncogene addiction (own expression), paralog buffering (De Kegel & Ryan 2019;
  Dede 2020; Ito 2021), CN-driven cutting toxicity (largely removed by Chronos, residual remains), complex
  co-dependency, lineage transcription factors, TP53 status (TP53-WT cells respond to DNA cuts).

## 2. Labels and metrics
- DepMap's probability-of-dependency uses unexpressed genes as its null, so it is partly circular with expression
  features. That file is also not in the local data. Gene effect thresholds (−0.5; −1 as sensitivity) avoid this.
- Pooling AUPRC over all (line, gene) pairs mostly rewards ranking common essentials. Report it per line (across
  genes) and per gene (across lines, for selective genes) separately.
- Hart CEGv2/NEGv1 (hart-lab/bagel, MIT) are standard core-essential and non-essential reference sets. DepMap uses
  overlapping lists to scale Chronos, so any set used to scale an external screen must be kept apart from the set
  used to evaluate it.

## 3. Non-transformed screens available
- **hTERT:** RPE1 dominates (Hart 2015 TKOv1; Brown 2019 TP53 WT/null; Olivieri 2020 RPE1 TP53−/− TKOv3, with
  untreated arms). No genome-wide KO fitness screen was found in BJ-hTERT, HME1, MCF10A or IMR-90.
- **ESC/iPSC:** Mair 2019 (H1, TKOv3), Yilmaz 2018 (haploid hESC), Ihry 2019 (hPSC).
- **iPSC-derived:** Tian 2021 neurons and Dräger 2022 microglia, but both are CRISPRi and must be analysed apart from
  knockout. This is the thinnest category.
- **Primary haematopoietic:** Shifrut 2018 primary CD8+ T (Brunello; short and partly FACS-based), Ting 2018 CD34+
  HSPC. These are short screens with a compressed dynamic range.
- **Inside DepMap:** 17 screened "Non-Cancerous" models (15 with RNA-seq), including 5 RPE1 aneuploid sublines,
  HA1E and immortalised prostate epithelium. These share DepMap's library and scoring, so they test transfer with
  no screen batch effect.
- **BioGRID ORCS** (MIT) aggregates curated screen scores with cell type, library and duration metadata. It may be
  used as a test-set source, never as a model input.

## 4. Harmonisation
- When counts are available, re-score with Chronos fitted fresh per library: DepMap's pretrained guide efficacies
  cover only Avana/KY. Then scale so common essentials have median −1 and non-essentials median 0. Then apply
  gene-wise location/scale correction (ComBat-like) fitted on anchor lines screened both in DepMap and by the
  external lab. Pacini et al. 2021 report that correction improved with anchor number; 8 was the smallest set they tested (not a recommended minimum), and 30–40 nearly matched full correction.
- Affine scaling and gene-wise ComBat preserve variation within a gene across lines. Quantile or rank transforms
  preserve only ranks, so screens available only as BF/MAGeCK scores support rank metrics and AUPRC only.
- **Anchors present in 26Q3 CRISPR data:** HAP1, HCT116, DLD-1 (Hart TKO), HT-29 and AML lines (Tzelepis KY),
  A-375 (Brunello), K-562 (Wang), RPE1-SS sublines. HEK293T, HeLa, MCF10A and BJ-hTERT are *not* in
  CRISPRGeneEffect.
- **Validating harmonisation:** CEG/NEG separation on held-out control genes; NEG genes centred at 0; same-cell-line
  recall across datasets (Dempster 2019: nAUC 0.53→0.91 after correction).

## 5. Data access and circularity
- 26Q3 is published only through the portal, which has a Cloudflare bot check (no Figshare mirror after 24Q4). Extra
  DepMap modalities need a manual download by the user: OmicsSignatures, Sanger MS proteomics, methylation (RRBS,
  CCLE legacy), metabolomics, CN segments, CRISPRScreenMap/ScreenSequenceMap, CRISPRGeneDependency, Achilles controls.
- **Open prior knowledge:** Ensembl Compara paralogs, HGNC, BioGRID physical interactions, HuRI, BioPlex 3.0,
  Reactome, MSigDB Hallmarks. CORUM is CC BY-NC and its download endpoint returned 503.
- **Excluded as circular or text-mined:** STRING, OmniPath's text-mined layers, BioGRID ORCS as a feature, Project
  Score matrices, co-essentiality networks, DepMap's CRISPRInferredCommonEssentials, synthetic-lethal labels learned
  from CRISPR. Any of these is acceptable only if rebuilt from training lines inside each fold.

## 6. Implications carried into DESIGN.md
1. Baselines: gene mean, shrunk lineage mean, and an expression-kNN baseline.
2. Primary model family: cheap linear models (multi-output kernel ridge on residuals), then gene-specific features
   (own gene, paralogs, complexes), then nonlinear learners only where they beat the linear model in CV.
3. Grouped-by-patient evaluation, with per-gene metrics split by gene class (all / variable / selective / common essential).
4. For non-transformed transfer: evaluate within-DepMap non-cancerous lines first (no batch), then harmonised
   external screens, with anchor-based validation and a cross-screen replicate ceiling.
