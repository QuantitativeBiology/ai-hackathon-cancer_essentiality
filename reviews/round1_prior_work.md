# Review round 1: use of prior work (DESIGN.md v0)

Reviewer: prior-work agent (gene-dependency prediction). Date: 2026-10-26.
Inputs read: PROMPT.md, DESIGN.md v0, LITERATURE.md, docs/lit/prior_methods.md. New sources were checked on the web and
are cited only where I saw the page or abstract.

## 1. Verdict

**Sound overall; needs revision before locking.** The model ladder follows the evidence well: linear models first,
expression-led, own-gene and paralog features, nested CV, and gene classes defined on training folds only. Five
problems remain:

- (a) No published recipe is re-run on our splits. Without that, any "beats prior work" claim compares different
  data.
- (b) The survey's own top recommendation, reduced-rank regression (matrix completion with side information), was
  dropped from §5.
- (c) Two of the three novelty gaps are overstated.
- (d) The nearest prior art for non-cancer transfer, GTEXDEPMAP, is not used, and Celligner is not offered.
- (e) Several cheap, documented feature effects are missing from M2.

## 2. Proposed changes

**1. §5, add rung "R: published-recipe reproductions" on the locked dev folds. (must)**
- R-RF, the DepMap/Dempster recipe: per gene, take the top 1,000 features by |r| inside the training fold, then a
  random forest with 100 trees, depth ≤ 8, min leaf 5. Run it on *variable* genes only (runtime-capped).
- R-EN, the Shi/TCGADEPMAP recipe: elastic net (α = 0.5) on expression only, and again on multi-omics.
- R-TrPLet: top-k correlation features with linear/ridge/SVR. Run it nested, and once non-nested to measure the
  optimism of best-of-60 selection.
- R-ridge, Chang & Zhang: default multi-output ridge on raw targets, not residuals.

Rationale: only these give an apples-to-apples comparison. The published numbers mix releases (CERES vs Chronos),
line sets, gene sets and Spearman vs Pearson (prior_methods §5.9, §6). Dempster 2020 (doi:10.1101/2020.02.21.959627);
Shi 2024 (doi:10.1038/s43018-024-00789-y); TrPLet (doi:10.1101/2024.10.24.620074); Chang & Zhang
(doi:10.1101/2023.11.29.569083).

**2. §5, add M1b: reduced-rank regression.** Factorise the training-fold residual Y at rank r ∈ {25…300}, predict the
factor scores from features with ridge, then reconstruct. Report it next to B2 (kNN). (should)
- Rationale: LITERATURE §1 claims similarity and matrix-completion transfer as "a gap we can fill cheaply". As
  written, DESIGN fills it only with kNN.
- prior_methods §7 ranked this third in its shortlist.
- It costs seconds.

**3. §4 Metrics, define how P1 aggregates across folds. (must)** State whether per-gene r is computed on the
out-of-fold predictions concatenated over the 5 folds, or averaged per fold.
- Recommendation: make the concatenated version P1, because it is the Dempster/Shi/TrPLet convention. Keep the
  per-fold mean for the SE in the selection rule.
- The two can differ materially for genes with few dependent lines.

**4. §4, add a "published-comparable" reporting block, kept outside model selection. (must)**
- Mean per-gene **Spearman** on the DeepDEP 1,298-gene set and the LEAP 1,539-gene set, using the supplement lists
  if obtainable.
- Mean r over all genes, and the number of genes with r ≥ 0.2 (TrPLet's form).
- The number of genes with r > 0.2 among genes with ≥5 dependent and ≥5 non-dependent lines (Shi's form).
- Label every external number "loosely comparable". The matched comparison is change 1. LEAP:
  doi:10.1038/s41598-026-56381-0.

**5. LITERATURE §1 and §6, and DESIGN §4: correct the novelty claims and add a bias-controlled metric. (must)**
- **A lineage null already exists.** Shi 2024 judged expression models against DepMap "confounder" models built
  from sex, Cas9 activity, age, lineage, primary/metastasis and others. Their r > 0.2 pass threshold corresponds to
  P < 0.03 in that null distribution (PMC11358024, seen).
- **Bias-controlled scoring already exists.** The EFFECT benchmark (Szalai et al. 2023,
  doi:10.1101/2023.10.02.560281, seen) scores with partial correlations that control for per-line and
  per-perturbation means. It also reports that, on a cell-line-exclusive CRISPR split, "few genes score above 0.3".
- **Per-gene thresholds already exist.** Shi used per-gene ≥5/≥5 criteria.
- What remains genuinely new: a *shrunk lineage mean as an explicit baseline, with Δ reported per gene*;
  PatientID-grouped CV; and per-gene AUPRC on selective genes.
- Restore *within-lineage per-gene r*, or the partial r given lineage, as a secondary metric. prior_methods §5.4
  recommended it and DESIGN dropped it.

**6. §3/§4, stratify evaluation by screen library or source, and add a cancer-side ceiling. (should)**
- Map each line to Avana, KY or other via CRISPRScreenMap. Report P1 and P2b per stratum. This is an evaluation
  stratum only, never a feature: Cas9 activity and similar screen-derived quantities are CRISPR data from the test
  line.
- For lines screened in two libraries, report agreement between those screens as the empirical ceiling. This
  mirrors §7.2.5.
- Rationale: Broad and Sanger agree at only about ρ 0.66 on shared lines (prior_methods §5.8). Pacini et al. 2021
  integrated the two datasets.

**7. §7.3 and §7.4, add prior-art expression alignment and pre-register the transfer model. (should)**
- Compare three expression alignments:
  - plain quantile normalisation;
  - quantile normalisation + ComBat, which is Shi's GTEXDEPMAP recipe for non-malignant tissue (PMC11358024);
  - Celligner (cPCA + MNN; Warren 2021, PMC7782593), run with and without removing the cPCs. Cancer-vs-normal axes
    may carry real biology.
- Select among them on §7.3 validation and on DepMap set A, which has no screen batch.
- State that GTEXDEPMAP was validated only qualitatively (toxicity-target essentiality, KRAS/BRAF in normal
  tissues), with no normal-cell screens. That is the advance §7 offers; cite it.
- Pre-register an **expression-only** variant of the finalist as a transfer model, as prior_methods §7 recommends.
  CN and mutation features are out of distribution in near-diploid cells.
- Also add B2 (kNN) as a transfer baseline. Note that B1 falls back to B0 in this setting.

**8. §5 M2, add features with documented effects. (must)**
- (a) An own-gene *not-expressed* indicator (TPM < 1) and a hinge term. A linear model cannot represent "no
  expression → effect ≈ 0". DepMap's dependency null is built from unexpressed genes, and own expression is the top
  RF feature (prior_methods §4).
- (b) Own-gene and paralog **deep-deletion** indicators, plus paralog damaging-mutation status. This covers
  collateral lethality (MTAP; Mavrakis 2016) and paralog loss. Relative log CN alone blurs homozygous loss.
- (c) For M4: paralog family size and maximum identity as gene covariates. These are the non-STRING parts of
  Kebabci 2026 (doi:10.64898/2026.01.19.700065) and De Kegel 2021 (doi:10.1101/2020.12.16.423022).

**9. §2, promote proteomics from "if the user supplies" to a planned ablation arm. (should)**
- Nusinow 2020 (about 375 lines) is openly downloadable now, per data_resources.md.
- Gonçalves 2022 (949 lines; doi:10.1016/j.ccell.2022.06.010, abstract seen) reports that the proteome predicts as
  well as the transcriptome. It also reports thousands of protein biomarkers of vulnerabilities that are not
  significant at transcript level.
- This is the modality a referee will ask about. Run it on the matched line subset (§6 rule).

**10. §2 Prior knowledge: replace "CORUM if obtainable" with hu.MAP 3.0 and EBI Complex Portal. (could)**
- hu.MAP 3.0 comes from co-fractionation MS; Complex Portal is manually curated. Neither is text-mined.
- CORUM is CC BY-NC and its download endpoint returned 503 (LITERATURE §5).

**11. §5, note the limits of gene-side encoders in M4. (could)** Zinchenko 2025 (doi:10.64898/2025.12.15.694331) found
that gene encoders help only for held-out genes. M4 should therefore justify itself as sharing strength across sparse
genes, not through prior-knowledge embeddings. Its gene-embedding ablation can stay small.

**12. §5 "Ideas considered": add the modifier-gene approach. (could)** Rosenski et al. 2023 (BMC Med Genomics,
doi:10.1186/s12920-023-01446-6, abstract seen) selects a few modifier genes per target with an ensemble of linear and
non-linear tests, and reports about 3,000 well-predicted genes. Cite it as precedent for M3's inner-fold feature
screening. It is not a separate rung.

No change needed:
- PatientID grouping. No reviewed paper does it, so the gap is real, though leakage is probably small in practice.
- Excluding STRING and ORCS.
- The threshold-based labels, which avoid the circular probability-of-dependency label.
- The plan to harmonise screens first and the anchor-based ceiling (§7.2).
