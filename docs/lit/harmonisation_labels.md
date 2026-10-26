# Harmonisation of CRISPR screens and essentiality labels

Scope: how to make an external, published CRISPR screen comparable to DepMap 26Q3 Chronos gene
effects, and how to define essentiality labels/metrics. All claims carry a DOI or URL; items I
could not verify are marked `[unverified]`. Compiled 2026-10-26.

## 1. Scoring algorithms and their comparability

| Method | Output | Units / scale | Comparable to Chronos gene effect? |
|---|---|---|---|
| Chronos [1] | gene effect (relative growth-rate change) | ~0 = no effect, −1 ≈ median common essential after scaling | yes, by construction |
| CERES [2] | gene effect | same convention after per-line scaling | approximately, after scaling |
| BAGEL / BAGEL2 [3,4] | log Bayes factor | monotone evidence score, not an effect size | rank only |
| MAGeCK RRA [5] | gene-level p/FDR/rank per contrast | rank only | rank only |
| MAGeCK-MLE (VISPR) [6] | beta score per condition | log-scale selection coefficient | rank, and loosely scale-matched |
| JACKS [7] | gene effect with guide-efficacy priors | relative to the screen set analysed jointly | rank; needs many screens |
| CRISPRcleanR [8] | corrected sgRNA counts and logFC | raw LFC units | pre-processing step, not a scorer |
| drugZ [9] | normZ for a drug-vs-control comparison | chemogenetic interaction, **not fitness** | no — exclude |

Practical consequences. (i) BAGEL BF and MAGeCK FDR are *evidence*, not effect size: they saturate,
depend on the reference training sets, and are not linear in fitness, so only rank-based or
re-thresholded use is defensible. (ii) CRISPRcleanR corrects gene-independent (mostly
copy-number-driven) responses by unsupervised segmentation of sgRNA logFC along the genome and emits
corrected counts, so it composes with BAGEL2 (which has no CN correction of its own [4]) and with
Chronos. (iii) Scale alignment is cheap and standard: Chronos's authors normalised CERES and BAGEL2
outputs so the per-line median of common essentials is −1 [1], following [2].

### Can Chronos be run on our external read counts?

Yes, and this is the single highest-value harmonisation move available to us. The code is
open source under BSD-3-Clause (github.com/broadinstitute/chronos [10]).

- **Inputs** (three data frames): a readcount matrix (sequenced entities × sgRNA); a *sequence map*
  with `sequence_id`, `cell_line_name` (must be `"pDNA"` for plasmid samples), `pDNA_batch`, and
  `days` post-infection; and a *guide map* with `sgrna` and `gene`, each sgRNA mapping to exactly one
  gene [10]. Optional but strongly recommended: `negative_control_sgrnas`, which enable readcount
  normalisation and per-screen estimation of the negative-binomial overdispersion [10].
- **Minimum design**: pDNA counts plus at least one later time point. The paper evaluates the
  single-screen regime explicitly — Chronos was run on each Achilles screen individually and
  outperformed competitors by NNMD, whereas "CERES is not designed to run on a single cell line" [1].
  Extra time points help (benchmarked on the 7-timepoint HT-29 series of [11]).
- **Transferring DepMap parameters**: initialise `Chronos(..., pretrained=True)` and call
  `model.import_model(dir)`, where `dir` holds `gene_effect`, `guide_efficacy.csv`, `t0_offset.csv`,
  `library_effect.csv` (and optionally `screen_delay.csv`); `chronos.fetch_parameters()` downloads
  the trained DepMap release parameters [10]. Caveat: the pretrained route only covers DepMap public
  libraries — README states "currently just Avana and KY" — and `load()` raises if the new data
  contain a library absent from the pretrained `library_effect` [10]. So TKOv1/v3, Brunello, GeCKO
  screens must be run *de novo* (pretrained=False), not through the DepMap parameter transfer.
- **Post-processing**: `chronos.alternate_cn(gene_effect, copy_number)` for CN correction (expects
  log2(relative CN + 1)); then scale so median common essentials = −1 and median nonessentials = 0.
  The README explicitly warns *against* per-cell-line normalisation of Chronos output [10].
- **Runtime (CPU)**: the packaged vignette trains in ~10 minutes on a laptop; the full Achilles
  dataset takes 3–4 hours on a 52 GB VM [10]. A single genome-wide external screen is therefore well
  inside our 12-core/16 GB budget.

## 2. Cross-dataset integration: what the batch effects actually are

**Broad vs Sanger (Dempster 2019 [12]).** 147 cell lines screened at both institutes. Gene means
agreed (Pearson 0.784 processed, 0.818 unprocessed). "Batch effects are driven principally by two
key experimental parameters: the reagent library and the assay length" [12]. Concretely: Broad used
21 days vs Sanger 14 days, and *late* dependencies are deeper at the Broad (median −0.402 vs
−0.269) while *early* dependencies match (−0.830 vs −0.781) [12] — i.e. screen length acts as a
gene-set-specific, not global, shift. Guide efficacy differences between libraries explain
individual gene disagreements (median estimated sgRNA efficacy, e.g. EIF3F, MDM2) [12]. ComBat on
unprocessed gene scores aligned gene means and variances (post-correction gene-mean r = 0.9997,
gene-SD r = 0.957) and raised recovery of a line's counterpart in the other dataset from nAUC 0.53
to 0.91 [12]. Screen quality (TPR for common essentials at 5% FDR calibrated on the nonessential
distribution) was a strong predictor of residual disagreement [12].

**Integrated Broad–Sanger (Pacini 2021 [13]).** 908 lines, 17,486 genes, 168 lines overlapping
across 16 tissues. Pipeline: per-screen pre-processing (CRISPRcleanR, or CCR+JACKS, or CERES) →
quantile normalisation → screen-quality adjustment by spline → **ComBat estimated on the overlapping
lines only and then applied to all screens** → second quantile normalisation → removal of PC1, which
the authors attribute to "technical confounders including assay length, guide library and media
conditions" [13]. Same-cell-line recall at k = 1 reached 90% with CERES input. Downsampling showed
as few as 8 shared lines already reduce batch effects, 30–40 give r > 0.995 with the full-set
correction, and convergence at ~90 [13]. **This is the key number for us: a handful of anchor lines
suffices to fit a usable gene-wise correction.**

**DepMap's own integration (26Q3-era pipeline [14]).** Achilles = Avana (Cas9) + Humagne-CD (Cas12),
combined with Sanger Project SCORE (KY Cas9). Chronos is run per library/screen-type to give
`ScreenGeneEffect`, then jointly across all library-screen types to give `CRISPRGeneEffect` using
"Chronos's innate batch correction"; CN correction is applied, then scaling to median common
essential = −1 / nonessential = 0; chromosome-arm median alignment corrects proximity bias [14].
Release notes: 26Q1 extended Chronos library correction from genes present in *all* libraries to
genes present in *more than one* screen batch [15]; 25Q3 split NovaSeq vs HiSeq Humagne screens into
separate Chronos batches [16]. Screen QC gates are explicit and worth reusing on external screens:
>70% Cas9/Cas12a activity, >185 mean reads per guide, screen-level NNMD < −1.25, cross-library read
fraction < 0.1 [14].

**Library/CN/Cas9 confounders.** Amplified regions produce gene-independent depletion [17,18], so
CRISPRcleanR [8] or Chronos's CN correction [10] is mandatory for aneuploid external lines. Library
generations differ markedly in essential/nonessential separation (GeCKOv2 ≪ Avana ≈ Brunello ≈
TKOv3) [19,20], so GeCKO-era screens have compressed dynamic range. Cas9 activity differences are a
per-screen quality term, handled by the NNMD gate rather than by rescaling.

**TP53 status.** Cas9 cutting activates p53 and selects against p53-competent cells, shown in
hTERT-RPE1 [22] and in hPSCs, where DSBs are p53-dependently toxic [23]; Cas9 expression itself
selects for p53-inactivating mutations [24]. Our external targets (RPE1, iPSC) are mostly TP53-WT
whereas most DepMap lines are not, so screen-level signal-to-noise can be lower and DDR/p53-pathway
genes can behave as apparent proliferation suppressors. Genome-scale screens *are* nonetheless
feasible in TP53-WT cells [25] — this is a quality/offset issue, not a disqualification.
Chromosome-arm-scale truncations after cutting [26] add a position-dependent artefact, which DepMap
handles by arm-median alignment [14].

## 3. Harmonising one external screen against DepMap

Ordered from weakest to strongest assumption; the first three preserve across-line variation needed
for per-gene correlation, the last two do not.

1. **Control-based affine scaling (default).** Within the external screen, compute median score of
   the common-essential control set (`m_E`) and of the nonessential set (`m_N`), then map
   `x → (x − m_N) / (m_N − m_E) × (−1)`. This is DepMap's own convention [14,2] and is applied by
   the Chronos authors to foreign-pipeline outputs [1]. It is monotone and linear, so it preserves
   both gene ranks and relative across-line differences; it corrects global depth/screen-length
   scale but not gene-specific library effects.
2. **Gene-wise ComBat using anchor lines (recommended addition).** Fit ComBat [27] with "study" as
   batch on the union of DepMap and external profiles, estimating gene-wise location/scale shifts
   from lines screened in *both* [13]. Preserves within-gene across-line variance up to an affine
   per-gene map, which is exactly the invariance per-gene Pearson correlation needs. Requires ≥8,
   preferably ≥30 anchors [13].
3. **Quantile normalisation to the reference distribution.** Maps each external profile onto the
   DepMap score distribution. Monotone *within* a profile; it therefore preserves per-line gene
   ranks and AUPRC exactly, and mostly preserves across-line ordering per gene, but it distorts
   magnitudes non-linearly, so per-gene correlations after QN are closer to rank correlations. Use
   as a secondary variant, as in [13].
4. **Per-screen z-scoring.** Removes screen-level mean/scale but cannot fix gene-specific library
   bias; acceptable only when control sets are missing or degenerate.
5. **Rank transformation per screen.** Robust to any monotone scoring (BAGEL BF, MAGeCK beta,
   normZ), and the only honest option for BF/FDR-only screens — but it destroys effect magnitude, so
   evaluate it with Spearman and AUPRC only, never with Pearson per gene.

Decision rule for our pipeline: **if raw counts + a pDNA/T0 sample exist, re-score with Chronos
(de novo for non-Avana/KY libraries) → CN-correct → control-scale → gene-wise ComBat on anchors.**
If only gene-level BF/MAGeCK output exists, use control-scaling where control sets are present, else
rank-transform and restrict the evaluation to rank metrics.

### Anchor lines screened both in DepMap and in external studies

Verified against the local `data/Model.csv` + `CRISPRGeneEffect.csv` (26Q3 assumption per DATA.md):

| Line | DepMap ModelID | In 26Q3 gene effect | External screen (library, duration) |
|---|---|---|---|
| hTERT-RPE1 (SS48, SS77 near-diploid controls) | ACH-002462, ACH-002463 | yes | Hart 2015 TKOv1, 18 d [28]; Hart 2017 TKOv1/TKOv3 incl. RPE1 and RPE1-dTP53 [29]; Brown 2019 TP53-WT/null [25] |
| HAP1 | ACH-002475 | yes | Hart 2017 TKOv1/TKOv3 [29] |
| HCT 116 | ACH-000971 | yes | Hart 2015 TKOv1 [28]; Hart 2017 [29] |
| DLD-1 | ACH-001061 | yes | Hart 2015 TKOv1 [28]; Hart 2017 [29] |
| HT-29 | ACH-000552 | yes | Tzelepis 2016 KY, 7-timepoint series [11] |
| MOLM-13, MV4-11, OCI-AML2, OCI-AML3 | ACH-000362, ACH-000045, ACH-000113, ACH-000336 | yes | Tzelepis 2016 KY [11] |
| A-375 | ACH-000219 | yes | Sanson 2018 Brunello vs GeCKOv2/Avana [20] |
| K-562 | ACH-000551 | yes | Wang 2015 Sabatini library [30] |
| HEK-293, HeLa, MCF 10A, BJ hTERT | ACH-001085, ACH-001086, ACH-001357, ACH-001018 | **no** | — (not usable as anchors) |

Caution: HEK293T and HeLa are *not* in the 26Q3 gene-effect matrix, so the often-quoted
"HEK293T/HeLa in both datasets" anchors do not exist for us. RPE1 anchors are the engineered
aneuploidy panel (near-diploid SS48/SS77 vs trisomic SS6/SS119/SS51), not parental hTERT-RPE1, and
they are Avana/Chronos — a genotype caveat for comparison against TKO RPE1 screens. Also available
are the 168 Broad/Sanger overlapping lines already harmonised in [13] and the cross-library
aggregations in PICKLES v3 (Avana + Score + TKOv3) [31] and iCSDB [32], usable as independent
cross-checks.

## 4. Essentiality labels and metrics

- **Gene-effect thresholds.** Conventional cut-offs on Chronos scores are −0.5 ("likely dependent")
  and −1 (≈ median common essential). −1 is anchored by the scaling convention [14,2]; −0.5 is a
  community convention, not a DepMap definition `[unverified as an official threshold]`.
- **Probability of dependency** (`CRISPRGeneDependency`): per gene × line, estimated by EM; the
  dependent component is taken from the essential-gene list and the null from unexpressed genes
  where expression data exist, else from the nonessential list [14]. P > 0.5 is the usual call.
  Chronos exposes the same idea via `get_probability_dependent()` and an empirical-p alternative
  whose resolution is capped at 1/|negative controls| [10] — relevant if an external screen is
  sub-genome-scale.
- **Common essentials.** `CRISPRInferredCommonEssentials` follows DepMap's definition: a gene that
  ranks in the top X most-depleting genes in ≥90% of lines, with X set empirically from the minimum
  of the distribution of gene ranks in each gene's 90th-percentile least-dependent line [12,14].
  (An alternative, benchmarked implementation of core-fitness-gene calling is CoRe [21].)
  The *control* sets used for scaling are different objects: `AchillesCommonEssentialControls` =
  intersection of the Hart [33] and Blomen [34] essentials; `AchillesNonessentialControls` = the Hart
  reference nonessentials [14].
- **CEGv2 / NEGv1.** CEGv2 = 684 core essential genes from BAGEL re-analysis of 17 screens
  (Hart 2017 G3) [29]; NEGv1 = reference nonessentials from Hart 2014 [33]. Download from
  github.com/hart-lab/bagel (`CEGv2.txt`, `NEGv1.txt`; 684 and 928 gene rows in the current copy),
  repository licence **MIT** [35]. These are independent of DepMap's scaling controls and so are the
  right *held-out* positive/negative sets for evaluating harmonisation.
- **How to compute AUPRC.** Three different questions, which must not be conflated:
  (a) **per line, pooled over genes** — label = CEGv2 vs NEGv1 (or dependency P > 0.5), score =
  −gene effect. This is the standard screen-quality/essentiality AUPRC and the one used in library
  and algorithm benchmarks [4,29,1]; (b) **per gene, across lines** — label = "dependent line",
  measures selective-dependency transfer, the quantity our model is actually for; (c) **pooled over
  all gene × line pairs** — dominated by common essentials and therefore optimistic. Report (a) and
  (b) separately; never report (c) alone.
- **Separating core-essential from context-specific transfer.** Stratify every metric by gene class:
  (i) CEGv2 / DepMap common essentials — expected to transfer and so a *harmonisation control*, not
  evidence of model skill; (ii) genes with high across-line variance in DepMap (e.g. SD > 0.2 and
  ≥5 dependent lines) — the real test; (iii) NEGv1/unexpressed genes — false-positive control.
  Per-gene correlation must be computed on (ii), because on (i) the across-line variance is mostly
  noise, and screen-quality differences alone can manufacture apparent correlation.

## 5. Recommendations for our design

**Harmonisation protocol (per external screen).**
1. Prefer raw counts: reconstruct a sequence map (`days`, `pDNA_batch`, pDNA/T0) and guide map, and
   re-score with Chronos de novo [10]. Record library, duration, Cas9 delivery and TP53 status.
2. CN-correct (`alternate_cn`, or CRISPRcleanR [8] upstream for non-Chronos routes).
3. Gate on DepMap's QC criteria (NNMD < −1.25, reads/guide, cross-library contamination) [14];
   record, do not silently drop.
4. Affine-scale to median common-essential controls = −1, nonessentials = 0 [14,2].
5. Fit gene-wise ComBat [27] on the §3 anchors and apply to all external screens [13]; keep the
   un-ComBat'ed version as an ablation, plus quantile-normalised and rank-only variants.
6. Evaluate only on genes assayed by both libraries (≥3 guides externally); report coverage.

**Labels.** Primary: DepMap dependency probability > 0.5 for the binary label and raw Chronos score
for regression; secondary thresholds −0.5 and −1 reported as sensitivity analysis. Gene strata: CEGv2
positives, NEGv1 negatives, DepMap common essentials, and selective genes as defined above.

**Controls that demonstrate harmonisation works.**
- *Positive*: CEGv2-vs-NEGv1 AUPRC and NNMD in the harmonised external screen should approach the
  DepMap distribution; essential-gene separation must not be created by the correction (so use held
  out CEG/NEG genes, not the control genes used for scaling).
- *Negative*: NEGv1 and unexpressed genes must stay centred at 0 after harmonisation; the correction
  must not shift them.
- *Round-trip anchor test*: re-score an externally published screen of a DepMap-screened line (HAP1,
  HCT116, DLD-1, HT-29, MOLM-13/MV4-11/OCI-AML2/OCI-AML3, A-375) through our pipeline and check the
  profile's nearest neighbour in DepMap is itself — the same-cell-line recall / k-NN metric of
  [12,13], with nAUC reported before and after correction.
- *Negative control for the anchor test*: profiles should **not** match unrelated lines; and a
  permutation of anchor labels must destroy the ComBat benefit.
- *Screen-length control*: early vs late common dependencies should show the [12] signature before
  correction and reduced separation after.

## References

1. Dempster JM et al. Chronos: a cell population dynamics model of CRISPR experiments that improves inference of gene fitness effects. *Genome Biol* 2021;22:343. doi:10.1186/s13059-021-02540-7
2. Meyers RM et al. Computational correction of copy number effect improves specificity of CRISPR-Cas9 essentiality screens in cancer cells. *Nat Genet* 2017;49:1779–1784. doi:10.1038/ng.3984
3. Hart T, Moffat J. BAGEL: a computational framework for identifying essential genes from pooled library screens. *BMC Bioinformatics* 2016;17:164. doi:10.1186/s12859-016-1015-8
4. Kim E, Hart T. Improved analysis of CRISPR fitness screens and reduced off-target effects with the BAGEL2 gene essentiality classifier. *Genome Med* 2021;13:2. doi:10.1186/s13073-020-00809-3
5. Li W et al. MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens. *Genome Biol* 2014;15:554. doi:10.1186/s13059-014-0554-4
6. Li W et al. Quality control, modeling, and visualization of CRISPR screens with MAGeCK-VISPR. *Genome Biol* 2015;16:281. doi:10.1186/s13059-015-0843-6 (see also MAGeCKFlute, *Nat Protoc* 2019, doi:10.1038/s41596-018-0113-7)
7. Allen F et al. JACKS: joint analysis of CRISPR/Cas9 knockout screens. *Genome Res* 2019;29:464–471. doi:10.1101/gr.238923.118
8. Iorio F et al. Unsupervised correction of gene-independent cell responses to CRISPR-Cas9 targeting. *BMC Genomics* 2018;19:604. doi:10.1186/s12864-018-4989-y
9. Colic M et al. Identifying chemogenetic interactions from CRISPR screens with drugZ. *Genome Med* 2019;11:52. doi:10.1186/s13073-019-0665-3
10. Chronos source and README/Vignette, github.com/broadinstitute/chronos (BSD-3-Clause; inputs, `pretrained`/`import_model`, `fetch_parameters`, `alternate_cn`, runtime statements) — inspected 2026-10-26
11. Tzelepis K et al. A CRISPR dropout screen identifies genetic vulnerabilities and therapeutic targets in acute myeloid leukemia. *Cell Rep* 2016;17:1193–1205. doi:10.1016/j.celrep.2016.09.079
12. Dempster JM et al. Agreement between two large pan-cancer CRISPR-Cas9 gene dependency data sets. *Nat Commun* 2019;10:5817. doi:10.1038/s41467-019-13805-y
13. Pacini C et al. Integrated cross-study datasets of genetic dependencies in cancer. *Nat Commun* 2021;12:1661. doi:10.1038/s41467-021-21898-7
14. DepMap CRISPR pipeline and analysis documentation. https://forum.depmap.org/t/crispr-pipeline-and-analysis/3393
15. DepMap 26Q1 release notes. https://forum.depmap.org/t/announcing-the-26q1-release/4606
16. DepMap 25Q3 release notes. https://forum.depmap.org/t/announcing-the-25q3-release/4476
17. Aguirre AJ et al. Genomic copy number dictates a gene-independent cell response to CRISPR/Cas9 targeting. *Cancer Discov* 2016;6:914–929. doi:10.1158/2159-8290.CD-16-0154
18. Munoz DM et al. CRISPR screens provide a comprehensive assessment of cancer vulnerabilities but generate false-positive hits for highly amplified genomic regions. *Cancer Discov* 2016. doi:10.1158/2159-8290.CD-16-0178
19. Doench JG et al. Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9. *Nat Biotechnol* 2016;34:184–191. doi:10.1038/nbt.3437
20. Sanson KR et al. Optimized libraries for CRISPR-Cas9 genetic screens with multiple modalities. *Nat Commun* 2018;9:5416. doi:10.1038/s41467-018-07901-8
21. Vinceti A et al. CoRe: a robustly benchmarked R package for identifying core-fitness genes. *BMC Genomics* 2021. doi:10.1186/s12864-021-08129-5
22. Haapaniemi E et al. CRISPR-Cas9 genome editing induces a p53-mediated DNA damage response. *Nat Med* 2018;24:927–930. doi:10.1038/s41591-018-0049-z
23. Ihry RJ et al. p53 inhibits CRISPR-Cas9 engineering in human pluripotent stem cells. *Nat Med* 2018;24:939–946. doi:10.1038/s41591-018-0050-6
24. Enache OM et al. Cas9 activates the p53 pathway and selects for p53-inactivating mutations. *Nat Genet* 2020;52:662–668. doi:10.1038/s41588-020-0623-4
25. Brown KR, Mair B, Soste M, Moffat J. CRISPR screens are feasible in TP53 wild-type cells. *Mol Syst Biol* 2019;15:e8679. doi:10.15252/msb.20188679
26. Lazar NH et al. High-resolution genome-wide mapping of chromosome-arm-scale truncations induced by CRISPR-Cas9 editing. *Nat Genet* 2024;56:1482–1493. doi:10.1038/s41588-024-01758-y
27. Johnson WE, Li C, Rabinovic A. Adjusting batch effects in microarray expression data using empirical Bayes methods. *Biostatistics* 2007;8:118–127. doi:10.1093/biostatistics/kxj037
28. Hart T et al. High-resolution CRISPR screens reveal fitness genes and genotype-specific cancer liabilities. *Cell* 2015;163:1515–1526. doi:10.1016/j.cell.2015.11.015 (TKOv1, 18 d; per-screen metadata via BioGRID ORCS screens 16/17/27/28/29)
29. Hart T et al. Evaluation and design of genome-wide CRISPR/SpCas9 knockout screens. *G3* 2017;7:2719–2727. doi:10.1534/g3.117.041277 (CEG2 = 684 genes; TKOv1 screens in DLD1, GBM, HAP1, HCT116, RPE1, RPE1-dTP53)
30. Wang T et al. Identification and characterization of essential genes in the human genome. *Science* 2015;350:1096–1101. doi:10.1126/science.aac7041
31. Novak LC et al. PICKLES v3: the updated database of pooled in vitro CRISPR knockout library essentiality screens. *Nucleic Acids Res* 2023. doi:10.1093/nar/gkac982
32. Choi A et al. iCSDB: an integrated database of CRISPR screens. *Nucleic Acids Res* 2021. doi:10.1093/nar/gkaa989
33. Hart T et al. Measuring error rates in genomic perturbation screens: gold standards for human functional genomics. *Mol Syst Biol* 2014;10:733. doi:10.15252/msb.20145216
34. Blomen VA et al. Gene essentiality and synthetic lethality in haploid human cells. *Science* 2015;350:1092–1096. doi:10.1126/science.aac7557
35. hart-lab/bagel repository (CEGv2.txt, NEGv1.txt; MIT licence). https://github.com/hart-lab/bagel
