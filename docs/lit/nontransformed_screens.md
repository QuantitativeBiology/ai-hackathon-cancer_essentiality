# CRISPR loss-of-function screens in non-transformed human cells

Catalogue compiled 2026-10-26 for the out-of-domain generalisation test of the DepMap 26Q3
essentiality model. Scope: published genome-wide (or large) CRISPR knockout screens
(CRISPRi secondary) in hTERT-immortalised lines, hESC/iPSC, iPSC-derived lineages and
primary haematopoietic cells, together with a matched bulk expression profile.
Items I could not confirm from a primary or curated source are marked `[unverified]`.

## Summary table

| # | Screen (first author, year) | Cell model | Category | TP53 | Library | Cas9/dCas9 delivery | Duration | Readout | Data released / where | Matched expression |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | Olivieri 2020, *Cell* | RPE1-hTERT *TP53*-/- Cas9 | hTERT line | KO | TKOv3 (genome-wide) | stable lentiviral Cas9 | ~18 d (T0→T18) | dropout, incl. untreated arms | per-guide read counts + supp raw data, Mendeley `10.17632/gfcn2wmrpf` (v2) | RPE1-SS sublines in DepMap 26Q3 |
| 2 | Brown 2019, *Mol Syst Biol* | RPE1-hTERT *TP53*-WT **and** *TP53*-null | hTERT line | both | genome-wide (TKOv3) | stable lentiviral Cas9 | `[unverified]` | dropout | GEO **GSE128210** (raw screen seq) | RPE1-SS sublines in DepMap 26Q3 |
| 3 | Hart 2015, *Cell* | hTERT-RPE1 (+ HCT116, HeLa, DLD1, GBM) | hTERT line | WT | TKOv1 (~90k guides) | stable lentiviral Cas9, MOI ~0.3 | 18 d timecourse | dropout, BAGEL BF | gene-level BF via BioGRID ORCS screen **28** (RPE1), **16**/**27** (HCT116); raw counts in *Cell* supp. tables `[unverified]` | RPE1-SS sublines in DepMap 26Q3 |
| 4 | Mair 2019, *Cell Rep* | H1 hESC, inducible Cas9, on feeders vs laminin | ESC | WT | TKOv3 | dox-inducible Cas9 (integrated) | `[unverified]` (~18–20 d) | dropout | Mendeley `10.17632/rhsw3nbjhz.1`; GEO **GSE156284** (CRISPRseq) | ENCODE/Roadmap H1 RNA-seq |
| 5 | Yilmaz 2018, *Nat Cell Biol* | h-pES10 haploid hESC | ESC | WT `[unverified]` | genome-wide, 18,166 genes | `[unverified]` | 23 d | dropout, FDR<0.05 | gene-level hits via BioGRID ORCS screen **597** / dataset **76**; follow-up RNA-seq GEO **GSE103846** | GSE103846 controls (same line); ENCODE hESC as proxy |
| 6 | Ihry 2019, *Cell Rep* | hPSC (H1 + iPSC `[unverified]`), inducible Cas9 | ESC/iPSC | WT | genome-scale (Novartis) | dox-inducible Cas9 | `[unverified]` | self-renewal/survival dropout + differentiation FACS | location `[unverified]` — check STAR Methods | ENCODE H1; HipSci for iPSC |
| 7 | Tian 2021, *Nat Neurosci* | WTC11 iPSC-derived glutamatergic neurons | iPSC-derived lineage | WT | genome-wide CRISPRi **and** CRISPRa (hCRISPRi/a-v2) | dox-inducible dCas9-KRAB / -VPH | ~10–14 d post-differentiation `[unverified]` | survival dropout | CRISPRbrain data commons (crisprbrain.org) | GEO for WTC11 neurons `[unverified]`; HipSci/Allen WTC11 profiles |
| 8 | Dräger 2022, *Nat Neurosci* | iPSC-derived microglia (iTF-microglia) | iPSC-derived lineage | WT | CRISPRi/a sublibraries `[unverified]` | dox-inducible dCas9-KRAB/-VPH | `[unverified]` | survival + phagocytosis FACS + scRNA | CRISPRbrain; PMC9448678 supp. | PMC9448678 / GEO `[unverified]` |
| 9 | Tian 2019, *Neuron* | WTC11 iPSC-derived neurons | iPSC-derived lineage | WT | CRISPRi "H1" sublibrary (~2,300 genes) | dox-inducible dCas9-KRAB | ~10 d `[unverified]` | survival dropout + scRNA + imaging | GEO **GSE124703** | as #7 |
| 10 | Shifrut 2018, *Cell* | primary human CD8+ T cells (healthy donors) | primary haematopoietic | WT | Brunello (genome-wide) | Cas9 protein electroporation (SLICE) | ~6 d (+4 d drug arm) | proliferation dropout / stimulation enrichment | gene-level via BioGRID ORCS screens **1107**, **1109**, **1110**; GEO accession `[unverified]` | donor T-cell RNA-seq in same paper `[unverified]`; DICE / Blueprint reference profiles |
| 11 | Ting 2018, *Nat Methods* | primary CD4+ T cells; primary CD34+ HSPC | primary haematopoietic | WT | genome-scale (Guide Swap) | Cas9 RNP + non-targeting guide, lentiviral sgRNA | `[unverified]` | depletion/enrichment; HSPC expansion | location `[unverified]` | Blueprint/HSPC reference RNA-seq |
| 12 | Bowden 2020, *eLife* | RPE1 *TP53*-WT vs *TP53*-KO, 2 clones each | hTERT line | both | custom DDR dual-sgRNA, 852 genes + 95 CEGs + OR controls | lentiviral Cas9 (monoclonal) | T0 d3 → d15, d19 | dropout, MAGeCK | Dryad `10.5061/dryad.2fqz612kr`; supp file 2 = guide counts, file 3 = MAGeCK gene scores | RPE1-SS sublines in DepMap 26Q3 |
| 13 | Leng 2022, *Nat Neurosci* | iPSC-derived astrocytes | iPSC-derived lineage | WT | CRISPRi sublibraries | dox-inducible dCas9-KRAB | n/a | **FACS/scRNA reactive-state**, not fitness | PMC9633461 supp.; CRISPRbrain `[unverified]` | n/a |
| 14 | Hart 2017, *G3* | TKOv3 design + reanalysis of 17 screens (CEG2, 684 genes) | reference, not a new non-transformed screen | — | TKOv3 | — | — | — | *G3* supplementary tables | — |

## Notes per screen

**1. Olivieri M, et al. (2020) "A genetic map of the response to DNA damage in human cells."
*Cell* 182(2):481–496.e21. doi:10.1016/j.cell.2020.05.040.** The single largest coherent body of
genome-wide TKOv3 screen data in an hTERT line. ~30 screens in RPE1-hTERT *TP53*-/- Cas9 against
27 genotoxins; the **untreated/mock arms are the usable fitness screens** for us. Mendeley deposit
contains a folder "CRISPR screen readcount files" (per-guide counts) plus "Supplementary Raw Data".
Caveat: *TP53* is knocked out, so p53-pathway dependencies are not interrogable; and the line is
pseudo-diploid but karyotypically drifted.

**2. Brown KR, Mair B, Soste M, Moffat J (2019) "CRISPR screens are feasible in TP53 wild-type
cells." *Mol Syst Biol* 15(8):e8679. doi:10.15252/msb.20188679.** Genome-wide screens in RPE1
with *TP53* WT and *TP53*-null backgrounds; GEO **GSE128210** confirmed (title matches). Bowden 2020
re-used five of its wild-type RPE1 screens. This is the best *TP53*-intact genome-wide hTERT-line
dataset and the natural control for #1.

**3. Hart T, et al. (2015) "High-resolution CRISPR screens reveal fitness genes and
genotype-specific cancer liabilities." *Cell* 163(6):1515–1526. doi:10.1016/j.cell.2015.11.015
(PMID 26627737).** The canonical hTERT-RPE1 TKOv1 screen. BioGRID ORCS curates it as hTERT-RPE1,
TKO v1, Cas9 knockout, 18-day timecourse, MOI ~0.3, BAGEL, significance BF > 6.84. Gene-level
Bayes factors are downloadable from ORCS; I could **not** verify a GEO/SRA accession for raw counts,
only the journal supplement — treat per-guide counts as `[unverified]`.

**4. Mair B, et al. (2019) "Essential gene profiles for human pluripotent stem cells identify
uncharacterized genes and substrate dependencies." *Cell Rep* 27(2):599–615.e12.
doi:10.1016/j.celrep.2019.02.041.** H1 hESC, inducible Cas9, TKOv3, two substrates (feeders and
laminin/feeder-free — the paper and the deposit word this differently). Mendeley
`10.17632/rhsw3nbjhz.1` holds raw data for all figures; GEO **GSE156284** is the CRISPRseq series.
Reported hPSC essential genes diverge from cancer-line essentials and depend on substrate — exactly
the context effect our test should detect.

**5. Yilmaz A, Peretz M, Aharony A, Sagi I, Benvenisty N (2018) "Defining essential genes for human
pluripotent stem cells by CRISPR-Cas9 screening in haploid cells." *Nat Cell Biol* 20(5):610–619.
doi:10.1038/s41556-018-0088-1.** Haploid hESC line h-pES10, library covering 18,166 genes, 23-day
screen, 2,553 hits at FDR<0.05 (ORCS screen 597). Haploidy removes buffering by the second allele,
so effect sizes are not on the same scale as diploid screens — a strong argument for rank-based
rather than absolute-score evaluation. GEO **GSE103846** is the siRNA-knockdown transcriptome
follow-up (SALL4, DSCC1, SEPHS1, VRTN), *not* the screen itself, but its controls give a same-line
expression profile.

**6. Ihry RJ, et al. (2019) "Genome-scale CRISPR screens identify human pluripotency-specific
genes." *Cell Rep* 27(2):616–630.e6. doi:10.1016/j.celrep.2019.03.043.** (The task sheet listed
this as *Genome Research*; it is *Cell Reports*, same issue as Mair.) Three genome-wide screens:
self-renewal/survival, survival after single-cell dissociation, differentiation capacity. The
survival screen is the usable fitness arm. Deposit location unresolved from public indexes —
needs the STAR Methods key-resources table.

**7–9, 13. Kampmann lab / CRISPRbrain.** Tian R, et al. (2021) *Nat Neurosci* 24(7):1020–1034,
doi:10.1038/s41593-021-00862-0 is the only **genome-wide** screen in an iPSC-derived lineage I could
confirm, and the survival arm is a true fitness screen. It is **CRISPRi**, so gene-level scores are
knockdown, not knockout — direction is comparable, magnitude is not. Tian 2019 (*Neuron*
104(2):239–255.e12) is a ~2,300-gene sublibrary; Dräger 2022 (*Nat Neurosci* 25(9):1149–1162,
doi:10.1038/s41593-022-01131-4) is microglia CRISPRi/a; Leng 2022 (*Nat Neurosci*
25(11):1528–1542) is a reactive-state FACS screen and is **not** usable as a fitness test.
All are indexed in **CRISPRbrain** (crisprbrain.org), a data commons for differentiated human cell
types; its programmatic/bulk download endpoint and licence are `[unverified]` (the site is a
JavaScript app and `/api/screens` returned 401), so budget manual export or a direct request.

**10. Shifrut E, et al. (2018) "Genome-wide CRISPR screens in primary human T cells reveal key
regulators of immune function." *Cell* 175(7):1958–1971.e15. doi:10.1016/j.cell.2018.10.024.**
SLICE: lentiviral sgRNA + Cas9 protein electroporation in primary CD8+ T cells, Brunello library,
MAGeCK; ORCS curates a 6-day timecourse screen and a 4-day drug-exposure screen. The short duration
means few population doublings, so only strong essentials will drop out — expect compressed dynamic
range relative to DepMap's 21-day Avana screens.

**11. Ting PY, et al. (2018) "Guide Swap enables genome-scale pooled CRISPR-Cas9 screening in human
primary cells." *Nat Methods* 15(11):941–946. doi:10.1038/s41592-018-0149-1.** Only confirmed
genome-scale screen touching primary CD34+ HSPC (expansion phenotype) plus CD4+ T cells. Data
availability `[unverified]`.

**Not found / negative results.** I found **no** genome-wide knockout fitness screen in BJ-hTERT,
HME1-hTERT, MCF10A or IMR-90 — targeted searches returned only the Hart/TKO reference papers.
Published iPSC-derived **cardiomyocyte** screens (Sapp et al. 2021 *Sci Rep*
doi:10.1038/s41598-021-92988-1, doxorubicin toxicity) and **endothelial** screens (protocols.io
WTC11Ci CRISPRi differentiation screens, CDH5 FACS) are drug-response or differentiation screens,
not baseline fitness; no hepatocyte screen found. Sanger/Wellcome **HipSci** genome-wide iPSC
fitness screens are an announced project ("Parts HipSci") and a Sanger thesis reporting
high replicate-to-replicate variability attributed to Cas9 DSB toxicity, but I found **no published
paper or accession** `[unverified]`. I found no genome-wide KO fitness screen in EBV-transformed
lymphoblastoid lines; note they are virally transformed and would not satisfy the
"non-transformed" criterion anyway.

## Aggregators

- **BioGRID ORCS** — version **2.0.18**, **2,217 screens from 418 publications**; downloads at
  `https://downloads.thebiogrid.org/BioGRID-ORCS` (per-organism zip archives of tab-delimited
  files, plus a `LATEST` alias for scripted download), **MIT licence**, academic and commercial use,
  no warranty. Screen records do carry cell line, cell type, library, enzyme, methodology, screen
  duration/timecourse, score type and the authors' significance threshold — i.e. exactly the
  harmonisation covariates we need. **Primary route for screen metadata; often gene-level scores
  only, not per-guide counts.**
- **iCSDB** (Choi et al., *Nucleic Acids Res* 2021;49:D956–D961; kobic.re.kr/icsdb) — 1,375
  genome-wide screens over 976 cell lines, built from DepMap + ORCS, with library batch effects
  removed and scores standardised. Useful as a cross-check on harmonisation, but inherits the
  cancer-line bias of its sources. Live status `[unverified]`.
- **PICKLES v3** (Lenoir et al., *Nucleic Acids Res* 2023) — 1,162 whole-genome screens,
  18,959 genes, Avana + Score + TKOv3; source code downloadable. Includes TKOv3 RPE1/HAP1
  alongside cancer lines.
- **Project Score** (Sanger) — cancer lines only; no use here.
- **GenomeCRISPR** — could not verify current availability `[unverified]`.

## DepMap 26Q3 "Non-Cancerous" models (verified locally)

`data/Model.csv` lists **152** models with `OncotreePrimaryDisease == "Non-Cancerous"`, of which
**17** have a profile in `CRISPRGeneEffect.csv`:

- five RPE1 single-/multi-trisomy sublines (ACH-002462 RPE1-SS48, ACH-002463 RPE1-SS77 — both
  annotated "near-diploid control"; ACH-002464 RPE1-SS6 trisomy 7; ACH-002465 RPE1-SS119
  trisomy 8; ACH-002467 RPE1-SS51 trisomies 7+22);
- HA1E immortalised embryonic kidney (ACH-001310);
- five immortalised/hyperplastic prostate epithelial lines (WPE1-NA22, BPH1, P4E6, SHMAC4, SHMAC5);
- HS729 fibroblast; C8166 and HG3 B-lymphoblast lines; ARH77 and IM9 "immortalized blood";
  one tubulovillous adenoma organoid (ACH-003722).

**15 of the 17 also have RNA-seq** in `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv`
(ARH77 and IM9 do not). In the literature the RPE1-SS panel is used as the aneuploidy-tolerance
series (near-diploid controls vs engineered trisomies) rather than as a "normal" reference, and
HA1E/prostate lines are immortalised but partly oncogene-transformed. **Decision: these 17 must be
held out of training** — they are the only internally-consistent (Avana, same pipeline) non-cancer
benchmark we have, and they give a library-matched sanity check before we face cross-library
external screens. They are not a substitute for the external test set: C8166, HG3, ARH77 and IM9
are virally/malignancy-derived, and the RPE1-SS lines are a single genetic background.

## Ranking by usability

Criteria: genome-wide + downloadable per-guide counts or gene-level scores + a matched expression
profile + knockout (not CRISPRi) + *TP53* status known.

1. **Olivieri 2020** RPE1-hTERT *TP53*-/- — genome-wide TKOv3, per-guide counts, many replicates.
2. **Brown 2019** RPE1 *TP53*-WT and -null — genome-wide, raw data in GEO, p53-matched pair.
3. **Mair 2019** H1 hESC (two substrates) — genome-wide TKOv3, Mendeley + GEO, ENCODE H1 expression.
4. **Hart 2015** hTERT-RPE1 TKOv1 — genome-wide, gene-level BF via ORCS; different library era, good
   for library-harmonisation stress-testing.
5. **Tian 2021** iPSC-derived neurons — genome-wide, but CRISPRi; the only genome-wide
   iPSC-derived-lineage fitness screen.
6. **Yilmaz 2018** haploid hESC — genome-wide; haploid scaling caveat; gene-level only.
7. **Shifrut 2018** primary CD8+ T cells — genome-wide Brunello, but short screen; ORCS gene scores.
8. **Ihry 2019** hPSC — genome-wide, but deposit unresolved.
9. **Bowden 2020** RPE1 ±*TP53* — only 852 genes, but clean per-guide counts and the best controlled
   p53 contrast; use as a focused validation, not a main test.
10. **Ting 2018** CD34+ HSPC / CD4+ T — genome-scale primary cells, data location unresolved.
11. **Dräger 2022** / **Tian 2019** microglia and neuron sublibraries — sublibrary scale, CRISPRi.
12. **Leng 2022** astrocytes — wrong readout; exclude.

### Proposed core test set (9 screens, all 3 required categories)

| Category | Screens |
|---|---|
| hTERT-immortalised | Olivieri 2020 (RPE1 *TP53*-/-); Brown 2019 (RPE1 *TP53*-WT **and** *TP53*-null, counted as two arms); Hart 2015 (hTERT-RPE1, TKOv1) |
| iPSC / ESC | Mair 2019 (H1, feeders + laminin = two conditions); Yilmaz 2018 (haploid hESC); Ihry 2019 (hPSC) |
| iPSC-derived lineage | Tian 2021 (neurons, genome-wide CRISPRi); Dräger 2022 (microglia, CRISPRi) |
| Primary haematopoietic | Shifrut 2018 (primary CD8+ T cells); Ting 2018 (CD34+ HSPC, if the deposit resolves) |

Hold Bowden 2020 back as a confirmatory p53-contrast set, and the 17 DepMap non-cancerous models as
the within-pipeline positive control.

### Harmonisation requirements implied by this catalogue

Libraries span TKOv1, TKOv3, Brunello and two CRISPRi libraries; durations span ~6 d (primary T
cells) to 23 d (haploid hESC); scores span BAGEL BF, MAGeCK, raw LFC and CRISPRbrain phenotype
scores; ploidy spans haploid, near-diploid and trisomic; *TP53* is WT in most but knocked out in the
largest dataset. Any conclusion about transfer must therefore be made on **within-screen ranks**
(per-gene Spearman across models; AUPRC against that screen's own essential/non-essential controls)
rather than on absolute gene-effect values, and the CRISPRi screens must be evaluated separately
from knockout screens.

### Immediate follow-ups

1. Resolve the Ihry 2019 and Ting 2018 deposits from the papers' data-availability sections.
2. Confirm whether Hart 2015 per-guide counts exist outside the *Cell* supplement.
3. Establish how to bulk-export CRISPRbrain and under what licence.
4. Find matched bulk RNA-seq for WTC11-derived neurons/microglia and for the Shifrut donor T cells;
   otherwise fall back on reference profiles (ENCODE, HipSci, DICE/Blueprint), and record that
   substitution as a limitation.
