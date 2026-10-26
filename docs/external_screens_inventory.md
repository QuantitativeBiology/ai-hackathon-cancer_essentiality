# External CRISPR screens and matched expression: inventory

Compiled 2026-10-26. This covers the out-of-domain test set (non-transformed human cells) and the
**anchor** screens (DepMap lines screened with the same external libraries, often by the same labs) that we'll
use to fit library/lab → DepMap harmonisation. Nothing here has been modelled yet.

**Reproduce:**

```
python scripts/fetch_external_screens.py            # downloads every file in configs/external_screens.json, writes sha256 manifest
python scripts/extract_orcs_nontransformed.py       # ORCS subset -> results/external_screens/orcs_*.csv (+ extracted tab files)
python scripts/check_external_screens.py            # sanity -> results/external_screens/sanity.csv
```

- All URLs are in `configs/external_screens.json`.
- Full per-file URL, bytes and sha256 are in `configs/external_screens_manifest.json` (82 files). The prefixes below come from that manifest.
- Data live under `data/external_screens/<deposit>/` and `data/external_expr/<source>/`. Both are git-ignored and treated as untrusted data: nothing in them is executed, and the one code archive offered (Doench 2016 MOESM9) was deliberately not kept.
- The Olivieri files' sha256 equal the hashes Mendeley publishes for them.

**Licence/terms legend:**
- **GEO**: NCBI places no restrictions on use; submitters may hold rights.
- **CC BY 4.0**: Mendeley, and the OA bundles from MSB, Cell Rep 2016, G3, Nat Commun and eLife.
- **Elsevier**: supplementary files from the publisher CDN under the Elsevier user licence. Fine for research use; don't redistribute.
- **Springer**: Springer Nature supplementary files, same caveat.
- **ORCS**: MIT. Allowed as a test-set source, never as a model input.

## 1. Screens obtained from original deposits

The sanity column gives the median gene score of CEGv2 vs NEGv1, then the AUROC (probability that a CEG gene is more depleted than a NEG gene). For count data the score is a crude T_end vs T0 log2FC; for score data it's the published score. Pass means AUROC ≥ 0.75. Full table: `results/external_screens/sanity.csv`.

| id | citation | category | TP53 | library | duration | data level | files (sha256 prefix) | terms | matched expression | status | sanity: CEG / NEG (AUROC) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| olivieri2020_TKOv3_NT (6 untreated arms S02,S04–S07,S12) | Olivieri 2020 *Cell* 182:481 | hTERT line (RPE1) | KO | TKOv3 (71,090 guides, 18,039 genes) | T0→T18 (S02: T15), 2 reps/arm, own T0 per arm | **counts** | Dataset_S2_readcounts.txt 109be820; Screen_Library.csv bf039af9; Zscores_31screens.csv 7f36e572 | CC BY 4.0 | DepMap RPE1-ss sublines | OK | −3.44 to −3.63 / 0.09 to 0.21 (0.97–0.99) |
| olivieri2020_TKOv2_NT (S08–S10) | same | hTERT line | KO | TKOv2 (70,697 guides) | T0→T18 | **counts** | Dataset_S1_readcounts.txt 4ecbf9a8 | CC BY 4.0 | DepMap RPE1 | OK | −2.15 to −3.30 / −0.58 to 0.15 (0.87–0.97) |
| brown2019_WT_TKOv3 (Durocher, PACE002, PACE006) | Brown 2019 *MSB* 15:e8679, GEO GSE128210 | hTERT line | **WT** | TKOv3 | Durocher: T0→"T3" (units not stated in file); PACE002: T0→T18, T30 (3 reps); PACE006: T0→T25 | **counts** | Durocher 84736b85; PACE002 9bde1b54 (bare-CR line endings); PACE006 26f636d1; TKOv3 lib 82508210; soft 58e33295 | GEO | DepMap RPE1 | OK | T18 −3.32/0.10 (0.99); T30 −4.41/−0.17 (0.99); T25 −5.45/−0.52 (0.97); Durocher −3.13/−0.01 (0.96) |
| brown2019_p53null_TKOv2 | same | hTERT line | KO | TKOv2 | T0→T18 | **counts** | deltaP53 4a2025af; TKOv2 lib 851d5ce2 | GEO | DepMap RPE1 | OK | −2.58 / 0.04 (0.97) |
| haapaniemi2018_RPE1_Brunello (found in the Brown 2019 supplement) | Haapaniemi 2018 *Nat Med* 24:927 | hTERT line | WT and null | Brunello (77,441) | D4→D28 | **counts** | EuropePMC_PMC6686785 zip 01fad33b (`haapaniemi_NatMed_June2018_rawReadCounts.txt`) | CC BY 4.0 | DepMap RPE1 | OK | WT −0.46/−0.21 (**0.66, fail**); null −0.81/−0.12 (0.91) |
| hart2015_RPE1_TKOv1 | Hart 2015 *Cell* 163:1515 (counts re-deposited in GSE128210) | hTERT line | WT | TKOv1 lib1 (91,320) | T0→T9,T12,T15,T18 | **counts** | Hart RPE1 readcounts 7f82cc62; TKOv1 lib 62099b4e | GEO | DepMap RPE1 | OK | T18 −2.95 / 0.18 (0.95) |
| bowden2020_RPE1_DDR | Bowden 2020 *eLife* 9:e55325 | hTERT line | WT (2 clones) and KO (2 clones) | custom dual-gRNA DDR (3,370 constructs, 968 genes) | D3→D15, D19 | **counts** | EuropePMC_PMC7244323 zip a9cbc271 (`supp2.tsv` counts, `supp1.csv` library) | CC BY 4.0 | DepMap RPE1 | OK (focused library) | WT −0.75/0.20 (0.99); KO −0.80/0.22 (0.99); only 95 CEG / 39 NEG |
| mair2019_H1_MEF / _laminin | Mair 2019 *Cell Rep* 27:599, GEO GSE156284 | ESC (H1 iCas9) | WT | TKOv3 | MEF: T0→T4,T8,T12 (T12 FACS Oct4-low/high); laminin: T0→T4,T9,T12 | **counts** | readCounts matrix 335ded34; TKOv3 lib c93610e3; soft 1bc87b7d; BF tables mmc2–mmc7 | GEO; Elsevier | GSE125509 (same H1-iCas9, both substrates) | OK | MEF T8 −2.10/0.07 (0.86); laminin T12 −2.17/−0.17 (0.92) |
| yilmaz2018_hpES10 | Yilmaz 2018 *Nat Cell Biol* 20:610, GEO GSE111309 | ESC (haploid) | WT [unverified] | Wang 2015 genome-wide (182,134 guides, 18,130 genes) | Day1→Day9, 23, 30 | **counts** | rawCounts.xlsx.gz 73e658d6; MOESM3–6 | GEO; Springer | GSE103846 — **state-mismatched**: per GEO that expression is from **H9** (diploid, Rex1-EGFP), not the screened haploid h-pES10 line | OK | D23 −2.24 / −0.27 (0.98) |
| yilmaz2020_germlayers | Yilmaz 2020 *Cell Stem Cell*, PMID 32735778, GEO GSE143369 | ESC-derived lineages (DE, mesoderm, neuroectoderm) | WT [unverified] | Wang 2015 | differentiated vs matched undifferentiated | **counts** | csv.gz caa828c7 | GEO | GSE103846 (undifferentiated only) | OK, but a differentiation contrast, **not** a fitness dropout | NE vs undiff −0.08/−0.04 (0.53): expected |
| ihry2019_hPSC | Ihry 2019 *Cell Rep* 27:616 | ESC (hPSC iCas9; H1 per paper [unverified]) | WT | Novartis genome-wide, 2 pools (93,105 guides, 18,626 genes) | no-dox D0 → dox D8, D14, D18 (+ plasmid, thaw) | **counts** | mmc2 41d11c1d (counts and DESeq2 guide LFC); mmc3–7 | Elsevier | GSE125509 H1 as proxy | OK | D18 −1.35 / −0.12 (0.92) |
| tian2021_neuron_CRISPRi | Tian 2021 *Nat Neurosci* 24:1020 | iPSC-derived lineage (WTC11 neurons) | WT | hCRISPRi-v2 genome-wide | per paper [unverified] | **gene-level phenotype** (no guide counts found) | MOESM3 82785da6 (survival ±AO, CRISPRi/a, FACS screens); MOESM6 sgRNA seqs | Springer | GSE124703 Quant-seq (WTC11 neurons) | OK (gene level only) | −1.94 / −0.01 (0.79) |
| tian2019_iPSC_neuron_CRISPRi | Tian 2019 *Neuron* 104:239 | iPSC (WTC11 D10) and iPSC-derived neurons (D14/21/28) | WT | CRISPRi H1 sublibrary (~2,300 genes) | as listed | **gene-level phenotype** | mmc2 f8002842; mmc3 149ec75d (validation guide phenotypes) | Elsevier | GSE124703 Quant-seq; Dräger iTF-iPSC RNA-seq | OK (sublibrary) | iPSC −1.29/0.06 (0.85); neuron D28 −0.84/0.07 (0.95) |
| drager2022_microglia_CRISPRi | Dräger 2022 *Nat Neurosci* 25:1149 | iPSC-derived lineage (iTF-microglia) | WT | CRISPRi sublibrary (5,109 rows / 2,319 genes) | per paper | **gene-level phenotype** | MOESM4 fca88546 (also has RNA-seq of iTF-iPSC and iTF-microglia) | Springer | in-file RNA-seq (Supp Table 1) | OK, weak | −0.08 / −0.01 (**0.64, fail**) |
| shifrut2018_CD8T_GW | Shifrut 2018 *Cell* 175:1958 | primary haematopoietic (CD8+ T, 2 donors) | WT | Brunello | stimulation → CFSE sort, dividing vs non-dividing | **MAGeCK-normalised guide counts** (per sample, inside sgRNA table) + gene-level MAGeCK | mmc4 ba238488 (GW screen); mmc6 34559f99 (drug arms); mmc2/3 | Elsevier | GSE204862 / GSE118165 CD8 | OK, but **sort-based**, not a T0 dropout | dividing vs non-dividing −0.15/0.02 (0.83) |
| ting2018_CD4T | Ting 2018 *Nat Methods* 15:941 | primary haematopoietic (CD4+ T) | WT | Novartis; deposited RPM covers one chip (13,003 guides, 2,483 genes) | D0 → D6 sorted CD4+CD45+CXCR4+ | **guide RPM** | MOESM3 86861ed5 | Springer | GSE118165 CD4 / GSE204862 | OK (partial library) | −0.48 / 0.04 (0.97) |
| ting2018_CD34_HSPC | same | primary haematopoietic (CD34+ HSPC) | WT | Novartis (64,399 guides) | expansion; CD34+ vs CD34− sort | **guide RPM**, sort-based, **not fitness** | MOESM4 a645cf23 | Springer | GSE74246 HSC/MPP | OK, unsuitable as a fitness test | 0.57: expected |

## 2. Anchor screens (category = anchor)

| id | citation | library | DepMap ModelIDs | data level | files (sha256 prefix) | terms | status | sanity |
|---|---|---|---|---|---|---|---|---|
| hart2015_panel_guideLFC | Hart 2015 *Cell* (matrix in the Brown 2019 supplement) | TKOv1 | DLD-1 ACH-001061; HCT 116 ACH-000971; HeLa ACH-001086 (expression only, **no CRISPR profile**); GBM (patient line, none); RPE1 → RPE1-ss sublines | **guide LFC** (no per-line raw counts found except RPE1) | Brown zip 01fad33b | CC BY 4.0 | OK | AUROC 0.93–0.98 |
| hart2015_panel_BF | Hart 2015 Table S2 | TKOv1 (+A375 GeCKO, HCT116 shRNA) | as above | BF | mmc3 860fb20c | Elsevier | OK | 0.97–0.99 |
| hart2017_g3_BF | Hart 2017 *G3* 7:2719 | TKOv1 / Wang / Yusa / GeCKO re-scored with BAGEL | 17 lines incl. HAP1 ACH-002475, HCT116, K562 ACH-000551, HT29 ACH-000552, KBM7 (none), Jiyoye (exact name not in Model.csv; ACH-002252 "JiyoyeP-2003" is a candidate) | BF (Tables S1, S3) | G3 zip 59f03132 | CC BY 4.0 | OK. **No TKOv3 HAP1/HEK293T counts in G3** | 0.97–1.00 |
| mair2019_HAP1_TKOv3 | Mair 2019 (HAP1 arm in GSE156284) | TKOv3 | HAP1 ACH-002475 | **counts** T0→T3,T6,T10,T14,T18 (+plasmid) | matrix 335ded34 | GEO | OK | T18 −4.32/0.20 (0.995) |
| tzelepis2016_AML / _HT29_HT1080 | Tzelepis 2016 *Cell Rep* 17:1193 | Yusa v1 (90,710) | MOLM-13 ACH-000362, MV4-11 ACH-000045, OCI-AML2 ACH-000113, OCI-AML3 ACH-000336, HL-60 ACH-000002 (no CRISPR profile), HT-29 ACH-000552, HT-1080 ACH-000054 | **counts** (plasmid ref.; HT-29 timecourse D7–D25) + Kallisto FPKM | Cell Rep zip ff55d06a (`mmc6.zip`) | CC BY 4.0 | OK | 0.96–1.00 |
| sanson2018_A375_Brunello | Sanson 2018 *Nat Commun* 9:5416 | Brunello | A-375 ACH-000219 | **counts** (pDNA ref.) | Nat Commun zip ccac41c2 (MOESM4) | CC BY 4.0 | OK | −2.35/0.14 (0.99) |
| doench2016_Avana | Doench 2016 *Nat Biotechnol* 34:184 | Avana (+GeCKOv1/v2) | HT-29, A-375 | **guide LFC** vs pDNA | MOESM8 zip 50657eb7 | Springer | OK | HT29 −2.28/−0.50 (0.99) |
| wang2017_AML_CS | Wang 2017 *Cell* 168:890 | Wang 2017 genome-wide | 14 AML lines (EOL-1, HEL, MOLM-13, MonoMac1, MV4-11, NB4, OCI-AML2/3/5, P31/FUJ, PL-21, SKM-1, TF-1, THP-1; ModelIDs in `orcs_anchor_screens.csv`) | gene CS (no guide counts in supplement) | mmc3 82398904; mmc7 cfe23b9e (focused, 42 heme lines) | Elsevier | OK | 0.995–0.998 |
| ORCS anchors (45 screens) | Hart 2015, Wang 2014/15/17, Tzelepis, Doench, Steinhart 2017, Bertomeu 2017, Chan 2022 | various | `results/external_screens/orcs_anchor_screens.csv` (KBM-7, Jiyoye, HEK293-ACE2-TMPRSS2 unmapped) | gene-level | in ORCS tarball | MIT | OK | mostly 0.95–1.00. Tzelepis ORCS "MaGeCK Score" and Bertomeu "RANKS" have unknown direction (no AUROC) |

MCF10A (ACH-001357; no expression either), BJ (ACH-003283 has expression; BJ hTERT ACH-001018 has none), HEK-293 (ACH-001085) and HeLa (ACH-001086) are in `Model.csv` but have no CRISPRGeneEffect row; HEK293T is not in `Model.csv`. The CRISPR row checks ran before `essml.guard` began blocking direct CRISPR loads; expression rows were re-checked after.

## 3. Non-transformed screens available only as ORCS gene-level scores (BioGRID ORCS 2.0.18)

- Metadata: `results/external_screens/orcs_nontransformed_screens.csv`. It has 193 rows for curated non-transformed cell lines, with `fitness_like` and `baseline_fitness` flags plus library, duration and score type.
- Screen files are in `data/external_screens/biogrid_orcs/extracted/`.
- 36 are baseline fitness screens. The ones not already covered above:

| ORCS id | cell / category | library, duration, score | sanity AUROC | note |
|---|---|---|---|---|
| 799 | hTERT-RPE1 (Behan 2019, Project Score) | Yusa v1.1, 14 d, scaled BF | 0.97 | Sanger pipeline |
| 1791 / 1792 | hTERT-RPE1 WT / TP53-KO (Xu 2021) | BARBEKO base-editor KO, 15 d, fitness score | 0.95 / 0.94 | base-editing knockout |
| 1932 / 1933 (+1934/35 3D) | hTERT-RPE1 WT / TP53-KO (Drainas 2020) | GeCKOv2, 7 d, mean FC | 0.93 / 0.93 (3D: 0.72 / 0.74) | short screen |
| 1487 / 1488 | hTERT-RPE1 (Wegner 2019) | 117-gene DUB library | n/a | too small |
| 572 / 573 | NSC-CB660 / NSC-U5 fetal NSC (Toledo 2015) | GeCKO, 21–23 d, BF | 0.92 / 0.95 | |
| 1030 / 1031 | HF7450 / HF6562 fetal NSC (MacLeod 2019) | TKOv3, 14 doublings, BF | 0.90 / 0.95 | |
| 613 / 614 | HMEC (Sack 2018) | Wang 2015, 11 doublings, p-value | n/a (direction not encoded) | expression: GSE109326 EV controls |
| 2248–2271 | HMEC (Tuano 2023) | targeted ~1,950-gene KO and CRISPRi, 2D and 3D | 0.82–1.00 (KO/i) | CRISPRa arms 2274–2285 fail, as expected |
| 2086 / 2087 | MCF-10A FBXW7-del background (Baxter 2023) | Yusa v1, 2 wk, Z | 1.00 | mutant background |
| 94 | HUES62 (Shalem 2014) | GeCKO, 14 d | 0.48 | only the top 1,000 genes deposited: unusable |
| 2089 | HUES62 (van der Weegen 2021) | TKOv3, ELOF1-KO synthetic-lethal Z | 0.20 | differential score, not fitness |
| 1905 (+1900–1904 conditions) | primary T cells (Carnevale 2022) | Brunello, 17 d, ZLFC | 0.79 | expression: GSE204862 (same study) |
| 1107 / 1109 / 1110 | primary CD8+ T (Shifrut 2018) | Brunello, MAGeCK | n/a | same screens as section 1 |
| 2443 | HUDEP-2 (Xu 2021) | 773-gene UPS library, 8 d | 0.97 | HPV16 E6/E7-immortalised: flag |
| 2485 / 2486 | primary small-intestinal organoid (Lin 2023) | 1,800-gene TF library, beta | 0.98 / 0.97 | |
| 2017 | stem-cell-derived islets (Sintov 2022) | Brunello, 10 wk, PBMC exposure | 0.55 | immune-selection screen, not fitness |

ORCS does **not** contain Mair 2019, Ihry 2019, Hart 2017 TKOv3 HAP1/HEK293T, Tian 2021, Dräger 2022, Ting 2018 or Aregger 2020.

## 4. Expression sources

| id | cell types | level / IDs | files (sha256 prefix) | terms | serves | status |
|---|---|---|---|---|---|---|
| DepMap 26Q3 (local) | RPE1-ss48, ss77 (near-diploid controls), ss6, ss119, ss111, ss51 = ACH-002462..2467; HAP1; all anchor lines | log2(TPM+1), RSEM/GENCODE, protein-coding: **same pipeline as training** | `data/OmicsExpression…csv` (see DATA.md) | DepMap terms | all RPE1 screens, anchors | **OK**: all 6 RPE1 rows have expression; 5 have CRISPR (not ss111 ACH-002466). BJ (ACH-003283; not 'BJ hTERT' ACH-001018), HEK-293, HeLa, HL-60 also have expression |
| mair2019_h1_icas9_gse125509 | H1-iCas9 WT / AAVS1 / DRICH1-KO on laminin and on MEFs (17 samples, + MEF-alone) | HTSeq counts, Ensembl versioned IDs (58,042 rows incl. `__` summary lines) | RAW.tar 878b14bc; soft 8d6a25f1 | GEO | Mair; proxy for Ihry, HUES62 | OK |
| yilmaz2018_hpes10_gse103846 | hESC siRNA controls (Luciferase, Control) + KDs (16 samples). **GEO sample characteristics: source cell line H9, Rex1-EGFP (diploid). The screened line (haploid h-pES10, Yilmaz 2018) and the profiled line may differ, so this screen–expression pairing is state-mismatched (different line and ploidy): proxy only** | raw counts, Ensembl (60,675) | rawCounts 1a2680c3 | GEO | Yilmaz 2018/2020 | OK |
| tian2019_wtc11_neuron_quantseq | WTC11 iPSC-derived neurons D14, D21, D28 (2 reps each) | Quant-seq 3′ counts, Ensembl (60,204). 3′ tags, so no length normalisation; not the same as full-length RSEM | 6 GSM files (a5fb22eb …) | GEO | Tian 2021/2019 neurons | OK |
| drager2022 Supp Table 1 (inside screen file) | iTF-iPSC D0, iTF-microglia D9/D14/D15 (± LPS), Brownjohn microglia | normalised expression (units per paper [unverified]), Ensembl versioned | MOESM4 fca88546 | Springer | Dräger microglia; Tian 2019 iPSC (iTF-iPSC proxy) | OK |
| carnevale2022_primary_t_gse204862 | primary human T cells (CD8/CD4 mixes; control vs RASA2-KO; stim 1 and 5; CAR/TCR) | gene counts, symbols. **Pre-filtered to ~13.5–14.4k genes** | 3 matrices 0c586621 / c5dd4d38 / a35b0526 | GEO | Shifrut, Ting CD4, Carnevale ORCS 1905 | OK |
| calderon2019_immune_gse118165 | 25 sorted immune cell types incl. CD8/CD4 naive/memory/effector, ± stimulation (166 samples) | gene abundance, Ensembl (20,206) | c9604e7c | GEO | activated CD8/CD4 T references | OK |
| corces2016_haem_gse74246 | HSC, MPP, LMPP, CMP, GMP, MEP, mono, T, B, NK, ery (81 samples) | gene counts, symbols (25,498) | fb6bc33b | GEO | Ting CD34+ HSPC | OK |
| sack2018_hmec_gse109326 | HMEC empty-vector EV-1/EV-2 (+ MYC/CCND1/SRC/KRTAP10-6/TP53-null) | voom-normalised counts, Ensembl versioned | RAW.tar b48f98bf | GEO | ORCS HMEC screens 613/614/2248+ | OK |
| tzelepis2016 Kallisto FPKM (in screen zip) | 5 AML lines, HT1080 (ENCODE), HT29 | FPKM, symbols | Cell Rep zip ff55d06a | CC BY | anchors (DepMap expression exists anyway) | OK |
| ENCODE (H1, H9, WTC11, primary CD8/CD4 T, CD34+ CMP, HMEC, MCF 10A, NPC) | – | RSEM gene quantifications on GENCODE v24, the closest match to DepMap | see DATA.md | ENCODE open | proxies + matched DepMap-line validation | **portal still down on a 2nd attempt**; resolved via ENCODE-DCC GEO deposits (deposit `encode_geo_polyA_rnaseq`: K562, A549, SK-N-SH, T47D, PANC-1, Jurkat, BE2C, Ishikawa, HepG2, MCF-7, H1-hESC, HMEC). No WTC11/H9/primary-T/CD34 human polyA sets in that GEO collection |
| HipSci, GTEx, DICE | – | – | – | – | – | not attempted: better-matched sources were found above |

Conversion to DepMap's scale (TPM over all genes with Ensembl 116 exon-union lengths, subset to DepMap's protein-coding Entrez columns, log2(TPM+1)), the harmonisation-method comparison and applicability-domain distances: `docs/expression_harmonisation.md`.

## 5. What failed and why

| item | outcome |
|---|---|
| BioGRID ORCS `LATEST` alias | Drops at ~56 MB and refuses byte ranges (curl 18/33). **Resolved**: versioned 2.0.18 object from the public bucket; md5 equals the ETag. |
| PMC supplement files (pmc.ncbi.nlm.nih.gov `/bin/`) | JavaScript proof-of-work bot challenge; **not circumvented**. Same files taken from the Elsevier and Springer CDNs. |
| Europe PMC `supplementaryFiles` for non-OA articles (Shifrut, Tian 2021, Aregger 2020, Doench 2016) | "not open access" or empty. europepmc.org `/bin/` returns 403. |
| Europe PMC bundle for Carnevale 2022 | Timed out at 900 s (25 MB received). Gene-level data are in ORCS 1900–1905 and expression in GEO. |
| Wang 2015 *Science* supplement | 403 (bot protection). Gene CS in ORCS 149–152. |
| ENCODE | Server outage (504). |
| CRISPRbrain | JS single-page app, no static bulk download; API needs auth. Tian 2021 and Dräger are covered by journal supplements, but only at gene level. |
| Hart 2015 non-RPE1 raw counts | Not deposited anywhere found; guide LFC + BF only. |
| Hart 2017 *G3* TKOv3 HAP1/HEK293T counts | G3 has BF only. HAP1 TKOv3 counts come from Mair 2019 instead; **no HEK293T TKOv3 counts**. |
| Aregger 2020 HAP1 | Not OA; no GEO located. |
| Tian 2021 / Dräger 2022 / Tian 2019 guide-level counts | Not found in supplements or GEO (Tian GEO series hold CROP-seq/Quant-seq only). |
| Ting 2018 T-cell guide RPM | Supplement covers one library chip (~2.5k genes), not the full genome-scale library. |
| HipSci iPSC fitness screens | No deposit (per catalogue). |

## 6. Caveats for the harmonisation step

- **T0 definitions differ:**
  - Cell-population T0: Olivieri, Brown, Hart, Mair, Ihry.
  - Early timepoint as T0: Day1 (Yilmaz), D3 (Bowden), D4 (Haapaniemi).
  - Plasmid/pDNA as T0: Tzelepis, Sanson, Doench.
- **Sort-based, not dropout:** Shifrut (dividing vs non-dividing), Ting CD34 (CD34+ vs CD34−) and Yilmaz 2020 (differentiated vs undifferentiated). The last two show no CEG/NEG separation and must not be scored as fitness.
- **CRISPRi screens** (Tian 2019/2021, Dräger) are knockdown phenotypes. Evaluate them separately from KO.
- **TP53-WT RPE1 screens differ by library and lab:**
  - Haapaniemi WT (Brunello, D4→D28) separates CEG/NEG weakly (AUROC 0.66 vs 0.91 in TP53-null).
  - The Brown/Durocher TKOv3 WT screens separate strongly (0.96–0.99).
- **Mair MEF arm:** lower AUROC (0.86) and fewer guides passing the T0 ≥ 30-read filter (52,874 of 71,090) even though T0 (ESC_1) has 44.7M reads. The T0 guide distribution is skewed, not shallow; the laminin T0 (34.1M reads) keeps 69,872.
- **Bowden 2020:** small NEG set (39 genes). Its AUROC is based on few genes.
- **Olivieri arms:** 31 drug arms share T0 with the NT arms. Only the NT arms are listed as fitness screens.
- **Duplication:** Brown 2019 supplement `matrix-rpe1_p53Null.txt` (columns mislabelled 'WT') and GEO deltaP53 TKOv2 are the same samples, re-counted: same 70,697 guides, Spearman 0.99 per column, counts about 4% lower. Count it once.
