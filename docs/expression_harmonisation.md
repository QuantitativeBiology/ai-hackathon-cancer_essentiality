# Expression harmonisation of the non-transformed profiles (DESIGN v3 §7.3)

Status: draft; method **(ii) quantile** adopted by PI decision D9 (2026-10-26) as an explicit override of the rule's outcome (none). The PI freezes the method. No screen data were loaded and no prediction model
was applied at any point.

**Reproduce:**

```
python scripts/fetch_external_screens.py --only encode_geo_polyA_rnaseq ensembl116_gtf
python scripts/harmonise_expression.py            # ~5.5 min, < 2 GB RAM
python -m pytest tests/test_expr_harmonise.py
```

**Inputs:**
- `configs/expression_harmonisation.json`: sources, conversion parameters, the selection rule (written before the
  comparison was run) and `chosen_method`.
- `configs/external_expr_profiles.csv`: the 117 samples used, with state, cell type, matched screens and DepMap match.

**Outputs:**
- `cache/external_expr/`:
  - `harmonised.parquet`: the chosen method (rows = profile id, columns = Entrez).
  - `harmonised_{none,quantile,quantile_location,rank}.parquet`, `log2tpm_unharmonised.parquet`.
  - `location_offsets.parquet`, `ensembl116_genes.parquet`.
- `results/expression_harmonisation/`:
  - `profiles.csv`
  - `method_comparison.csv`, `method_summary.csv`, `selection.json`
  - `count_path_check.csv`, `length_bias_check.csv`, `source_qc.csv`, `ad_reference.json`

## 1. Conversion to DepMap's scale
DepMap gives log2(TPM+1) on 19,218 protein-coding Entrez columns. Its TPM is computed over all genes and then
subset, so per-line protein-coding sums are about 0.94M (dev median 938,204). We follow the same steps for every source.

1. **Parse** the deposited table. Ensembl versions are stripped and duplicate IDs summed. Symbols are mapped to
   Ensembl (unique Ensembl 116 `gene_name`, then HGNC approved, then a unique previous/alias symbol).
2. **TPM over all genes:**
   - **Counts:** count / exon-union length (Ensembl 116 GTF, all transcripts' exons merged), normalised to 1e6 over
     every gene that has a length. Between 99.2% and 100% of counts fall in genes with a length.
   - **FPKM:** FPKM / sum × 1e6.
   - **ENCODE RSEM TPM:** ERCC spike-ins dropped, then renormalised (spike-ins take 0.1–9% of TPM).
   - **Quant-seq 3′ tags:** CPM with no length correction. Flagged `tpm_compatible = False`.
3. **Subset** to DepMap's Entrez columns (Ensembl → Entrez via HGNC, with the DepMap symbol as fallback; several
   Ensembl genes mapping to one Entrez are summed), then take log2(TPM+1). A gene missing from a full table stays NaN.
   Feature code already fills absent expression with the training mean. Between 96.7% and 100% of the 19,218
   columns are observed.
4. **Restricted-universe tables:** Carnevale is pre-filtered to about 14k genes and Calderon covers 20,206 genes. Neither
   can give TPM over all genes. As the closest compliant option, their protein-coding subset is rescaled to the
   DepMap dev median sum, and absent genes are set to 0 because they were filtered for low expression. Both are flagged
   `universe = restricted`.

**Effective length.** Ensembl 116 exon-union length, with no fragment-length correction (`counts_fragment_length = 0`).

**Check on ENCODE** (`count_path_check.csv`, 23 files). Our TPM from the same files' RSEM `expected_count` agrees
with ENCODE's own RSEM TPM at mean r = 0.963, MAE = 0.46 log2 units and mean shift −0.10. Subtracting a 200-bp
fragment length makes it worse (r = 0.958, MAE = 0.57, shift −0.36). So 0 was kept. This fixes the count path's error
level: the count-based sources sit about 0.46 log2 units of MAE away from an RSEM quantification of the same reads.

**Unit check for Dräger Supp. Table 1** (`length_bias_check.csv`). Every column sums to exactly 1e6 over 60,232
genes. We regressed log2(profile / ENCODE H1 TPM) on log2 gene length:

| profile | slope |
|---|---|
| ENCODE H1 replicate (TPM) | 0.03 |
| Mair H1 count path | −0.18 |
| Quant-seq neurons | +0.52 |
| ENCODE cancer lines (TPM) | −0.03 to +0.37 |
| Dräger iTF-iPSC | +0.30 |

Dräger's slope sits within the range of known TPM tables and far from the ≈ +1 expected for CPM. Treated as TPM
**[unverified: the paper's units were not checked]**.

## 2. Validation set: matched external RNA-seq of DepMap lines
encodeproject.org still timed out on a second attempt. The ENCODE DCC's own GEO deposits of the uniform pipeline's
GRCh38 gene quantifications ("polyA mRNA RNA-seq from …") were used instead, together with the Kallisto FPKM table
in the Tzelepis 2016 supplement. One validation profile = one line × source, averaged over replicates by mean TPM.

| stratum | lines (DepMap ModelID) | source |
|---|---|---|
| **exact dev (selection population, n = 12)** | K562, A549, SK-N-SH (untreated), T47D (DMSO vehicle), PANC-1, Jurkat E6-1 | ENCODE (2 replicates each) |
| | MOLM-13, MV4-11, OCI-AML2, OCI-AML3, HT-1080, HT-29 | Tzelepis 2016 (Kallisto; the HT-1080 and HT-29 reads come from ENCODE and GSE41586) |
| near-match dev (reported only) | BE2C → SK-N-BE(2); ENCODE Ishikawa → DepMap "Ishikawa (Heraklio) 02 ER-" | ENCODE |
| test lines (reported separately, DESIGN §7.2) | HepG2, MCF-7 | ENCODE |
| no split role | HL-60 (has DepMap expression, never screened) | Tzelepis |

## 3. Candidate comparison (`method_summary.csv`, `method_comparison.csv`)
**How each metric is computed:**
- Nearest-neighbour (NN) pool: all 1,022 DepMap dev lines. For a test-line profile the pool is the dev lines plus
  that test line.
- Correlation: Pearson on the top 5,000 variable genes (variance across dev lines).
- Profile r and MAE: computed against the line's own DepMap profile over all observed genes. Centred r subtracts the
  DepMap dev mean first.
- Leakage control for the fitted candidates:
  - The quantile reference excludes the matched line's own DepMap profile.
  - (iii) offsets are fitted on exact-dev states only, leave-one-matched-line-out (LOO).
  - Leave-one-source-out (LOSO) is a sensitivity analysis.
  - The fitting code asserts that the held-out line never enters its own fit (also unit-tested).

Exact dev matches, n = 12, "reported" path (ENCODE RSEM TPM, Tzelepis FPKM):

| candidate | recall@1 | recall@5 | recall@10 | recall@1 (centred) | profile r | centred r | MAE |
|---|---|---|---|---|---|---|---|
| (i) none | 0.833 | 0.917 | 1.00 | 0.833 | 0.951 | 0.750 | 0.537 |
| (ii) quantile | 0.833 | 0.917 | 1.00 | 0.833 | 0.951 | 0.760 | **0.506** |
| (iii) quantile + gene-wise location, LOO-line | 0.833 | 0.917 | 1.00 | 0.833 | 0.956 | 0.769 | 0.513 |
| (iii) same, LOSO | 0.833 | 0.917 | 1.00 | 0.833 | 0.944 | 0.725 | 0.603 |
| (iv) rank (both sides ranked) | 0.833 | 0.917 | 1.00 | 0.833 | 0.944 | 0.724 | (rank units) |
| (v) Celligner-style | deferred (DESIGN v3) | | | | | | |

**Count path.** The same comparison with the ENCODE profiles rebuilt from `expected_count` (as for most
non-transformed sources) gives:

| candidate | MAE | profile r |
|---|---|---|
| none | 0.625 | 0.934 |
| quantile | 0.595 | 0.934 |
| (iii) LOO | 0.562 | 0.949 |
| (iii) LOSO | 0.689 | 0.930 |

Recall is identical across candidates.

**Other strata.** Near-matches, test lines and HL-60 rank first under every candidate.

**The two misses are the same under every candidate:**
- ENCODE SK-N-SH: true rank 6–8. Its nearest line is **SH-SY5Y**, a subclone of SK-N-SH.
- Tzelepis MV4-11: true rank 2. Its nearest line is OCI-AML3.

**Pre-registered rule outcome: (i) none. Adopted: (ii) quantile, by PI override (RESEARCH_TREE D9).**
- All candidates tie on the primary (recall@1).
- Profile r differs by less than 0.01.
- The rule then picks the simplest candidate (`selection.json`).
- (iii) is not excluded by the LOSO guard, because its recall does not drop, but under LOSO it is *worse* than none
  on every agreement metric. Gene-wise offsets learnt from one lab do not carry over to another, and every
  non-transformed source is a new lab.

**Why the PI overrode the rule** (D9; `pi_override` in the config, rule and outcome left unchanged; the script
re-checks that the rule still returns `none` and refuses to run if it does not):
- The primary and secondary criteria are invariant to a global log shift and are saturated. The scale-sensitive MAE
  favours quantile normalisation by about 6%.
- Quantile normalisation specifically repairs library-composition shifts that "none" leaves in place. In those shifts,
  non-coding structural RNAs (RN7SL, RN7SK, RMRP, rRNA, MT pseudogenes) absorb TPM:
  - Tzelepis HT-1080: protein-coding TPM sum 0.41M; MAE 0.745 → 0.523 and centred r 0.53 → 0.65 under quantile.
  - Yilmaz hESC: protein-coding sum 0.71M.
  - ENCODE H1 replicate 2: 0.77M.
- DESIGN v3 specifies no tie-break. Quantile was chosen on MAE (an existing metric in the comparison table) and for
  the library-composition repair, before any prediction model touched external data. `selection.json` records both
  `rule_outcome` (none) and `chosen` (quantile). `cache/external_expr/harmonised.parquet` = quantile;
  `harmonised_none.parquet` keeps the rule's outcome.

## 4. Applicability domain (`profiles.csv`)
**AD distance:** mean (1 − r) to the 10 nearest DepMap dev lines on the top 5,000 variable genes, in the quantile-
normalised space (pool = DepMap dev lines as supplied). Nearest lineage = majority OncotreeLineage among those 10.

**Reference** (`ad_reference.json`): each dev line against the other dev lines.

| quantile | AD distance |
|---|---|
| median | 0.212 |
| 95% | 0.317 |
| 99% | 0.359 |
| maximum | 0.409 |

ENCODE matched cancer-line profiles sit at 0.20–0.27, except HepG2 (0.38) and MCF-7 (0.28); Tzelepis profiles sit at 0.14–0.25.

State-level profiles (mean TPM over replicates; per-sample rows are also in the table):

| profile (state) | screens served | TPM-comp. | AD | dev percentile | nearest lineage (n/10) | nearest line |
|---|---|---|---|---|---|---|
| DepMap RPE1 sublines (native) | Olivieri, Brown, Hart, Haapaniemi, Bowden | yes | 0.183 | 30 | CNS/Brain (3) | Mel285 |
| ENCODE HMEC (1 sample) | Sack, Tuano HMEC | yes | 0.210 | 48 | Head and Neck (9) | HSC-3 |
| Sack HMEC EV | Sack, Tuano HMEC | yes | 0.231 | 64 | Head and Neck (9) | HSC-3 |
| Carnevale T cells, input | Carnevale; Shifrut; Ting CD4 | yes (restricted) | 0.328 | 96 | Lymphoid (9) | KE-37 |
| Carnevale T cells, ctrl stim 1 | same | yes (restricted) | 0.311 | 95 | Lymphoid (9) | MJ |
| Calderon CD8+ T, unstim / stim | Shifrut | yes (restricted) | 0.399 / 0.353 | 99.8 / 98.5 | Lymphoid (9 / 10) | SMZ-1 / MEC-1 |
| Calderon effector CD4+ T, unstim / stim | Ting CD4 | yes (restricted) | 0.496 / 0.358 | 100 / 98.7 | Lymphoid (9 / 10) | CI / MJ |
| Calderon naive CD4+ Teff, unstim / stim | Ting CD4 | yes (restricted) | 0.418 / 0.369 | 100 / 99.5 | Lymphoid (9 / 10) | KE-37 / MEC-1 |
| Mair H1-iCas9, laminin | Mair laminin; Ihry (proxy) | yes | 0.387 | 99.6 | Testis (3) | NTERA-2 cl.D1 |
| Mair H1-iCas9, MEF | Mair MEF | yes | 0.383 | 99.6 | Testis (3) | SUSA |
| ENCODE H1-hESC | Mair / Ihry (proxy) | yes | 0.351 | 98.4 | Testis (3) | SUSA |
| Yilmaz hESC control siRNA (H9 per GEO; state-mismatched) | Yilmaz 2018 / 2020 (proxy) | yes | 0.433 | 100 | Testis (3) | SUSA |
| Dräger iTF-iPSC D0 | Tian 2019 iPSC (proxy) | yes | 0.365 | 99.5 | Testis (3) | NTERA-2 cl.D1 |
| Tian WTC11 neurons D14 / D21 / D28 | Tian 2019 / 2021 | **no** (3′) | 0.391 / 0.403 / 0.407 | 99.7–99.9 | PNS (10) | CHP-134 |
| Dräger iTF-microglia D9 / D14 PBS / D15 | Dräger 2022 | yes | 0.436 / 0.406 / 0.445 | 99.9–100 | CNS/Brain (2–3) | SF539 |
| Corces HSC / MPP / LMPP / CMP | Ting CD34 (sort, not fitness) | yes | 0.416 / 0.384 / 0.382 / 0.352 | 98.5–100 | Myeloid (8–10) | MUTZ-8 |

**What the AD table shows:**
- **Inside the domain:** RPE1 and HMEC.
- **At or beyond the dev maximum (0.409):** every stem-cell, neural, immune and HSPC profile lies at the 94th–100th
  dev percentile. Their transfer results must be read as extrapolation.
- **Nearest neighbours make biological sense:**
  - pluripotent cells → embryonal-carcinoma lines (NTERA-2, SUSA), which supports the planned EC leave-out arm;
  - neurons → neuroblastoma;
  - T cells → T/B lymphoid lines;
  - HSPC → AML;
  - HMEC → squamous carcinoma.
- **Microglia** have no lineage majority (at most 3 of 10): they are the furthest from any DepMap lineage.
- **Stimulated T cells** are consistently closer to DepMap than resting ones (proliferating lines).
- AD distances under the rule's outcome (none) differ from these by at most 0.02 for state profiles (up to 0.06
  for single Calderon samples).

## 5. Problems and caveats
- **The criterion does not discriminate.** With 12 matched lines, recall@1 is 10/12 for every candidate. The
  comparison cannot separate the candidates on the pre-registered primary. Only MAE does (see the PI flag).
- **Count path vs RSEM:** count-based profiles (most non-transformed sources) carry about 0.46 log2 units of MAE
  relative to RSEM on the same reads (HTSeq union counting vs RSEM multi-mapping and isoform handling). This is not
  correctable with the matched set, which is TPM/FPKM.
- **Restricted universes** (Carnevale, Calderon): protein-coding renormalisation is an approximation. Calderon's
  "gene abundance" is non-integer normalised counts, used as counts.
- **3′ Quant-seq** (Tian neurons): not TPM-compatible. These profiles are harmonised and given an AD distance, but are
  flagged. No full-length WTC11 neuron RNA-seq was found.
- **Proxies:**
  - GSE103846 is **H9 Rex1-EGFP** (diploid) per GEO, not the haploid h-pES10 line Yilmaz screened.
  - Ihry hPSC is matched only via H1 profiles.
  - Tian 2019 iPSC is matched via Dräger's iTF-iPSC (same WTC11 background, a different engineered line).
  - Primary T and HSPC references come from different donors and labs than the screens.
- **Mair series-level file** `GSE125509_Moffat_11_…htseqReadCounts.txt.gz` is not tied to a GSM and was not used.
- **Disclosure:** while listing the sheets of the Dräger workbook, the first five rows of its CRISPRi phenotype sheet
  (Supp. Table 3) were printed to the console. They were not stored or used. The parser reads only
  `Supplementary Table 1`.
