# DATA — sources, versions, checksums, licences

Raw data live under `data/` (git-ignored) and are never committed. Derived caches live under `cache/` (git-ignored).
Checksums for local DepMap files: `configs/depmap_local_sha256.txt` and `configs/depmap_local_md5.txt`.

## Assumptions (explicit)
- A1: The DepMap files in `data/` were supplied by the user as **DepMap Public 26Q3**. The CSVs carry no release tag,
  and the portal file index (which lists release and md5) cannot be fetched by script because of a Cloudflare bot
  check. Release identity rests on the user's statement plus file dates (2026-10-25/26).

## DepMap Public 26Q3 (local)
Source: https://depmap.org/portal/data_page/ (downloaded manually by the user). Licence: see the licence section below.

| File | Rows × cols (incl. header/meta) | sha256 (first 16) |
|---|---|---|
| CRISPRGeneEffect.csv | 1362 × 18435 | 1e0c492944a7297c |
| Model.csv | 2555 × 48 | 65eb1895800e71a9 |
| ModelCondition.csv | 3386 × 17 | a1edda78cd60718e |
| OmicsCNGeneWGS.csv | 1847 × 19345 | 31cace66f202adaa |
| OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv | 2315 × 19224 | 252d1fedc7014953 |
| OmicsFusionFiltered.csv | 102820 × 15 | 985ac6383a95c39f |
| OmicsSomaticMutationsMatrixDamaging.csv | 3808 × 19738 | 736b7a82b4d1f234 |
| OmicsSomaticMutationsMatrixHotspot.csv | 3808 × 171 | 0cbdacdce2355d2d |

## External datasets
(added as they are acquired)

Downloaded by `scripts/download_external.py` from `configs/external_sources.json`; full sha256, size and
retrieval date in `configs/external_manifest.json`. Files under `data/raw_external/` and `data/external_omics/`
(git-ignored).

| Resource | Version | Source | sha256 (first 16) | Licence | Used by |
|---|---|---|---|---|---|
| ensembl_paralogs | Ensembl 116 | https://ftp.ensembl.org/pub/release-116/tsv/ensembl-compara/homologies/homo_sapiens/Compara.116.protein_default.homologies.tsv.gz | fc101ab6e9096609 | Ensembl: no restrictions (open) | priors.py (paralogs) |
| hgnc | download date | https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt | 12363916be36f8f6 | CC0 | ID mapping (priors, features, genesets) |
| biogrid | BioGRID 5.0.262 (BIOGRID-ALL-5.0.262.tab3.txt) | https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-ALL-LATEST.tab3.zip | 7221dad43bd90ec2 | MIT | priors.py (physical partners; ORCS not used) |
| huri | HuRI 2020 | https://interactome-atlas.org/data/HuRI.tsv | 4cdf94366ff6e23e | CC BY 4.0 | priors.py (partners) |
| bioplex_293t | BioPlex 3.0 293T | https://bioplex.hms.harvard.edu/data/BioPlex_293T_Network_10K_Dec_2019.tsv | 5a98def087fc77c7 | academic use (see site) | priors.py (partners) |
| reactome_gmt | download date | https://reactome.org/download/current/ReactomePathways.gmt.zip | 8c1dbc8578431da5 | CC BY 4.0 | not yet used |
| msigdb_hallmark | MSigDB 2024.1.Hs | https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/h.all.v2024.1.Hs.symbols.gmt | ee2463540042078b | CC BY 4.0 (Hallmarks) | features.py (Hallmark scores) |
| hart_ceg2 | CEGv2 (Hart 2017) | https://raw.githubusercontent.com/hart-lab/bagel/master/CEGv2.txt | 1b1df46f829e2c47 | MIT (bagel repo) | genesets / harmonisation |
| hart_neg | NEGv1 (Hart 2017) | https://raw.githubusercontent.com/hart-lab/bagel/master/NEGv1.txt | 18a6b37a8f182a1b | MIT (bagel repo) | genesets (non-essential class) |
| deepdep_genes_1298 | Prep4DeepDEP @ f713c84 (2021-08-24), inst/extdata/default_dep_genes_1298.RData (gene list only) | https://raw.githubusercontent.com/chenlabgccri/Prep4DeepDEP/master/inst/extdata/default_dep_genes_1298.RData | 38e185e3bc166b68 | GPL-3.0 (repo); gene names only | G3 published-comparable block (scripts/g3_published_genesets.py) |
| leap_genes_1539 | owkin/leap @ 7c709df (2026-07-09), data/list_of_perturbations_task_1.csv | https://raw.githubusercontent.com/owkin/leap/main/data/list_of_perturbations_task_1.csv | c349c8a89caac3fa | repo licence file (not SPDX-detected); gene names only | G3 published-comparable block |
| nusinow2020 | Nusinow et al. 2020 Cell (CCLE TMT proteomics, normalised protein quant; file dated 2020-02-06) | https://gygi.hms.harvard.edu/data/ccle/protein_quant_current_normalized.csv.gz | b72a9ff39b123c80 | no licence stated on download site; published data of Nusinow et al. Cell 2020 (doi:10.1016/j.cell.2019.12.023); research use, not redistributed | features.build_proteomics (P arm) |

## Licence notes
(filled from the data-resources survey)

## External CRISPR screens (non-transformed test set + anchors) and matched expression
Downloaded by `scripts/fetch_external_screens.py` from `configs/external_screens.json` into `data/external_screens/` and `data/external_expr/` (git-ignored). sha256/bytes/URL per file: `configs/external_screens_manifest.json` (82 files, retrieved 2026-10-26). Sources: GEO, Mendeley (CC BY 4.0), Europe PMC OA bundles (CC BY 4.0), Elsevier/Springer supplementary files (publisher terms; do not redistribute), BioGRID ORCS 2.0.18 (MIT; test-set use only, never model input). Per-screen source, data level, licence and failures: `docs/external_screens_inventory.md`.

## External expression for harmonisation (DESIGN v3 §7.3)
Downloaded by `scripts/fetch_external_screens.py --only encode_geo_polyA_rnaseq ensembl116_gtf` (deposits in `configs/external_screens.json`); full sha256, bytes and URL per file in `configs/external_screens_manifest.json`. Retrieved 2026-10-26. Used by `src/essml/expr_harmonise.py` / `scripts/harmonise_expression.py`. Expression sources that were already on disk (GEO deposits, supplementary tables) are listed in the section above and in `docs/external_screens_inventory.md` §4.

| Resource | Version | Source | sha256 (first 16) | Licence | Used by |
|---|---|---|---|---|---|
| ensembl116_gtf | Ensembl 116, GRCh38 primary-assembly GTF | https://ftp.ensembl.org/pub/release-116/gtf/homo_sapiens/Homo_sapiens.GRCh38.116.gtf.gz | ed992f0eac7197d9 | Ensembl: no restrictions (open) | gene exon-union lengths (counts -> TPM), gene_name/biotype for symbol mapping |
| encode_geo_polyA_rnaseq | ENCODE uniform RNA-seq pipeline (STAR + RSEM, GENCODE v24, GRCh38) gene quantifications, ENCODE-DCC deposits in GEO (2016-2017) | GEO sample supplements, e.g. https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2308nnn/GSM2308418/suppl/GSM2308418_ENCFF768TKT_gene_quantifications_GRCh38.tsv.gz (encodeproject.org timed out on 2026-10-26) | per file below | ENCODE: open, no restrictions (cite the ENCODE Consortium); GEO: no NCBI restrictions | matched validation of DepMap lines; H1-hESC and HMEC proxies |

ENCODE gene-quantification files (GSM_ENCFF; line; GEO series; sha256 first 16):

| file | biosample (series, ENCODE experiment) | sha256 (first 16) |
|---|---|---|
| GSM2308418_ENCFF768TKT, GSM2308419_ENCFF172GIN | K562 (GSE86660, ENCSR000CPH) | e52cfd3e2a6fcfa4, 733082dcd05ffa9a |
| GSM2308412_ENCFF369ZNM, GSM2308413_ENCFF627QMV | A549 (GSE86657, ENCSR000CON) | 20ecf02bce4307b2, 65bcb51b866d91c9 |
| GSM2453357_ENCFF067ZMG, GSM2453358_ENCFF390KQP | SK-N-SH untreated (GSE93449, ENCSR000BYL) | 55ac915043841939, 1d49ea433cbeb675 |
| GSM2453343_ENCFF601IKR, GSM2453344_ENCFF172ATN | T47D, DMSO 4 h (GSE93442, ENCSR000BYE) | 9385709e7b91786c, 1a3ba5ada8ccc5b8 |
| GSM2453359_ENCFF248YCR, GSM2453360_ENCFF890DEQ | PANC-1 (GSE93450, ENCSR000BYM) | 8bc3764e76dffee4, 9a1dbeb0612cd3f2 |
| GSM2453329_ENCFF558JTV, GSM2453330_ENCFF489SJY | Jurkat E6-1 (GSE93435, ENCSR000BXX) | 352520f7f8b7b35e, bf30cb6ae0a00fce |
| GSM2453355_ENCFF238HUZ, GSM2453356_ENCFF396XRF | BE2C (GSE93448, ENCSR000BYK) | 6e9179477cfafe4a, f360d7c3eaf413a9 |
| GSM2453335_ENCFF394VCU, GSM2453336_ENCFF681NAQ | Ishikawa, DMSO 4 h (GSE93438, ENCSR000BYA) | 3134366cfd7fe0b0, 72529254312a104b |
| GSM2400160_ENCFF226LQB, GSM2400161_ENCFF458OHS | HepG2 (GSE90230, ENCSR000CPF) — test-role line, expression only | c5b9ca510fc77027, 1947ff5003e1d5db |
| GSM2400174_ENCFF009GDJ, GSM2400175_ENCFF885LEQ | MCF-7 (GSE90237, ENCSR000CPT) — test-role line, expression only | 903c1a10312b5be4, 0a4b757f75a22b02 |
| GSM2400152_ENCFF653XHG, GSM2400153_ENCFF915AUQ | H1-hESC (GSE90225, ENCSR000COU) | f02c5ea543188149, bfc6f0730dd5b47e |
| GSM2072420_ENCFF380GBC | mammary epithelial cell / HMEC (GSE78586, ENCSR000COX) | 64b8a2ee3b52a621 |
| GSE*_family.soft.gz (12 series) | sample metadata | see manifest |

Also read, from files already on disk: Tzelepis 2016 `mmc6.zip/RNA-seq_5AMLs-HT1080-HT29_FPKM_Kallisto.txt`, read inside `EuropePMC_PMC5081405_supplementaryFiles.zip` (sha256 ff55d06a…, CC BY 4.0) without extracting it; only the expression table, never the guide-count files. Dräger 2022 MOESM4: only the sheet `Supplementary Table 1` (expression).
