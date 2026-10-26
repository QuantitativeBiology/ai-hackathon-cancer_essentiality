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
| nusinow2020 | Nusinow et al. 2020 Cell (CCLE TMT proteomics, normalised protein quant; file dated 2020-02-06) | https://gygi.hms.harvard.edu/data/ccle/protein_quant_current_normalized.csv.gz | b72a9ff39b123c80 | no licence stated on download site; published data of Nusinow et al. Cell 2020 (doi:10.1016/j.cell.2019.12.023); research use, not redistributed | features.build_proteomics (P arm) |

## Licence notes
(filled from the data-resources survey)
