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

## Licence notes
(filled from the data-resources survey)
