# Data resources catalogue

Checked 2026-10-26 with `curl -sIL -m 20` (HTTP status given). Note: "403" from cog.sanger.ac.uk means the guessed file name does not exist (it is a missing-key response), not that data is restricted. Raw data are never committed; DATA.md records source/version/checksum.

## 1. DepMap Public release

**Release status.** The DepMap download API (`https://depmap.org/portal/api/download/files`) returns HTTP 200 with a ~4.8 kB body, i.e. the Cloudflare challenge page, not file data; `https://depmap.org/portal/data_page/` behaves the same (4.8 kB). We did not try to bypass it. 26Q3 release notes could not be fetched; the latest announcement I could read is 26Q1 (https://forum.depmap.org/t/announcing-the-26q1-release/4606). The local files in `data/` are the ground truth for 26Q3. The file inventory below comes from the last open mirror, **figshare "DepMap 24Q4 Public"** (article 27993248, licence **CC BY 4.0**, metadata via `https://api.figshare.com/v2/articles/27993248` = 200). DepMap staff stated that later releases are portal-only (forum thread "Availability of DepMap 25Q2 on figshare"). Figshare file hosts (`ndownloader.figshare.com/files/<id>`) answered 403 to HEAD and 302 to a ranged GET: usable as a **24Q4-only** fallback (e.g. OmicsSignatures id 51065726, README id 51065795), not a 26Q3 source. Broad GCS bucket `storage.googleapis.com/depmap-external-downloads/` = 403 (no public listing). File names in 26Q3 may differ from 24Q4; verify on the portal.

**Licence.** Figshare record says CC BY 4.0 (attribution; redistribution allowed, but we do not commit). The 26Q1 portal page only shows a Broad copyright and a Terms and Conditions link; treat as CC BY 4.0 and cite release plus original-source papers.

**Files of interest (24Q4 sizes; L = already local)**

| File | What | Size | Notes / circularity |
|---|---|---|---|
| CRISPRGeneEffect.csv (L) | Chronos gene effect | 429 MB | Target; training lines only as input |
| CRISPRGeneDependency.csv | P(dependent) | 421 MB | Target-derived; labels for AUPRC, training lines only |
| CRISPRInferredCommonEssentials.csv | pan-line common essentials | 21 kB | Derived from ALL lines incl. test: circular unless recomputed on training lines |
| AchillesCommonEssentialControls.csv / AchillesNonessentialControls.csv | Hart-style control sets | 17 / 11 kB | Prior knowledge, low risk |
| AchillesHighVarianceGeneControls.csv | high-variance genes | 7 kB | Screen-derived: circular |
| CRISPRScreenMap.csv, ScreenSequenceMap.csv | screen to model/sequence IDs | 36 / 347 kB | Needed to map screens to models and libraries (Avana/KY/Humagne) |
| Avana/KY/Humagne GuideMap, RawReadcounts, LogfoldChange | guide-level data | 16 MB to 3.4 GB | Only for non-cancer harmonisation |
| CRISPRInferred{LibraryEffect,ModelEfficacy,ModelGrowthRate,GuideEfficacy}.csv, CRISPRGeneEffectUncorrected.csv | Chronos internals | up to 418 MB | Growth rate/efficacy are screen-derived covariates: circular for test lines |
| ScreenGeneEffect/Dependency.csv | screen-level target | ~500 MB | Alternative granularity |
| Model.csv (L), ModelCondition.csv (L) | metadata | 0.6 MB | SangerModelID 1220/2554, COSMICID 977, CCLEName 2122 |
| OmicsProfiles.csv, OmicsDefaultModelProfiles.csv, OmicsDefaultModelConditionProfiles.csv | profile-to-model maps | 0.1-0.25 MB | Needed to choose one profile per model |
| OmicsSignatures.csv, OmicsSignaturesProfile.csv | mutational signatures | 132 / 194 kB | Not local; needed for signature ablation |
| OmicsCNGene.csv / OmicsAbsoluteCNGene.csv; OmicsCNGeneWGS.csv (L) | gene-level CN | 1.4 GB / 239 MB | |
| OmicsCNSegmentsProfile.csv, OmicsAbsoluteCNSegmentsProfile.csv | CN segments | 41 / 34 MB | Derive arm-level CN and ploidy (no dedicated ploidy file) |
| OmicsLoH.csv | gene-level LOH | 237 MB | Optional |
| OmicsMicrosatelliteRepeatsProfile.csv | MSI repeats | 93 MB | MSI derivable; optional |
| OmicsExpressionProteinCodingGenesTPMLogp1.csv | expression | 507 MB | local copy is the "HumanProteinCodingGenes" variant |
| OmicsFusionFiltered.csv (L); OmicsSomaticMutationsMatrixDamaging/Hotspot (L) | | | Full OmicsSomaticMutations.csv (339 MB) adds other variant classes |

**Methylation, proteomics, metabolomics are not in the 24Q4 figshare listing.** The 26Q1 notes list new `harmonized_Sanger_MS_2022.csv` (Goncalves 2022 proteomics) and `harmonized_Olink_2023_best_dilution.csv` (Olink Explore HT, 161 lines), both portal-only. Methylation (CCLE RRBS TSS1kb, Ghandi 2019) and metabolomics (Li 2019, Nat Med 25:850; 225 metabolites, 928 lines; legacy name "CCLE_metabolomics_20190502") were on the portal's CCLE/legacy page; my guessed `data.broadinstitute.org/ccle_legacy_data/...` URLs returned 404, so current locations are **unverified**.

## 2. Non-DepMap profiles of the same lines

ID mapping: Model.csv has SangerModelID (SIDM...) for 1220 of 2554 models, COSMICID for 977, CCLEName for 2122.

| Resource | URL | HTTP / size | Licence | Circularity |
|---|---|---|---|---|
| CMP model list | https://cog.sanger.ac.uk/cmp/download/model_list_latest.csv.gz | 200, 145 kB (csv: 954 kB) | CMP states data are publicly available, API free for non-commercial use; no explicit redistribution licence found. Do not commit | none |
| CMP mutations | https://cog.sanger.ac.uk/cmp/download/mutations_all_20250318.zip | 200, 223 MB | as above | none |
| CMP RNA-seq | https://cog.sanger.ac.uk/cmp/download/rnaseq_all_20250117.zip | 200, 1.1 GB | as above | none |
| CMP WES CN (PureCN) | https://cog.sanger.ac.uk/cmp/download/WES_pureCN_CNV_genes_20221213.zip | 200, 776 MB | as above | none |
| CMP CRISPR (Project Score) | https://cog.sanger.ac.uk/cmp/download/essentiality_matrices.zip | 200, 241 MB | as above | **CRISPR data of overlapping lines: not allowed as test-line input; leakage** |
| CMP proteomics (Goncalves 2022, 949 lines, DIA-MS) | exact file name not resolved (guesses 403); see https://cellmodelpassports.sanger.ac.uk/downloads (JS page) | - | paper open access; verify file terms | none |
| Methylation, Iorio 2016 / GDSC (450k) | GEO GSE68379 (200); ArrayExpress E-MTAB-3610 (200) | large | public repositories | none |
| Nusinow 2020 TMT proteomics (~375 lines) | https://gygi.hms.harvard.edu/data/ccle/Table_S2_Protein_Quant_Normalized.xlsx; page https://gygi.hms.harvard.edu/publications/ccle.html | 200, 60 MB | open-access Cell paper | none; ID is CCLE name |
| Ghandi 2019 CCLE supplement (RRBS, RPPA) | portal legacy page only | unverified | DepMap terms | none |

## 3. Prior-knowledge resources

| Resource | URL (status) | Licence | Notes / circularity |
|---|---|---|---|
| CORUM 5.0 (7193 complexes; NAR 2025 53:D651) | https://mips.helmholtz-muenchen.de/corum/ ; download endpoints returned 503, WebFetch hit a certificate error; corum.uni-saarland.de unreachable | **CC BY-NC 4.0** (per OmniPath resource metadata): research use, do not redistribute | manually curated; low circularity. Likely needs browser download. Fallbacks: Reactome complexes, BioPlex, hu.MAP 3.0 (https://humap3.proteincomplexes.org, 200) |
| BioGRID (physical) | https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-ALL-LATEST.tab3.zip (200); MV-Physical file 200 | MIT | literature-curated, not NLP-mined; keep physical experimental interactions only; drop genetic interactions (many from CRISPR screens) |
| BioGRID ORCS | https://downloads.thebiogrid.org/Download/BioGRID-ORCS/Latest-Release/BIOGRID-ORCS-ALL-homo_sapiens-LATEST.screens.tar.gz (200) | MIT | **Circular**: contains CRISPR hits including DepMap/Project Score. Do not use as features |
| HuRI / HI-union (Y2H) | https://www.interactome-atlas.org/data/HuRI.tsv (http 301 to https, 200); HI-union.tsv same pattern | CC BY 4.0 (OmniPath metadata) | experimental |
| BioPlex 3.0 | https://bioplex.hms.harvard.edu/data/BioPlex_293T_Network_10K_Dec_2019.tsv (200, 12.9 MB); BioPlex_interactionList_v4a.tsv (200, 5.5 MB) | open (confirm terms) | AP-MS; expression-biased by 293T/HCT116 |
| IntAct | https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip (200, 1.4 GB) | CC BY 4.0 | curated; filter non-human and any text-mining-type evidence |
| Reactome | https://reactome.org/download/current/ReactomePathways.gmt.zip (200, 298 kB) | CC BY 4.0 | expert-curated |
| Gene Ontology | https://current.geneontology.org/annotations/goa_human.gaf.gz (200, 15 MB); go-basic.obo (200, 32 MB) | CC BY 4.0 | exclude evidence code IEA; remaining annotations are literature curation, some derived from essentiality papers (mild circularity) |
| HGNC | https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt (200, 17 MB); json variant (200, 33 MB); gene groups via https://www.genenames.org/download/ (200). `genes_in_groups.txt` path 404 | CC0 | symbol harmonisation, families |
| MSigDB Hallmark 2024.1.Hs | https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/h.all.v2024.1.Hs.symbols.gmt (200, 49 kB) | CC BY 4.0 for Hallmarks; some other collections (e.g. KEGG) carry restrictions | Hallmarks only |
| Ensembl Compara paralogs | https://ftp.ensembl.org/pub/current/tsv/ensembl-compara/homologies/homo_sapiens/ (200). BioMart martservice requires POST (GET 405) | Ensembl: open | paralog type, % identity; preferred over BioMart for reproducibility |
| De Kegel and Ryan paralog pairs (Cell Syst 2021) | verified direct URL not found (supplement / Ryan lab) | journal CC BY (verify) | **Circular**: SL labels inferred from DepMap across all lines. Use only pair features (identity, shared PPI), or recompute on training lines |
| Hart CEGv2 / NEGv1 | shipped in AchillesCommonEssentialControls / NonessentialControls; Hart 2017 G3 supplement | open | derived from older cell-line CRISPR screens, independent of 26Q3 test lines; fine as prior and for label definition |
| OmniPath | https://omnipathdb.org/interactions?genesymbols=1&datasets=omnipath (200); resources https://omnipathdb.org/resources | per resource | **Flag**: aggregates text-mined sources (REACH, Sparser, RLIMS-P via ProtMapper) and HPRD (own licence) in extra datasets. Use only core curated `omnipath` dataset after checking sources, otherwise avoid |
| Co-essentiality networks (FIREWORKS, Wainberg 2021, DepMap-derived GIN) | - | - | **Circular** unless rebuilt from training lines only; do not use |
| STRING and other text-mining networks | - | - | excluded by project rules |

## 4. Recommended downloads

**P1 (core design)**
- Ensembl Compara paralogs: https://ftp.ensembl.org/pub/current/tsv/ensembl-compara/homologies/homo_sapiens/
- HGNC: https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt
- BioGRID: https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-ALL-LATEST.tab3.zip
- HuRI: https://www.interactome-atlas.org/data/HuRI.tsv
- BioPlex 3.0: https://bioplex.hms.harvard.edu/data/BioPlex_293T_Network_10K_Dec_2019.tsv
- Reactome: https://reactome.org/download/current/ReactomePathways.gmt.zip
- GO: https://current.geneontology.org/annotations/goa_human.gaf.gz and https://current.geneontology.org/ontology/go-basic.obo
- MSigDB Hallmarks: URL in section 3
- CORUM (browser download, CC BY-NC)
- Manual DepMap files (below)

**P2 (ablation)**
- Sanger CMP model list and mutations; Goncalves 2022 proteomics (DepMap harmonized file or CMP); Nusinow 2020 xlsx; GDSC methylation (GEO GSE68379); IntAct; hu.MAP 3.0.

**P3**
- Sanger RNA-seq and PureCN CN (replication, non-cancer profiling); OmicsLoH; MSI repeats; Humagne/KY raw counts for library harmonisation.

**Manual download from the DepMap portal (https://depmap.org/portal/data_page/, release 26Q3)** - names as in 24Q4, confirm against the 26Q3 list:
1. OmicsSignatures.csv, OmicsSignaturesProfile.csv
2. OmicsProfiles.csv, OmicsDefaultModelProfiles.csv, OmicsDefaultModelConditionProfiles.csv
3. OmicsAbsoluteCNSegmentsProfile.csv or OmicsCNSegmentsProfile.csv (arm-level CN, ploidy)
4. CRISPRScreenMap.csv, ScreenSequenceMap.csv, CRISPRGeneDependency.csv, AchillesCommonEssentialControls.csv, AchillesNonessentialControls.csv; CRISPRInferredCommonEssentials.csv (only with training-only restriction)
5. harmonized_Sanger_MS_2022.csv (proteomics); harmonized_Olink_2023_best_dilution.csv (optional)
6. CCLE RRBS TSS1kb methylation (Ghandi 2019) and CCLE_metabolomics_20190502.csv (Li 2019) from the legacy CCLE page
7. OmicsSomaticMutations.csv (full), OmicsLoH.csv, OmicsMicrosatelliteRepeatsProfile.csv (optional)
8. README.txt for the release
9. AvanaGuideMap.csv, AvanaRawReadcounts.csv (and KY/Humagne equivalents) only if guide-level harmonisation is needed

## Licence concerns
- DepMap: CC BY 4.0 per figshare; portal terms not explicit. Cite release and source papers.
- CORUM is CC BY-NC 4.0: research use only, no redistribution.
- CMP downloads: no explicit open licence found; do not redistribute.
- OmniPath includes text-mined and HPRD-licensed sources.
- KEGG-derived MSigDB sets restricted; use Hallmarks, Reactome, GO.
