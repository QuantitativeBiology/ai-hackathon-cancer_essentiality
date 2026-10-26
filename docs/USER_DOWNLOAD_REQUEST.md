# Files that need a manual download from the DepMap portal

The DepMap portal (https://depmap.org/portal/data_page/) is behind a Cloudflare bot check, so these files cannot be
fetched by script. If you can download them from **DepMap Public 26Q3**, place them in `data/` unchanged. They enter
the pipeline automatically (checksums are added to DATA.md). File names follow the 24Q4/26Q1 conventions; the 26Q3
names may differ slightly.

## Priority 1 (screen metadata: library-stratified reporting and harmonisation)
- `CRISPRScreenMap.csv` — screen → model, library
- `ScreenSequenceMap.csv` — screen → pDNA batch, days (screen duration)
- `AchillesCommonEssentialControls.csv`, `AchillesNonessentialControls.csv` — DepMap's scaling controls

## Priority 2 (extra feature modalities for the ablation)
- `OmicsSignatures.csv` (mutational signatures, MSI, ploidy, LOH summaries)
- `OmicsCNSegmentsProfile.csv` or `OmicsAbsoluteCNSegmentsProfile.csv`
- `harmonized_Sanger_MS_2022.csv` (or the 26Q3 proteomics file)
- Methylation (CCLE RRBS gene/promoter-level file) and metabolomics (`CCLE_metabolomics_20190502.csv`), from the
  CCLE legacy section

## Priority 3
- `CRISPRGeneDependency.csv` (probability of dependency; sensitivity labels only)
- `OmicsDefaultModelProfiles.csv`, release README
