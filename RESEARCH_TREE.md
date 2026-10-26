# Research tree — CRISPR essentiality prediction (DepMap 26Q3)

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding

## Status summary (latest first)

### 2026-10-26 — M0: project start
Workspace holds 8 DepMap files supplied by the user (CRISPRGeneEffect, Model, ModelCondition, OmicsCNGeneWGS,
OmicsExpressionTPMLogp1HumanProteinCodingGenes, OmicsFusionFiltered, OmicsSomaticMutationsMatrixDamaging/Hotspot).
CRISPRGeneEffect covers 1361 models × 18434 genes. 17 screened models are annotated "Non-Cancerous" (RPE1 sublines,
HA1E, immortalised prostate, lymphoblastoid). Compute: 12 CPU cores, 16 GB RAM, no GPU, ~54 GB free disk.

The DepMap portal file API sits behind a Cloudflare bot check, so it cannot be scripted. Extra DepMap modalities
(methylation, proteomics, metabolomics, signatures, screen metadata) need a manual download by the user or an open
mirror. Next: literature and dataset survey (parallel agents), plus a data audit of the local files.

## Decision log
- D0 (2026-10-26): Treat the user-supplied files as the DepMap Public 26Q3 release. The user states this and the file
  dates fit, but the files carry no release tag; it is recorded as an assumption in DATA.md. Do not bypass the portal bot check.

## Tree
- [~] 0. Setup: PROMPT.md, tree, LESSONS, DATA.md, checksums
- [ ] 1. Literature and dataset survey → LITERATURE.md
- [ ] 2. DESIGN.md draft → review rounds (≤3) → frozen design
- [ ] 3. Data pipeline and verification (toy data, controls)
- [ ] 4. Baselines (per-gene mean, per-gene lineage mean)
- [ ] 5. Models and feature ablation
- [ ] 6. Non-transformed generalisation (assemble screens, harmonise, apply)
- [ ] 7. REPORT.md
