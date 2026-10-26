# Research tree — CRISPR essentiality prediction (DepMap 26Q3)

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding

## Status summary (latest first)

### 2026-10-26 — M1: literature done, design v0 in review round 1
LITERATURE.md synthesises four agent surveys: prior methods, non-transformed screens, harmonisation/labels and data
resources. The decision-relevant points:
- expression dominates as a feature, and linear models (ridge/LASSO) match or beat deep models on held-out lines;
- the literature lacks lineage-mean baselines, patient-grouped CV and selective-dependency AUPRC;
- non-transformed KO screens are concentrated in RPE1 and hPSC, and iPSC-derived screens are CRISPRi-only;
- anchors (DepMap lines screened externally) make screen harmonisation testable.

The data audit found 1214 screened lines with full E+C+M omics and 48 shared-patient groups. The locked split
(1022 dev / 255 test / 15 non-cancerous) is committed. A dev-only probe was run, but review found it FLAWED (λ tuned on validation; imputation used validation lines), so its numbers are not cited; see tree item 1.3. External screen acquisition (including anchor screens) runs in parallel with design review round 1
(7 reviewers).

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
- D1: Locked the test split before any modelling. Population = screened models with expression, excluding
  Non-Cancerous models and their patients. Grouped by PatientID, stratified by lineage; 20% test. Why: 104 screened
  lines share patients, and the test set must never inform design choices.
- D2: Primary essentiality label = gene effect < −0.5 (−1 as sensitivity), not probability of dependency.
  Why: that file is not local, and it is partly circular with expression.

## Tree
- [x] 0. Setup: PROMPT.md, tree, LESSONS, DATA.md, checksums, external resource downloader
- [x] 1. Literature and dataset survey → LITERATURE.md (+ docs/lit/*)
  - [x] 1.1 Data audit (docs/data_audit.md) and parquet cache (src/essml/io.py)
  - [x] 1.2 Locked split (configs/splits.csv)
  - [-] 1.3 Dev-fold probe (scripts/probe/): ridge 0.15/0.31 vs lineage 0.09/0.18 mean per-gene r. FLAWED (R1 agent-failure
    review): λ chosen on the validation fold, NaNs mean-filled using validation lines, gene set ≠ P1. Order-of-magnitude
    calibration only; never cite as a result.
- [~] 2. DESIGN.md draft → review rounds (≤3) → frozen design
  - [~] 2.1 Round 1: external PI, statistics, prior work, data, algorithms/compute, agent failure modes, biology
- [~] 6.0 External screen and expression acquisition (docs/external_screens_inventory.md)
- [ ] 3. Data pipeline and verification (toy data, controls)
- [ ] 4. Baselines (per-gene mean, per-gene lineage mean)
- [ ] 5. Models and feature ablation
- [ ] 6. Non-transformed generalisation (assemble screens, harmonise, apply)
- [ ] 7. REPORT.md
