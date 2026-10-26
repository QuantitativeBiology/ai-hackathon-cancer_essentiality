# Research tree — CRISPR essentiality prediction (DepMap 26Q3)

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding

## Status summary (latest first)

### 2026-10-26 — M2: design frozen (DESIGN v3)
Three review rounds ran: round 1 with 7 specialist reviewers, round 2 with the same 7 re-checking, and round 3 with
one fresh reviewer focused on consistency and cuts. Each round moved the design toward something decisive and
enforceable. The main changes:
- P1 is computed per fold, with skill over the lineage mean as the headline;
- NaN-aware fits that respect DepMap's library-structured missingness;
- a transfer model T that uses only features available externally;
- non-transformed DepMap lines tiered by immortalisation mechanism;
- anchor-based screen harmonisation that excludes non-transformed lines and p53/DDR genes;
- pre-registered transfers/fails/inconclusive rules per biological background;
- code-level safeguards: test guard bound to committed finalists, frozen evaluation hashes, a run registry and gate
  evidence files.

In progress: G0 infrastructure and baselines (G1), feature/prior-knowledge matrices with an open proteomics arm,
and external screen acquisition (non-transformed + anchor screens). Next: verify G0, freeze the evaluation code
(tag `freeze-eval`), then the model ladder M1 → M2 → M3/M4 → M5 with ablations on dev CV.

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
- D3 (M2): P1 is computed per fold, then averaged, not from pooled out-of-fold predictions. Why: pooling makes the
  gene-mean baseline's r spuriously negative (statistics review); prior-work's pooled variant is kept only as a
  published-comparable number.
- D4 (M2): Screen-specific dependency = dependency fraction ≤ 0.5 in cancer dev (headline), ≤ 5% (secondary).
  Why: compromise between the statistics and biology reviews; it removes core-essential inflation.
- D5 (M2): Kept M4 (neural, shared across genes) despite the round-3 suggestion to cut it. Why: the goal is the best
  method; it is gated by a cheap pre-check and time-boxed.
- D6 (M2): Set-A tiers come from metadata, never from those lines' CRISPR outcomes. Why: tiering on outcomes
  conditions the evaluation on the target.

## Tree
- [x] 0. Setup: PROMPT.md, tree, LESSONS, DATA.md, checksums, external resource downloader
- [x] 1. Literature and dataset survey → LITERATURE.md (+ docs/lit/*)
  - [x] 1.1 Data audit (docs/data_audit.md) and parquet cache (src/essml/io.py)
  - [x] 1.2 Locked split (configs/splits.csv)
  - [-] 1.3 Dev-fold probe (scripts/probe/): FLAWED (λ tuned on validation; imputation used validation lines).
    Numbers deliberately not recorded; superseded by G1/G2 runs.
- [x] 2. DESIGN.md draft → review rounds (≤3) → frozen design (v3)
  - [x] 2.1 Round 1: external PI, statistics, prior work, data, algorithms/compute, agent failure modes, biology
  - [x] 2.2 Round 2: same reviewers re-check v1 (all minor/accept)
  - [x] 2.3 Round 3: fresh consistency/over-engineering review → v3 frozen
- [~] G0 evaluation infrastructure (metrics, guard, registry, tests) + G1 baselines
- [~] Features and priors (src/essml/features.py, priors.py; Nusinow proteomics)
- [x] 6.0 External screen and expression acquisition (docs/external_screens_inventory.md)
  - 16 non-transformed screens (12 per-guide), 8 anchor sets, 9 expression sources, 36 ORCS non-transformed screens.
  - Control-gene sanity AUROC: ≥0.86 for RPE1/hPSC KO screens. Fails: Haapaniemi RPE1-WT Brunello (0.66) and
    Dräger microglia CRISPRi (0.64).
  - Gaps: iPSC-derived screens are CRISPRi gene-level only; primary blood screens are sort-based (Shifrut, Ting CD34+).
- [~] G5a Screen harmonisation (anchors + controls only) → configs/harmonisation_frozen.json, screen_plan.csv
- [~] G5b Expression harmonisation (matched DepMap dev lines) → configs/expression_harmonisation.json
- [ ] 3. Data pipeline and verification (toy data, controls)
- [ ] 4. Baselines (per-gene mean, per-gene lineage mean)
- [ ] 5. Models and feature ablation
- [ ] 6. Non-transformed generalisation (assemble screens, harmonise, apply)
- [ ] 7. REPORT.md
