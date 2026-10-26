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
- D7 (G1): Anchor/external comparability. With 5-fold CV, each dev line is out-of-fold for exactly one model, so
  "4 non-training models" (DESIGN v3) is impossible. Set A/B skill is computed per fold model and averaged, matching
  the anchors' single out-of-fold model; the 5-model ensemble is reported as secondary.
- D8 (G1): rare_selective gains a ≥1% floor (≥3 lines and 1% ≤ fraction < 5%). Why: ≥3 lines alone admitted 8,281
  mostly-noise genes. Decided after baselines only, before any model comparison. Frozen hashes regenerated.
- D9 (G5b): Expression harmonisation = quantile normalisation to DepMap dev, adopted as a recorded override of the
  rule outcome ("none"). Why: the pre-registered criteria tied (recall@1 0.833 and r ≈0.95 for all candidates);
  MAE favours quantile (0.506 vs 0.537), and quantile corrects library-composition shifts. Decided before any model
  was applied to external data.
- D10 (G5a): Per-screen affine targets = DepMap dev medians of the same control sets (scaling set −1.468, NEG
  0.007), not literal −1/0. Why: the intent was DepMap's scale; our scaling set is deeper than DepMap's control list,
  and literal −1/0 compressed screens 1.47× (anchor agreement 0.164 vs 0.200).
- D11 (G5a): Chronos refit for TKOv3 allowed (time-boxed). It is adopted only if HAP1 anchor agreement rises by ≥0.02
  over median LFC. Caveat: HAP1 is the single TKOv3 anchor.
  Outcome: Chronos raised HAP1 agreement 0.061 → 0.074 (+0.013), so it is **not adopted**. SSMD improved on 10/11.
- D12 (G5a): Wang-library gene-wise correction not applied to Yilmaz 2018. Why: the anchors are published CS gene
  scores, while Yilmaz is count-scored by us, so the pipelines differ. (Rationale restated in D13 as pipeline
  mismatch only; the replicate-agreement drop was an observation, not the criterion.)
- D13 (pre-freeze-transfer): transfer criteria rebuilt after an independent review found the draft was not a true
  pre-registration (wrong MDE quantity, pseudoreplication over 18 correlated RPE1 arms, undefined anchor gate,
  undefined B1 for set B).
  - Study = unit (one consensus profile per study and TP53 status). Hierarchical bootstrap over studies and paralog
    blocks. MDE = 2 × SE of ceiling skill.
  - "Transfers" = lower CI > 0 and skill ≥ 50% of the replicate-ceiling skill. The anchor gate (library-matched)
    is used only to attribute "fails".
  - Q3a primary label: scale-free top-10% rank. Q3b comparator: B1 with a frozen nearest DepMap lineage.
  - Scaling-half control genes are excluded from headline metrics. H1–H3 are operational; H4 is inconclusive by
    declaration.

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
- [x] G0 evaluation infrastructure (metrics, guard, registry, tests; 78 pass) + G1 baselines
  - [x] B0: P1 = 0 exactly (sanity). B1 (κ≈10): P1var 0.184±0.013, P2a 0.906. B2 kNN (k≈20): P1var 0.253±0.009,
    ΔP1 vs B1 +0.069, beats B1 on 81% of variable genes. Runs 20261026-151647_b1, -152737_b2.
  - [x] Red-team review of evaluation code (reviews/redteam_g0/). The independent metric code reproduces the B1
    registry numbers exactly (P1var 0.183737, P2a 0.906274, P2b 0.086748). Fixed: guard bypass via a substitute
    splits file; control-gene reads for test lines; committed audit tables with test-line CRISPR stats (now dev-only;
    the older commit b868f4d still holds them in history; no model used them).
  - [ ] Pending after G2/G5a finish:
    - registry hashes of eval_genes/gene_classes;
    - `freeze-eval` tag;
    - B0/B1/B2 re-run on a clean tagged tree;
    - guarded reader for cache/harmonised.
- [~] G2 M1/M1b kernel ridge, R1, transfer model T, modality ablation, M2 local stage
- [~] Features and priors (src/essml/features.py, priors.py; Nusinow proteomics)
- [x] 6.0 External screen and expression acquisition (docs/external_screens_inventory.md)
  - 16 non-transformed screens (12 per-guide), 8 anchor sets, 9 expression sources, 36 ORCS non-transformed screens.
  - Control-gene sanity AUROC: ≥0.86 for RPE1/hPSC KO screens. Fails: Haapaniemi RPE1-WT Brunello (0.66) and
    Dräger microglia CRISPRi (0.64).
  - Gaps: iPSC-derived screens are CRISPRi gene-level only; primary blood screens are sort-based (Shifrut, Ting CD34+).
- [x] G5a Screen harmonisation (anchors + controls only) → configs/harmonisation_frozen.json, screen_plan.csv
  - QC pass: RPE1 (Olivieri, Brown, Hart, Haapaniemi-null, Bowden) and hPSC (Mair, Ihry, Yilmaz).
    QC fail: Haapaniemi-WT, Shifrut CD8 T, CRISPRi (Tian, Dräger).
  - Anchor deviation agreement: TKO ≈0.09 (gene-wise not applied); Yusa 0.26 and Wang 0.29 after gene-wise.
    Simulated batch undone (RMSE 1.17 → 0.09). Chronos for TKOv3 not adopted (D11).
  - [!] Replicate ceilings, independent same-background screens used as predictors:
    - RPE1: skill over the pan-cancer gene mean ΔAP 0.031 (SE 0.051); Δrecall vs B1-nearest 0.002.
    - hPSC: −0.098 (other hPSC screens predict Mair worse than the cancer gene mean).
    - Cross-study variability in untransformed screens is as large as any cell-type-specific signal. This bounds
      what any model can show for Q3b.
  - [x] Transfer criteria frozen (D13), tag `freeze-transfer`. MDE: RPE1 Q3a 0.101, Q3b 0.162; hPSC 0.082/0.082.
    Nearest lineages frozen: RPE1 → CNS/Brain, hPSC → Testis.
- [x] G5b Expression harmonisation (matched DepMap dev lines) → configs/expression_harmonisation.json
  - 12 exact dev matches (ENCODE via GEO, Tzelepis). Recall@1 0.83 and recall@10 1.0 for every candidate → quantile
    adopted (D9).
  - [!] Applicability domain:
    - inside the cancer distribution: RPE1 (0.18) and HMEC (0.21–0.24), against a dev median of 0.21;
    - at or beyond the edge (94th–100th dev percentile): hESC/iPSC 0.34–0.45 (nearest embryonal carcinoma),
      neurons 0.39–0.40 (neuroblastoma), microglia 0.40–0.44, T cells 0.31–0.48, HSPC 0.35–0.41.
    - Expectation: transfer is most plausible for RPE1/HMEC.
- [ ] 3. Data pipeline and verification (toy data, controls)
- [ ] 4. Baselines (per-gene mean, per-gene lineage mean)
- [ ] 5. Models and feature ablation
- [ ] 6. Non-transformed generalisation (assemble screens, harmonise, apply)
- [ ] 7. REPORT.md
