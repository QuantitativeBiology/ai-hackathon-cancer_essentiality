# Implementation check of accepted review items (RETROSPECTIVE)

**This check is retrospective.** DESIGN §9 required it at G0 and at G4, and it was done at neither. It was carried
out on 2026-10-28, after G4 and G6 had run and been independently reproduced, against HEAD 49b741a. It cannot have
influenced any model, threshold or verdict. Where an item was implemented only after a gate it applied to, the
note says so.

**Method.**
- **Scope:** every item marked accepted (A) or accepted with modification (M) in reviews/round1_response.md,
  round2_response.md and round3_response.md, plus every fixed or partly fixed row of transfer_freeze_response.md.
  Rejected or declined items are listed at the end of each section and were not checked.
- **What counts as evidence:** code, configs, tests, result-file keys, the run registry, git tags and git history.
  DESIGN.md, RESEARCH_TREE.md and docs/ describe intent, so they were used only to locate code. The exception is an
  item that is itself a text change (e.g. "LITERATURE.md amended"), which is verified in that file.
- **Who checked:** three read-only agents, one per block of items, each re-deriving what the review asked from the
  review file itself. The PI spot-checked a sample of their citations (guard.py:113, kernel_ridge.py:48,
  transfer.py:442, the tags, the gate files, the inventory columns) and corrected one finding (R1-ST-15).
- **Hygiene:** no target, CRISPR or expression values were read.
- **Line numbers** refer to HEAD 49b741a.

**Statuses.**
- **implemented:** the cited line does the work, or a test asserts it.
- **partially implemented:** what is missing is named.
- **not implemented:** what was searched for is named.
- **superseded:** a later decision replaced the item; the decision is cited.

**Facts used throughout.**
- **Tags:** `freeze-eval` is at 62aab5b, set after G2/G5a rather than at G0; G2 rows were re-run on the tagged tree.
  `freeze-transfer` is at 987bdec, and `freeze-finalists` at 84efb5c.
- **Gate files:** results/gates/ holds G0, G1, G2, G3, G4 and G6. **There is no G5.json**, and scripts/gate.py
  defines no G5. G5 is evidenced only by the commit and tag of the frozen criteria.
- **Process artefacts:** results/scope_ledger.csv, results/exclusions.csv and this file were first written on
  2026-10-28, in the same commit.

## Round 1 (reviews/round1_response.md)

### External PI
| id | item (what was asked) | verdict | status | evidence | note |
|---|---|---|---|---|---|
| R1-PI-1 | Core-essential test on observed labels, paired with B0; CE recall/precision; lost/gained essentials | A | partially implemented | src/essml/transfer.py:15-17; scripts/evaluate_transfer.py:441-447, 558-575 | no per-screen CE recall/precision or CE-subset AP; lost/gained only as counts/fractions |
| R1-PI-2 | RPE1/non-transformed lines removed from anchors; anchors predicted out-of-fold | A | implemented | scripts/harmonise_screens.py:263-267 (aborts unless every anchor is a dev line); scripts/evaluate_transfer.py:396 | single OOF model per anchor (D7) |
| R1-PI-3 | Pre-registered transfers/fails/inconclusive rule | A | superseded | D13; replacement at src/essml/transfer.py:426-464; configs/transfer_criteria.json; tag freeze-transfer | DiD rule replaced by ceiling-MDE + library-matched anchor gate (D13, D27 m5) |
| R1-PI-4 | Expression-swap negative control; centroid DepMap profile | A | partially implemented | scripts/evaluate_transfer.py:472-484; src/essml/transfer.py:451; tests/test_transfer.py:168 | no centroid-profile prediction |
| R1-PI-5 | Transfer model T; full model with C/M/F neutralised in dev CV; pre-registered lineage mapping | A | partially implemented | src/essml/models/kernel_ridge.py:755-764; configs/runs/t_E_status.json; criteria nearest_depmap_lineage | no dev-CV run of the full model with C/M/F neutralised |
| R1-PI-6 | Ordered work plan with gates; early acquisition; time boxes | A | implemented | scripts/gate.py; results/gates/*.json; scripts/evaluate_transfer.py:948, 986-987 | no G5.json |
| R1-PI-7 | Set A tiers; RPE1 sublines = one background; no per-gene r across sublines | A | partially implemented | configs/transfer_criteria.json set_A.tier_A1; src/essml/transfer.py:496-512 | only A1 operationalised (D21f) |
| R1-PI-8 | Within-lineage per-gene r; ΔP1 vs B1 as headline | A | implemented | src/essml/evaluate.py:286-288, 311-314; scripts/evaluate_test.py:419-466 | |
| R1-PI-9 | P1 per fold then averaged; paired 1-SE | A | implemented | src/essml/evaluate.py:240-246; scripts/g2_report.py:73-85, 224-245 | D3 |
| R1-PI-10 | P2a as Δ vs B0 | A | implemented | src/essml/evaluate.py:289-291 | |
| R1-PI-11 | Rank-only metrics for BF screens; CRISPRi descriptive | A | implemented | scripts/harmonise_screens.py:663-669; configs/screen_plan.csv | |
| R1-PI-12 | Applicability-domain distance per external profile | A | implemented | src/essml/expr_harmonise.py:517-543; results/expression_harmonisation/profiles.csv | not joined to G6 outputs |
| R1-PI-13 | User download list; Nusinow proteomics | A | implemented | docs/USER_DOWNLOAD_REQUEST.md; src/essml/features.py:235; registry m1_E_Psub, m1b_E_P_Psub | |

### Statistics
| id | item | verdict | status | evidence | note |
|---|---|---|---|---|---|
| R1-ST-1 | P1 within fold, never pooled; unit test B0 = 0 | A | implemented | src/essml/evaluate.py:240-246; tests/test_evaluate.py:32-44 | |
| R1-ST-2 | Scale-aware skill 1−SSE/SSE_B0 and r of residuals from B1 | A | partially implemented | scripts/evaluate_test.py:419-443 (lineage/library-centred r) | no 1−SSE/SSE_B0 skill; centred r is a proxy only |
| R1-ST-3 | Fold-wise classes for training; frozen dev class list before unblinding | A | implemented | configs/gene_classes.csv (5181e51, before 84efb5c); kernel_ridge.py:24 | gene_classes.csv is sha256-recorded per run but not in configs/frozen_hashes.json |
| R1-ST-4 | Normalised AP + AUROC; prevalence ≤ 0.5 headline; ties; constant-predictor test | A | implemented | src/essml/metrics.py:100-127; tests/test_metrics.py:40 | |
| R1-ST-5 | Per-line ΔAP vs B0, paired | A | implemented | src/essml/evaluate.py:289-291 | |
| R1-ST-6 | Paired 1-SE, complexity order, δ = 0.005 | A | implemented | scripts/g2_report.py:27, 210-245; scripts/g3_report.py:101-106 | |
| R1-ST-7 | One named headline before unblinding | A | implemented | configs/finalists.json; tag freeze-finalists 84efb5c | |
| R1-ST-8 | Declared ablation family; patient-bootstrap CIs; Holm; "helps" needs CI > 0 and \|Δ\| > δ | A | partially implemented | scripts/g2_report.py:45-53, 88-98, 181-191 | Holm on fold-paired t-tests, not a patient bootstrap; no \|Δ\| > δ condition |
| R1-ST-9 | No gene-independent tests; permuted-null calibration; lineage-stratified patient bootstrap (2000) | A | implemented | src/essml/heldout.py:211; scripts/evaluate_test.py:38, 408-416 | null shown beside results, not subtracted |
| R1-ST-10 | Power statement; per-lineage only for ≥ 20 test lines | A | implemented | scripts/evaluate_test.py:100, 445-466 | |
| R1-ST-11 | Set A: line as unit, percentile within cancer test | A | implemented | src/essml/transfer.py:489-512; scripts/evaluate_transfer.py:731-754 | gene-block bootstrap (one patient) |
| R1-ST-12 | Operational "context-specific survive": recall at fixed NEG FPR + nulls | A | implemented | src/essml/transfer.py:18-21; src/essml/metrics.py:130-148 | comparator B1-nearest (D13) |
| R1-ST-13 | Affine validation by medians/SSMD; 20 seeded half-splits; CE outside CEGv2 | A | partially implemented | src/essml/harmonise.py:183-230; scripts/evaluate_transfer.py:506-513 | one fixed-seed split only, no 20-split spread |
| R1-ST-14 | Anchor agreement on deviations; EB shrinkage; LOAO; dev-only assertion | A | implemented | src/essml/harmonise.py:234-252, 339-376; scripts/harmonise_screens.py:255-267 | LOAO from ≥ 2 anchors, not ≥ 3 |
| R1-ST-15 | Near-duplicate test–dev flag + sensitivity; dev-only fusion filter and quantile reference; training-fold NaN fill | A | partially implemented | src/essml/features.py:441-454; src/essml/expr_harmonise.py:19; near-duplicate check only in reviews/g4_review/REPORT.md:232 | the near-duplicate check exists only in the independent G4 reviewer's analysis; no pipeline output and no sensitivity run (PI correction) |
| R1-ST-16 | Disclose the Chronos joint fit | A | implemented | DESIGN.md:24-26; D26 item 5 | final disclosure belongs in REPORT.md (G7, not yet written) |

### Prior work
| id | item | verdict | status | evidence | note |
|---|---|---|---|---|---|
| R1-PW-1 | Rung R: ridge, EN, RF, TrPLet-like (nested + non-nested) | M | partially implemented | registry r1_E, g3_r2_sub150, g3_r4_sub500, g3_r4k_sub500, g3_r4_oracle_sub500 (leaky=1); src/essml/models/published.py | R3 random forest never coded or run; R2 on 150 genes (D19) |
| R1-PW-2 | RRR / matrix completion as diagnostic | M | superseded | reviews/round3_response.md must 8 (cut); DESIGN §5 | |
| R1-PW-4 | Published-comparable block | A | implemented | scripts/evaluate_test.py:468-480; results/g3_nonlinear_summary.json; scripts/g3_published_genesets.py | |
| R1-PW-5 | Correct novelty claims; restore within-lineage r | A | implemented | LITERATURE.md:28-29 (48bd50d); scripts/evaluate_test.py:419-443 | |
| R1-PW-6 | Library-stratified reporting via NaN-group proxy; ScreenMap requested | M | implemented | scripts/make_library_groups.py:28; configs/library_groups.csv; docs/USER_DOWNLOAD_REQUEST.md | |
| R1-PW-7 | Expression alignment variants; E-only transfer model; kNN transfer baseline | A | implemented | src/essml/expr_harmonise.py:18-21; scripts/evaluate_transfer.py:88 | Celligner deferred (round 3) |
| R1-PW-8 | Not-expressed hinge, deep deletion, paralog mutation, family size/identity | A | implemented | src/essml/priors.py:143-159, 229-308, 338-340; src/essml/models/m2_local.py:46-47 | family size/identity only feed M4, which was dropped (D19) |
| R1-PW-9 | Proteomics ablation arm | A | implemented | registry m1_E_Psub, m1b_E_P_Psub | |

Rejected: R1-PW-3 (concatenated OOF P1 for selection; kept only in the published-comparable block).

### Data
| id | item | verdict | status | evidence | note |
|---|---|---|---|---|---|
| R1-DA-1 | No mean-fill; pattern-wise masked fits; ≥ 50/≥ 40 observed; report by NaN group | A | implemented | src/essml/models/kernel_ridge.py:271-298; scripts/make_eval_genes.py:16; tests/test_kernel_ridge.py:41 | |
| R1-DA-2 | Provenance fields ablation-only; age missing flag; check that the model predicts the NaN group | A | partially implemented | src/essml/features.py:68-69, 361-381, 530-531; tests/test_features.py:252 | NaN-group predictability check not found (within-library P1, D16, partly covers it) |
| R1-DA-3 | E-only transfer model | A | implemented | src/essml/models/kernel_ridge.py:755-764 | |
| R1-DA-4 | CN transform, clip, all-NaN drop, absent flag, autosomes | A | implemented | src/essml/features.py:74-92, 369; tests/test_features.py:13-32 | |
| R1-DA-5 | Fusion filter (self, non-coding, ≥ 5 lines) | A | implemented | src/essml/features.py:109-155, 374-375; tests/test_features.py:40 | |
| R1-DA-6 | TPM on all genes then subset; QuantSeq flagged; same-line validation | A | implemented | src/essml/expr_harmonise.py:8-19, 201-217; tests/test_expr_harmonise.py:49-83 | |
| R1-DA-7 | Set A per background; verify Hs 729; viral lines out of headline | A | partially implemented | configs/transfer_criteria.json set_A.tier_A1 | Hs 729 origin never verified; only A1 evaluated (D21f) |
| R1-DA-8 | M as 0/1/2 dosage; ≥ 5 training lines; log burden | A | implemented | src/essml/features.py:29-30, 94-96, 371-373 | |
| R1-DA-9 | Entrez/HGNC mapping; coverage; RPE1 with/without CEGv2; hu.MAP fallback | A | partially implemented | src/essml/priors.py:58, 218-225; results/features_build_report.json; scripts/evaluate_transfer.py:506-513 | neither hu.MAP nor CORUM used |
| R1-DA-10 | Dev-only anchors, OOF; HepG2/MCF-7 expression-only; shrinkage applied only if LOAO improves | A | implemented | scripts/harmonise_screens.py:263-267, 353-398; src/essml/harmonise.py:55-58, 234-252 | |

### Algorithms & compute
| id | item | verdict | status | evidence | note |
|---|---|---|---|---|---|
| R1-AC-1 | eigh + closed-form LOO, wide λ grid, grouped λ | A | implemented | src/essml/models/kernel_ridge.py:48, 181-232, 271-359, 657-661; tests/test_kernel_ridge.py:61, 111 | |
| R1-AC-2 | Leave-line-out lineage residuals; κ in grouped CV | A | implemented | src/essml/baselines.py:94; kernel_ridge.py:385, 645-652; tests/test_baselines.py:42 | model of record later uses the B0 residual (D14) |
| R1-AC-3 | M2 on OOF M1, batched normal equations | A | implemented | src/essml/models/m2_local.py:75-112, 160-180, 195-210; tests/test_m2_nested.py:37 | |
| R1-AC-4 | Factorised M4 behind an MLP-vs-ridge pre-check | A | implemented | src/essml/models/factorised.py:166-200; configs/runs/g3_m4pre_* | pre-check failed, full M4 never ran (D19) |
| R1-AC-5 | Nested stacking, non-negative per-class weights, plain average baseline | A | implemented | src/essml/models/m5_stack.py:94-128, 214-245; tests/test_m5_stack.py:84, 105 | |
| R1-AC-6 | NaN handling: pattern-wise exact fits | M | implemented | kernel_ridge.py:101, 271; tests/test_kernel_ridge.py:41; tests/test_g3_models.py:83, 97 | |
| R1-AC-7 | Thread env vars before numpy; threadpoolctl; 1-thread workers; main guards | A | implemented | scripts/run.py:61-82, 93, 420-422; src/essml/models/lgbm.py:39 | |
| R1-AC-8 | Numerics (float64 eigh, clip, intercept, constant r = 0) | A | implemented | kernel_ridge.py:196-197, 619-620, 639-640; src/essml/metrics.py:54; tests/test_metrics.py:40-65 | |
| R1-AC-S1 | M1b per-modality kernels | A | implemented | kernel_ridge.py:472-550, 624-643; configs/runs/m1b_*, m1bt_* | |
| R1-AC-S2 | Reduced-rank regression diagnostic | A | superseded | reviews/round3_response.md must 8 (cut) | |
| R1-AC-S3 | M3 go/no-go with in-fold screening | A | implemented | src/essml/models/lgbm.py:24, 39, 208; configs/runs/g3_m3n.json | |
| R1-AC-S4 | Runtime/memory logged by the harness | A | implemented | results/registry.csv fit_s, predict_s, total_s, peak_rss_mb | |
| R1-AC-S5 | λ re-tuned per ablation | A | implemented | kernel_ridge.py:686-688 | |

### Agent failure modes
| id | item | verdict | status | evidence | note |
|---|---|---|---|---|---|
| R1-AF-1 | Frozen metric modules hashed, failing test | A | implemented | configs/frozen_hashes.json (first in 5181e51); tests/test_frozen.py; tag freeze-eval 62aab5b | tag set after G2, not at G0; G2 re-run on the tagged tree |
| R1-AF-2 | Test guard: env unlock + committed finalist hash; log every access | A | implemented | src/essml/guard.py:58-60, 113-138, 202-222; tests/test_guard.py:42-106 | |
| R1-AF-3 | Set A/B models frozen before prediction | A | implemented | guard.py:15-19; tag freeze-finalists (84efb5c) precedes G6 stage A (c5b2ec8) | fold-model hashes recorded only after the fact (results/transfer/g6_input_hashes.json, D27) |
| R1-AF-4 | Harmonisation fixed blind and committed before set-B predictions | A | implemented | configs/harmonisation_frozen.json; tag freeze-transfer; guard.py:310-352; tests/test_no_direct_cache_harmonised.py:72 | one unguarded column-name read (D21d) |
| R1-AF-5 | Numeric pre-registration of verdicts and failure attribution | A | implemented | configs/transfer_criteria.json verdicts/hypotheses; src/essml/transfer.py:426-464 | rebuilt in D13 |
| R1-AF-6 | Registry by construction; numbers traced; dirty runs not citable | A | partially implemented | scripts/run.py:51-52, 353, 436, 535-547; scripts/check_report_numbers.py | no seed or training-ID hash column; the checker does not refuse dirty-tree runs; REPORT.md not yet written |
| R1-AF-7 | Independent reproduction with own metric code, 1e-6 tolerance, blinded prompts stored, red team per finding | A | partially implemented | reviews/redteam_g0, g2_review, g4_review, g6_review (own code, e.g. g6_review/reproduce_g6.py) | no docs/review_prompts/; G4 tolerance 1e-4 (g4_review/compare.py:13); red teams per gate, not per finding |
| R1-AF-8 | Probe numbers quarantined; probe bugs become tests | A | implemented | RESEARCH_TREE tree 1.3 [-]; tests/test_run_leakage.py:1-3, 62 | |
| R1-AF-9 | Stopping rules | A | implemented | applied in D19 and D23f | process rule; not enforced in code |
| R1-AF-10 | Literature verification policy; re-tag lit files; intake records; Pacini fix | A | partially implemented | LITERATURE.md:4, 66; configs/external_screens_manifest.json; DATA.md | docs/lit/harmonisation_labels.md and data_resources.md not re-tagged |
| R1-AF-11 | Post-freeze bug protocol | A | implemented | baselines re-run 54df618; scripts/evaluate_transfer.py:92-96 (KNOWN_POST_FREEZE, D21c) | process rule |
| R1-AF-12 | Exclusion accounting per filter; paired comparisons assert identical sets | A | partially implemented | results/harmonisation/exclusions.csv; results/features_build_report.json; src/essml/evaluate.py:62-63 | global results/exclusions.csv did not exist at G0/G4; written retrospectively 2026-10-28 |
| R1-AF-13 | Fixed status block (hashes, access count, open defects, pending reviews) | A | partially implemented | RESEARCH_TREE.md:5 | narrative log without the fixed fields |

### Biology
| id | item | verdict | status | evidence | note |
|---|---|---|---|---|---|
| R1-BI-1 | Rare-selective class | A | implemented | src/essml/genesets.py:109; configs/gene_classes.csv rare_selective | 1% floor (D8) |
| R1-BI-2 | Set A tiers by transformation mechanism | A | partially implemented | scripts/transfer_criteria.py:399; criteria set_A.tier_A1 | no A2/A3 assignment in code or configs; only A1 evaluated (D21f) |
| R1-BI-3 | TP53 hypothesis split (learnable vs diagnostic) | A | implemented | criteria hypotheses.H1 / H1_prime; src/essml/transfer.py:515-528 | H1 scored by sign, not rank/slope (D13) |
| R1-BI-4 | Isogenic TP53 counterfactual (H2) | A | implemented | scripts/evaluate_transfer.py:577-606; src/essml/transfer.py:531-547 | |
| R1-BI-5 | Profile-swap null + nearest-lineage comparator | A | implemented | scripts/evaluate_transfer.py:358, 472-483; src/essml/transfer.py:439-449 | |
| R1-BI-6 | Recall at fixed NEG FPR as headline | A | implemented | src/essml/transfer.py:44-46, 114, 173 | |
| R1-BI-7 | EC leave-out for pluripotency (H3) | A | implemented | configs/runs/t_E_status_noEC.json; src/essml/transfer.py:550-566 | 4 EC dev lines (D20d) |
| R1-BI-8 | Post-mitotic divergence split (proliferation vs maintenance); fitness-only arms | A | partially implemented | fitness rule: results/harmonisation/exclusions.csv (X_* rows) | no proliferation/maintenance split anywhere in code |
| R1-BI-9 | KO↔CRISPRi anchors if obtainable, else CRISPRi descriptive | M | implemented | configs/screen_plan.csv CRISPRi rows ("no K562 CRISPRi anchor obtained") | fallback branch; QC only, no descriptive CRISPRi output exists |
| R1-BI-10 | ENCODE-matched expression correction; rank variant; state matching | A | partially implemented | src/essml/expr_harmonise.py:1-30; configs/expression_harmonisation.json (D9) | rank only as a harmonisation candidate, not a transfer model; config `status` still "draft" |
| R1-BI-11 | Short-screen QC gating | A | implemented | src/essml/harmonise.py:183-215; results/harmonisation/screen_qc.csv | |
| R1-BI-12 | p53/Rb composites incl. viral-oncoprotein flag | A | partially implemented | src/essml/features.py:170-212 | no viral-oncoprotein flag |
| R1-BI-13 | Oncogene addiction as negative control; CN-collateral null | A | implemented | scripts/evaluate_transfer.py:97-98, 630-643 | |
| R1-BI-14 | No per-gene r across sublines; pooled T1+T2; aneuploidy contrast | A | partially implemented | criteria set_A.tier_A1.statistics | no pooled T1+T2, no aneuploidy contrast |
| R1-BI-15 | Haematopoietic positives, unexpressed-gene FPR, declared TCR failures | A | partially implemented | criteria hypotheses.H4; src/essml/transfer.py:569 | H4 declared inconclusive; no unexpressed-gene FPR control in code |
| R1-BI-16 | Missing-modality convention for primary cells; karyotype; zero-fill vs asserted-WT | A | not implemented | searched configs/, scripts/, src/ for karyotype / asserted-WT | convention only in DESIGN §2; moot for T, relevant to the m5_avg secondary (D25d) |
| R1-BI-17 | Record Cas9 delivery mode as a covariate | A | not implemented | docs/external_screens_inventory.md has no delivery column; nothing in harmonisation_screens.json or screen_plan.csv | the response says it is recorded in the inventory; it is not |

## Round 2 (reviews/round2_response.md)
| id | item | verdict | status | evidence | note |
|---|---|---|---|---|---|
| R2-CONFLICT | Screen-specific: fraction ≤ 0.5 headline, ≤ 5% secondary | adopted | implemented | src/essml/transfer.py:44-45, 179; scripts/evaluate_transfer.py:83-85 | D4 |
| R2-MOD-1 | Tiers from metadata; functional p53 check descriptive only | adopted (functional tiering rejected) | implemented | scripts/evaluate_transfer.py:762-765; setA summary p53_functional_check | D6 |
| R2-MOD-2 | Test uses all-dev refit; set A/B use the 5 fold models | decided | implemented | scripts/evaluate_transfer.py:139-162; src/essml/transfer.py:312-319 | skills averaged per fold model (D7) |
| R2-§2-a | CN transform log2((CN+0.1)/1.1) | adopted | implemented | src/essml/features.py:60, 74-80; tests/test_features.py:13 | |
| R2-§2-b | Instability on raw \|CN−1\| > 0.3 | adopted | implemented | src/essml/features.py:62, 83-91; tests/test_features.py:24 | |
| R2-§2-c | Expanded p53/Rb composites | adopted | implemented | src/essml/features.py:66-67, 170-205 | |
| R2-§3-a | A1 labelled a single-background case study | adopted | implemented | criteria set_A.tier_A1 | |
| R2-§4-a | Fixed eval_genes.csv; mask-only test counts | adopted | implemented | scripts/make_eval_genes.py:20-27; src/essml/guard.py:226; tests/test_guard.py:141 | |
| R2-§4-b | Prevalences over observed lines | adopted | implemented | scripts/make_gene_classes.py:26-29 | |
| R2-§5M1-a | NaN-pattern merging | adopted | implemented | kernel_ridge.py:101-155, 271-298; tests/test_kernel_ridge.py:84 | |
| R2-§5M1-b | Patterns < 50 rows get the offset only | adopted | implemented | kernel_ridge.py:145-147, 272 | |
| R2-§5M1-c | Pattern sub-kernels of one fold K | adopted | implemented | kernel_ridge.py:19-20, 304; tests/test_kernel_ridge.py:41 | |
| R2-§5M1-d | M1b weight tuning (≤ 25 eigh; defaults reported) | adopted | implemented | kernel_ridge.py:468-552, 678-706; tests/test_kernel_ridge_solver.py:106 | D15, D20a |
| R2-§5M1-e | R2 time box 2 h, fallback 500 genes | adopted | superseded | D19 / D-G3-3: registry g3_r2_probe killed (time box), g3_r2_sub150 ran | 150-gene fallback replaced the 500-gene one |
| R2-§5T-a | No lineage in primary T | adopted | implemented | kernel_ridge.py:755-765 | |
| R2-§5-a | Frozen-model protocol; T and no-EC frozen at G4 | adopted | implemented | configs/finalists.json (84efb5c) | |
| R2-§7.1-a | K562-only CRISPRi anchor | adopted | not implemented | screen_plan.csv CRISPRi rows: anchors "none"; harmonisation O4 | data not obtained; CRISPRi stays QC-only |
| R2-§7.1-b | Neuron proliferation vs maintenance split | adopted | not implemented | no gene-set split in scripts/ or src/ | CRISPRi stratum never evaluated at G6 |
| R2-§7.1-c | Committed screen plan | adopted | implemented | configs/screen_plan.csv (083578a); sha256 in criteria inputs | |
| R2-§7.2-a | Location-only below 8 anchors | adopted | implemented | src/essml/harmonise.py:56-57, 282; tests/test_harmonise.py:145 | |
| R2-§7.2-b | Leave-anchor-out anchor skill | adopted | implemented | src/essml/harmonise.py:320-335; tests/test_harmonise.py:177 | anchors only gate "fails" (D13) |
| R2-§7.2-c | p53/DDR genes excluded from correction | adopted | implemented | src/essml/harmonise.py:67-74, 259-268, 305-318; tests/test_harmonise.py:157 | before/after H1 report moot: no gene-wise correction on any RPE1/hPSC verdict screen |
| R2-§7.3-a | Celligner-style candidate | adopted | not implemented | configs/expression_harmonisation.json ("DEFERRED"); src/essml/expr_harmonise.py:21 | deferral accepted in round 3, never built |
| R2-§7.3-b | GTEXDEPMAP cited | adopted | implemented | DESIGN.md:375-376 | text item |
| R2-§7.45-a | Primary statistics named | adopted | implemented | src/essml/transfer.py:16-20, 163-180; scripts/evaluate_transfer.py:380-384, 462-469 | D13 form (rank label; B1-nearest) |
| R2-§7.45-b | B2 comparator | adopted | implemented | scripts/evaluate_transfer.py:88, 525-543; *_secondary.json *_B2 | B2 scored against B0/B1n; T−B2 not computed directly |
| R2-§7.45-c | USP7 dropped from H1 | adopted | implemented | criteria hypotheses.H1.excluded | |
| R2-§7.45-d | CN null expectation | adopted | implemented | scripts/evaluate_transfer.py:97-98, 630-645 | |
| R2-§7.45-e | H4 re-anchored on Shifrut T cells | adopted | superseded | D13: H4 inconclusive by declaration (Shifrut fails QC); src/essml/transfer.py:569-570 | |
| R2-§7.6-a | "Transfers" defined | adopted | implemented | src/essml/transfer.py:442; tests/test_transfer.py:151 | ≥ 50% of the replicate ceiling (D13) |
| R2-§7.6-b | MDE condition | adopted | implemented | src/essml/transfer.py:411, 443 | MDE from ceiling skill (D13) |
| R2-§7.6-c | Block bootstrap | adopted | implemented | src/essml/transfer.py:65-75, 241-270; tests/test_transfer.py:212 | paralog families only (D13.6) |
| R2-§7.6-d | Count (majority) rule | adopted | superseded | D13/D21f: every verdict background must agree; src/essml/transfer.py:466-480 | |
| R2-§7.6-e | Swap-null cap | adopted | implemented | src/essml/transfer.py:444, 450-453; tests/test_transfer.py:168 | |
| R2-§9-a | Guard bound to finalist hashes and a clean tree | adopted | implemented | src/essml/guard.py:113-137; tests/test_guard.py:53, 60 | |
| R2-§9-b | Set A in the guard | adopted | implemented | guard.py:7-8, 129 | |
| R2-§9-c | Static test against direct CRISPR reads | adopted | implemented | tests/test_no_direct_crispr_reads.py:12-34 | scans src/ and scripts/ only (D26.8) |
| R2-§9-d | Freeze tags + test that the diff against the tag is empty | adopted | partially implemented | tags exist; tests/test_frozen.py:10-14 (hashes) | no diff-vs-tag test; harmonisation scripts not hash-checked; post-freeze read-path change to harmonise_screens.py (8a613f0) has no D-entry (D21c names transfer_criteria.py only) |
| R2-§9-e | Harmonisation freeze file required before set-A/B access | adopted | implemented | configs/harmonisation_frozen.json; guard.py:52, 129; scripts/evaluate_transfer.py:947 | file `_status` still "DRAFT" |
| R2-§9-f | Gate evidence files; each gate checks the previous | adopted | partially implemented | results/gates/G0–G4, G6; scripts/gate.py:61, 117, 163, 203 | no G5.json; G6 requires G4, not G5 |
| R2-§9-g | Registry row at run start | adopted | implemented | scripts/run.py:431-438, 535, 547 | |
| R2-§9-h | Registry columns (dirty, splits hash, roles, training-ID hash) + leaky flag | adopted | partially implemented | registry header; scripts/run.py:433-436, 464, 482-489 | no training-ID hash (count only); dirty runs not refused by check_report_numbers.py |
| R2-§9-i | Scope ledger | adopted | not implemented | absent from git history before 2026-10-28 | written retrospectively now: results/scope_ledger.csv |
| R2-§9-j | Implementation check at G0 and G4 | adopted | not implemented | absent from git history before 2026-10-28 | this file, retrospective |
| R2-FIX-a | LITERATURE.md Pacini statement | fixed | implemented | LITERATURE.md:66 | |
| R2-FIX-b | Untagged-note warning | fixed | implemented | LITERATURE.md:4 | |
| R2-FIX-c | Stale probe numbers removed | fixed | implemented | RESEARCH_TREE.md:149, 495 | |

## Round 3 (reviews/round3_response.md)
| id | item | verdict | status | evidence | note |
|---|---|---|---|---|---|
| R3-MUST-1 | One frozen dev-wide gene_classes.csv for reporting; fold-wise where it drives training | M | implemented | configs/gene_classes.csv (5181e51); scripts/run.py:45, 435; kernel_ridge.py:657-661 | |
| R3-MUST-2 | eval_genes_external | A | implemented | src/essml/transfer.py:156; scripts/evaluate_transfer.py:745 | |
| R3-MUST-3 | Guard exceptions (control genes; per finalist × screen) | A | implemented | src/essml/guard.py:15-24, 102-138, 263-352; tests/test_guard.py:106, 121 | |
| R3-MUST-4 | Anchors predicted by the 4 non-training fold models | A | superseded | D7; src/essml/transfer.py:299-304 | impossible with 5-fold CV |
| R3-MUST-5 | Thresholds committed in transfer_criteria.json at G5 | A | implemented | criteria harmonisation/bootstrap/hypotheses; eb23a79, tag freeze-transfer, before G6 | no G5.json |
| R3-MUST-6 | Verdicts per biological background | A | implemented | criteria backgrounds; src/essml/transfer.py:426, 466 | |
| R3-MUST-7 | Scale-free rank label; absolute AUPRC vs ceiling as primary Q3a evidence | A | superseded | D13: primary is ΔAP vs B0 judged against ceiling skill; rank label at src/essml/transfer.py:42, 163 | absolute AUPRC kept as a secondary |
| R3-MUST-8b | Cut reduced-rank / low-rank diagnostics | A | implemented | DESIGN §5; no RRR config | (cutting M4 was rejected) |
| R3-SH-1 | Median guide LFC default; Chronos only where needed | A | implemented | configs/harmonisation_frozen.json scoring; D11 | |
| R3-SH-2 | Set A before set B | A | implemented | scripts/evaluate_transfer.py:987, 1292; commits c5b2ec8 then f348ffb | |
| R3-SH-3 | Defer M1b tuning, R3, Celligner | A | implemented | no R3 config; src/essml/expr_harmonise.py:22; M1b tuning ran later (D15) | |
| R3-SH-4 | Per-gate compute budgets | A | partially implemented | DESIGN §12 time boxes for G3 only; per-run boxes in configs | no per-gate budgets or total |

Rejected: cutting M4 (must 8a), dropping the proteomics arm, the 2-SE selection rule.

## Transfer-freeze review (reviews/transfer_freeze_response.md)
| id | item | response | status | evidence | note |
|---|---|---|---|---|---|
| TF-F1 | MDE = 2 × SE of ceiling skill | Fixed | implemented | src/essml/transfer.py:368-414; scripts/transfer_criteria.py:206 | |
| TF-F2 | Hierarchical bootstrap (studies → blocks → predictor pools) | Fixed | implemented | src/essml/transfer.py:47, 241-281, 400-405; tests/test_transfer.py:190, 202, 248 | |
| TF-F3 | One study-level population; Bowden H2-only; minimum n | Fixed | implemented | src/essml/transfer.py:183-187, 574-595; criteria population.minimum_n; tests/test_transfer.py:226 | |
| TF-F5 | Anchor statistic written out; library-matched; gate attributes "fails" only | Partly fixed | implemented | criteria verdicts.anchor_gate; scripts/evaluate_transfer.py:387-417; src/essml/transfer.py:418-422, 443 | accepted part only |
| TF-F6 | One primary per question; secondaries marked; IUT | Fixed | implemented | src/essml/transfer.py:42-46; scripts/evaluate_transfer.py:80-87 | |
| TF-F7 | Scaling-half genes dropped everywhere | Fixed | implemented | src/essml/transfer.py:151-160; scripts/evaluate_transfer.py:410, 745; tests/test_transfer.py:77 | |
| TF-F8 | B1-nearest comparator, hashed; B0 secondary | Fixed | implemented | scripts/transfer_criteria.py:37, 146-164, 182 | |
| TF-F9 | Swap set enumerated, "any" aggregation, matched-minus-swap secondary | Partly declined | implemented | criteria expression_swap_profiles; scripts/evaluate_transfer.py:472-485 | accepted part only |
| TF-F10.1-5 | Weighting, minimum n, boundaries, tie rules, overall/mixed | Fixed | implemented | src/essml/transfer.py:52-54, 210-218, 284-296, 442-449, 466-480; tests/test_transfer.py:151, 298 | |
| TF-F10.6 | AP tie handling = sklearn; constant AP = π | Fixed | implemented | src/essml/transfer.py:79-147; tests/test_transfer.py:52, 71 | |
| TF-F11 | Criteria void without committed finalists; tag split in two | Fixed by declaration | partially implemented | finalists committed 84efb5c before G6; scripts/evaluate_transfer.py:946 | tags freeze-transfer-criteria / freeze-transfer-models never created (D27 m4) |
| TF-F12 | Per-fold skill averaging, ensemble secondary | Partly fixed | implemented | src/essml/transfer.py:299-319; tests/test_transfer.py:278 | |
| TF-F13 | H1 RPE1 two-point counterfactual; sanity-check label; hPSC descriptive | Fixed | partially implemented | src/essml/transfer.py:515-528; scripts/evaluate_transfer.py:577-588 | the promised descriptive hPSC H1 is in neither criteria nor G6 output |
| TF-F14 | H2 ≥ 10 genes; Spearman + sign; raw secondary | Fixed | implemented | src/essml/transfer.py:531-547; scripts/evaluate_transfer.py:589-606 | |
| TF-F15 | H3 observed-side condition; TDGF1 excluded | Partly fixed | implemented | src/essml/transfer.py:550-566; criteria H3.not_in_depmap | accepted part only |
| TF-F16 | Harmonisation narrative corrected | Fixed | partially implemented | docs/harmonisation_validation.md §7 corrected | §10 of the same doc still shows the pre-D13 draft rules, unmarked |
| TF-F17 | D12 restated as pipeline mismatch; Yilmaz predictor-only | Fixed | implemented | configs/harmonisation_frozen.json:495; criteria hPSC yilmaz2018 role | |
| TF-F18a | Chronos joint fit recorded | recorded | implemented | criteria decisions.D11 | |
| TF-F18b | Remove D9's Yilmaz citation | open (left to PI) | not implemented | configs/expression_harmonisation.json:58 | never edited |
| TF-F19 | Representative arm by SSMD | Moot | implemented | src/essml/transfer.py:574-589 (all QC-passing arms enter) | moot |
| TF-E-1 | `_status` flipped at the freeze | fixed | partially implemented | criteria `_status` | still "PRE-FREEZE" and names tags that do not exist; harmonisation_frozen.json `_status` "DRAFT" |
| TF-E-2 | D8, D9, D11 recorded | fixed | implemented | criteria decisions | |
| TF-E-3 | Inputs sha256-hashed | fixed | implemented | criteria inputs_sha256; scripts/evaluate_transfer.py:951-962 | |
| TF-E-4 | QC per arm and per consensus | fixed | implemented | criteria backgrounds.*.harmonisation_validation; scripts/evaluate_transfer.py:436-439 | |
| TF-E-5 | Guard API for the set-B write path | Not done | not implemented | src/essml/guard.py has only the reader load_harmonised (310) | acknowledged in the response |

Declined: F4 (complexes in gene blocks); the declined parts of F5, F9 and F15.

The response's three "points for the PI" are disclosures, not review items, and are excluded from the counts. All
three are carried in the outputs:
- the ratio-condition-vacuous flag (src/essml/transfer.py:445-446);
- the fixed-pool MDE sensitivity (scripts/transfer_criteria.py:223-224);
- the 4 RPE1 resampling units (criteria RPE1.studies).

## Summary

| block | items checked | implemented | partially implemented | not implemented | superseded |
|---|---|---|---|---|---|
| Round 1 | 90 | 62 | 23 | 2 | 3 |
| Round 2 | 47 | 36 | 3 | 5 | 3 |
| Round 3 | 12 | 9 | 1 | 0 | 2 |
| Transfer freeze | 25 | 19 | 4 | 2 | 0 |
| **Total** | **174** | **126** | **31** | **9** | **8** |

**Not implemented (9):**
1. R1-BI-16: missing-modality convention for primary cells; karyotype; zero-fill vs asserted-WT comparison.
2. R1-BI-17: Cas9 delivery mode as a covariate. The response claims the inventory records it; it does not.
3. R2-§7.1-a: K562 CRISPRi anchor. Not obtained (open item O4).
4. R2-§7.1-b: neuron proliferation vs maintenance split. The CRISPRi stratum was never evaluated.
5. R2-§7.3-a: Celligner-style expression candidate. Deferred in round 3 and never built.
6. R2-§9-i: scope ledger at G0/G4. Written retrospectively now.
7. R2-§9-j: implementation check at G0/G4. This file, retrospective.
8. TF-F18b: D9's Yilmaz citation in configs/expression_harmonisation.json. Never edited.
9. TF-E-5: guard API for the set-B harmonised write path.

**Partially implemented (31):**
- **Round 1 (23):** PI-1, PI-4, PI-5, PI-7, ST-2, ST-8, ST-13, ST-15, PW-1, DA-2, DA-7, DA-9, AF-6, AF-7, AF-10,
  AF-12, AF-13, BI-2, BI-8, BI-10, BI-12, BI-14, BI-15.
- **Round 2 (3):** §9-d, §9-f, §9-h.
- **Round 3 (1):** SH-4.
- **Transfer freeze (4):** F11, F13, F16, E-1.

**Superseded (8):**
- **Round 1:** PI-3 (D13), PW-2 and AC-S2 (round-3 cut).
- **Round 2:** §5M1-e (D19), §7.45-e (D13), §7.6-d (D13/D21f).
- **Round 3:** MUST-4 (D7), MUST-7 (D13).

The gaps that recur across reviews:
- no G5 gate file;
- only set-A tier A1 evaluated (D21f);
- no descriptive CRISPRi analysis;
- the process artefacts of DESIGN §9 were written only now;
- frozen-file `_status` fields never flipped;
- no training-ID hash in the registry, and the "dirty runs not citable" rule is not enforced.
