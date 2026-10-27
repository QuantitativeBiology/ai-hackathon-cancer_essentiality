# Research tree — CRISPR essentiality prediction (DepMap 26Q3)

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding

## Status summary (latest first)

### 2026-10-26 — M3: first model rung (G2) done
- **Result:** the kernel-ridge rung (G2) more than doubles the lineage-mean baseline's per-gene skill, but it is not
  better than a plain ridge regression on expression.
- **Numbers:** dev-CV P1 (variable genes) 0.308 vs B1 0.184 and B2 (expression kNN) 0.253. P2b nAP 0.19 vs 0.087.
- **Cost:** about 80–100 s per fold at 12 threads.
- **Residual target:** the pre-registered B1-residual target was worse than the B0 residual. I switched it (D14).
- **Modality ablation:** with default block weights, every non-expression modality hurts. Copy number carries real
  signal in leave-one-out, so block-weight tuning comes before any negative claim (D15).
- **Controls:** the permuted-feature control reproduces B1 exactly. Provenance metadata behaves as a batch proxy and
  stays excluded.
- **Reproduction:** 19 configs were re-run independently and are bit-identical.
- **Independent review** (reviews/g2_review): no leakage. The solver is exact, and all P1/Δ values reproduce with
  independent code.
- **Non-essential genes** get P1 0.11. This is not predicted line offsets: it is broad library- and lineage-linked
  structure.
- **Within-library P1** of the best model is 0.279 (Δ vs B1 +0.108). It will be reported next to the headline (D16).

Next (G3):
- independent code/result review of G2;
- M1b weight tuning;
- line-offset diagnostic;
- M2 pass-through fix;
- M3 LightGBM and M4 neural go/no-go;
- M5 stacking;
- leave-one-lineage-out and no-EC variants.

Also pending:
- provenance fixes: registry hashes, `freeze-eval` tag, B0/B1/B2 re-run on a clean tree;
- disk housekeeping: cache/models is 15 GB and C: is 97% full.

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

- D14 (post-G2): M1's residual target switches from the shrunk lineage mean (B1) to the gene mean (B0).
  - Dev-CV evidence: Δ +0.0084 (SE 0.0017), 5/5 folds, which clears δ = 0.005 and 1 SE. The test set is still
    blinded.
  - This is a model-selection choice, not a change to metrics, labels or controls. The B1-residual run stays in the
    registry and is reported.
  - Consequence: the E-only M1 (B0 residual, grouped λ) equals R1 (raw ridge, single λ) within noise (Δ +0.0008,
    n.s.). It also equals T (Δ +0.0001).
  - So the current best is effectively ridge regression on expression. The honest framing: our ladder has not yet
    beaten the published linear baseline.
  - All later M1b/M2/M3/M5 runs use the B0 residual.
- D15 (post-G2): the default 1/√p_b block weights make every add-one modality look harmful. Leave-one-out shows C
  (+0.009) and E (+0.021) carry signal.
  - "Genomics adds nothing" may not be claimed until the pre-registered block-weight tuning (DESIGN §5: training
    LOO, coordinate descent) has been run.
  - The provenance arm (+0.014) is a batch proxy and is excluded, as pre-registered. It is never a candidate.
- D16 (post-G2, revised after the independent review reviews/g2_review/): the line-offset hypothesis is refuted.
  - Per-line centring keeps 99–100% of P1, and an offset-only model scores ≤ 0.024.
  - Non-essential-gene skill is broad. It is partly library group and partly lineage, and it scales with every
    model's general skill: the NEG/variable P1 ratio is about 0.35 for B1, B2 and M1.
  - Diagnostic adopted instead: report within-library and within-lineage P1 next to the headline. DESIGN §4 already
    promises a by-library column.
  - Dev numbers:
    - within-library: M1-E-resB0 0.279, ΔP1 vs B1 +0.108;
    - within-lineage: M1-E 0.237 vs B2 0.177, so skill beyond lineage is real.
- D17 (post-review): the current best model is M1-E with the B0 residual (`m1_e_resb0`). It is reported with its
  run ID from the clean re-run on the `freeze-eval` tree.
  - R1 and T are equivalent within noise. R1 is the published comparator, and T is the transfer model.
  - The review's added controls are adopted:
    - expression permuted within lineage, and within lineage × library;
    - a permuted B0-path run.
  - A Nadeau–Bengio corrected SE is shown alongside the fold SE. This is disclosure only; the selection rule is
    unchanged.
- D18 (pre-G3): M2 did not clear δ, so M3 (LightGBM) is fitted on the residuals of the current best rung.
  - The current best is M1-E-resB0, or M2 if the fixed M2 (D14 base, pass-through) clears δ first.
  - This follows the ladder's intent (each rung corrects the previous one) without inheriting a failed rung.
  - G3 runs in two parallel tracks with separate files:
    - kernel track: M1b weight tuning, M2 fix, report columns and controls, secondary splits, no-EC variant;
    - nonlinear track: M4 pre-check and M4, M3, R2/R4 and the published-comparable block.
  - Runtimes during G3 are indicative only. Finalists are re-timed alone at G4 (lesson 14).
- D19 (G3 nonlinear track; docs/g3_nonlinear_results.md):
  - M4 is dropped at its pre-registered pre-check: an MLP on PCs scored 0.3028 against ridge on the same PCs at
    0.3068, losing on all 5 seeds.
  - M3 is no-go. The version of record has nested early stopping (deviation D-G3-1, decided before any result):
    ΔP1 +0.0018 (SE 0.0002) < δ.
  - The un-nested sensitivity run (+0.0043) is not used anywhere. Its training-row features came from models that
    saw the validation targets. If M5 needs an M3 ingredient, nested M3 is the only admissible one; choosing the
    other would be post hoc.
  - R2 deviation D-G3-3 is accepted. The pre-declared 500-gene fallback projected to about 25 h, so the elastic net
    ran on the first 150 genes of the same pre-declared random subset. R2 is a comparator, not a candidate, and the
    150 genes were fixed by the earlier random draw, not by results.
  - R4's top-k choices pile up at the grid edge (k = 200). A k extension is queued so the prior-method comparison is
    fair; until then R4's gap is reported as an upper bound.
  - M3 peaked at 4.2 GB, above the 4 GB guide in §10. Logged, not a blocker.
- D20 (G3 kernel track; docs/g3_kernel_results.md):
  - **(a) Selection:** no M1b or M2 candidate clears δ. The best is M2 + paralogs at +0.0034. Under tuned weights,
    genomics (C, M, H, F, Meta) carries real but tiny signal: at most +0.0015, 5/5, Holm p < 0.01. D15's condition
    is met, so the finding stands as "genomics adds no selectable gain over expression". It is not "no signal".
  - **(b) Control misses accepted as interpreted.** Both are in the conservative direction or informative:
    - within-lineage permutation gives P1 0.127, below the expected band;
    - lineage × library permutation leaves 0.071 within-lineage P1, so library groups explain within-lineage
      variation.
    - The lineage × library P1 column (headline 0.199) is reported next to the headline as the strictest
      decomposition.
  - **(c) λ criterion:** a declared diagnostic is added for M5/G3 close-out. The headline is re-fitted with
    leave-one-out MSE and with grouped inner-CV r instead of leave-one-out r. If one beats the headline by > δ,
    it is treated like D14 (dev-CV model selection; test still blind). Otherwise it is reported as robustness only.
  - **(d) Erratum, EC lines:** DESIGN §3/§7.5 and configs/transfer_criteria.json say "7 germ-cell/EC dev lines".
    Only 4 are dev lines (NCCIT, NTERA-2 cl.D1, PA-1, SUSA): GCT27 is test, and 1411H and SCH are unscreened.
    - The frozen criteria file is not edited (tag `freeze-transfer`). This erratum fixes the wording only.
    - The H3 no-EC variant is defined operationally as "trained without the EC dev lines", which is 4.
    - No threshold or label changes.
  - **(e) Prediction recorded before any set-B access:** T will probably fail H3 (hPSC-specific dependencies),
    because POU5F1 and PRDM14 are not predicted in EC lines even with them in training (−0.10 vs −0.88 observed).
    This is stated now so a failed H3 cannot be explained away later.
  - **(f) Caches:** keep cache/m2_inner and cache/runs until M5 and G4 finish; delete after G4.
- D21 (pre-G4, held-out evaluation scripts; scripts/evaluate_test.py, evaluate_transfer.py; proposed, PI to
  confirm). All choices were made on toy data and a dev-only dry run (dev fold 0 standing in for the test set),
  before any test, set-A or set-B access.
  - **(a) Access-log protocol:** results/test_access_log.csv is tracked, and the guard appends a row on each access.
    So the scripts commit that file alone before every further access and at the end. Otherwise the second
    finalist would be refused with "dirty tree". The alternative is to untrack the log.
  - **(b) P2b CI method:** resampled per-gene AP is biased upward, because a duplicated positive line credits its own
    copy. On the dev dry run, nAP was +0.017 and the percentile CI excluded the point. The P2b CI is therefore
    point ± 1.96 × bootstrap SE, with the percentile interval and the bias also reported. P1, P2a and AUROC are not
    affected. The frozen transfer bootstrap has the same, much smaller effect (estimated at a few thousandths) on Q3a AP. It is
    disclosed as boot_mean_minus_point in the G6 output, and the frozen code is unchanged.
  - **(c) Post-freeze change acknowledged:** scripts/transfer_criteria.py differs from its freeze-transfer hash. In
    commit 8a613f0 only the read path changed (guard.load_harmonised); computations are unchanged (independent
    review). G6 accepts exactly the pinned content (sha256 87ecfaf5…, LF-normalised).
  - **(d) Unguarded read, logged:** while restricting the anchor readers, GSE156284 (Mair 2019 counts, which include
    the hESC set-B arms next to the HAP1 anchor) was loaded once to list its column names. No values were inspected
    or used. The G6 anchor reader now parses only anchor columns.
  - **(e) Dry-run observation:** with the real anchor screens, no library-matched anchor has positive T skill (all
    lower CIs ≤ 0 on dev-fold stand-in models). If this holds at G6, the anchor gate is not informative and
    "fails" is unreachable. Verdicts would then be transfers or inconclusive only.

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
  - [x] Provenance and freeze (after G2/G5a):
    - [x] registry columns eval_genes_sha256, gene_classes_sha256, tree_clean (old rows empty);
    - [x] guarded reader `guard.load_harmonised` for cache/harmonised (direct readers routed; static test);
    - [x] disk: models persist only with `--save-models` / `save_models` (default off); cache/models and
      regenerable cache/runs cleared (about 19 GB);
    - [x] `freeze-eval` tag at 62aab5b (guard.py hash refreshed for load_harmonised);
    - [x] Re-run on frozen evaluation code, 6 threads, tree_clean=1. Every number is identical to the earlier rows.
      - B0 20261026-220531; B1 -220921 (0.18374); B2 -221456 (0.25267);
      - m1_E_resB0 -222117 (0.3081; headline candidate, D17); t_E_status -224125 (0.3082); m1_E -225039
        (0.2997);
      - r1_E -223037 (0.3073) registered tree_clean=0, most likely because the new test file was briefly untracked.
        Its numbers are identical.
    - Note: transfer_criteria.json is committed, so the guard now refuses full-profile re-runs of the
      replicate-ceiling/transfer-criteria scripts. This is intended: they are frozen.
- [x] G2 M1/M1b kernel ridge, R1, transfer model T, modality ablation, M2 local stage (docs/g2_results.md)
  - [x] M1-E (B1 residual) 0.2997; R1 0.3073; M1-E B0 residual 0.3081; T 0.3082. All 5/5 folds above B2 (+0.047
    to +0.055). Runtime about 400–490 s for 5 folds.
  - [-] Add-one modalities under default weights: all ≤ 0 (C −0.013, H −0.024, F −0.019). [!] Leave-one-out: C +0.009
    and E +0.021 (Holm p < 0.001). Weight tuning is needed (D15).
  - [-] M2 local stage: paralogs +0.0048 (< δ), and it lowers P1-all (0.150 → 0.143). Not selected; pass-through
    fix in G3.
  - [-] P arm (253 lines): −0.006, underpowered.
  - [x] Negative control (permuted E) = B1 exactly. Provenance arm +0.014 = batch proxy, excluded.
  - [x] Independent review (reviews/g2_review, 8cf9d42): no leakage; solver exact; numbers reproduce. Provenance:
    every G2 run had a dirty tree, so finalists need a re-run on the clean tagged tree.
  - [!] Non-essential genes P1 0.110 (B1 0.068). Not line offsets (centring keeps 99%). It is library- and
    lineage-linked; within-library P1 0.279 (D16).
- [ ] G3 model ladder continued
  - [x] Re-run m1_E_resB0, r1_E, t_E_status and m1_E on frozen evaluation code (run IDs above, under G0)
  - [x] Report columns (scripts/g3_kernel_report.py; docs/g3_kernel_results.md): within-library and
    within-lineage P1 (reviewer definitions; reproduce the review exactly), within lineage × library P1, NB SE
    (≈1.5× fold SE). m1_E_resB0: 0.279 / 0.245 / 0.199.
  - [x] Stronger negative controls. Permuted B0 path: P1 0.003. Within lineage: within-lineage P1 −0.004.
    Lineage × library: within-lin×lib −0.007, but within-lineage 0.071 (> expected 0.03), so library groups carry
    within-lineage structure. [!] On the B0 path the LOO-r criterion → −1 as λ → ∞ (D3 artefact), so with no
    signal it picks the smallest λ. Possible small-λ bias on weak genes (decision pending).
  - [x] Solver tests moved into tests/ (p>n, LOO with B1 offset, real-path held-out feature perturbation), plus
    tests for M1b tuning, M2 nesting, exclude_from_training/LOLO.
  - [-] M1b block-weight tuning (B0 residual; ≤25 eigh): all tuned sets beat E-only, but only by +0.0004 (E+C)
    to +0.0015 (full), 5/5, Holm p < 0.01. Below δ. Tuned θ: C 1/32–1/64, M 1/16, H/F ≤ 1/64, Meta 1/4.
  - [-] M2 on the D14 base (pass-through, nested local λ): none +0.0019, paralogs +0.0034, partners +0.0025.
    Below δ. P1-all is no longer hurt. Paralogs add +0.010 within lineage × library cells.
  - [x] Leave-one-lineage-out (18 lineages ≥ 20 lines): m1_E_resB0 0.203, T 0.205, B2 0.120, B1 0. M1 > B2
    in 18/18. EC leave-out: no effect on P1. POU5F1/PRDM14 are not predicted even with EC lines in training
    (−0.09 vs observed −0.88). Only 4 of the "7 EC lines" are dev lines.
  - [x] H3 no-EC variant: t_E_status_noEC (20261027-000415), dev CV P1 0.3081; to be frozen at G4.
  - [-] M4: the pre-check MLP lost to ridge on PCs (0.3028 vs 0.3068, 5/5 seeds). Dropped (D19).
  - [-] M3 LightGBM on m1_E_resB0 residuals, nested: ΔP1 +0.0018 ± 0.0002 (< δ). No-go. Gains are threshold genes
    (NRAS hotspot; X-buffered RPS4X, EIF1AX; NAA10). Permuted control −0.004.
  - [x] R2 elastic net (150 genes): −0.030 vs ridge on the same genes, at about 45 s per gene.
  - [x] R4 TrPLet-like (500 genes): nested −0.038; non-nested optimism +0.018.
  - [x] R4 top-k extension (k ≤ 2000): ΔP1 −0.0196 vs ridge on the same 500 genes, with choices again at the edge.
    R4 improves monotonically toward ridge as it keeps more features, so feature selection buys nothing over
    shrinkage here. Not extended further.
  - [x] Published-comparable block (loosely comparable): m1_E_resB0 mean Spearman 0.284 on DeepDEP genes and
    0.303 on LEAP genes; 4,362 genes with pooled r ≥ 0.2.
  - [ ] λ-criterion diagnostic on the headline: LOO-MSE and grouped inner-CV r (D20c)
  - [ ] M5 stacking (inputs: m1_E_resB0, m2n_E_paralogs, g3_m3n nested); plain averaging is the baseline
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
