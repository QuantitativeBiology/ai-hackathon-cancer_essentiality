# Round 1 — PI response (DESIGN v0 → v1)

Reviews: `round1_{external_pi,statistics,prior_work,data,algorithms_compute,agent_failure_modes,biology}.md`.
Legend: **A** accepted · **M** accepted with modification · **R** rejected (reason given).

## External PI (major revision)
1. Core-essential test scored against observed screen labels, paired with B0; lost/gained essentials — **A** (§7.4).
2. RPE1/non-transformed lines removed from anchors; anchors predicted out-of-fold — **A** (§7.2). Anchors are dev
   lines, so they get CV out-of-fold predictions.
3. Pre-registered transfers/fails/inconclusive rule — **A** (§7.6, `configs/transfer_criteria.json`).
4. Expression-swap negative control — **A** (§7.5, §7.6).
5. Transfer model T — **A** (§5).
6. Ordered work plan with gates; acquisition early; time boxes — **A** (§12). Acquisition was already running.
7. Set A tiers; RPE1 sublines counted as one background — **A** (§3).
8. Within-lineage r; ΔP1 vs B1 as headline — **A** (§4).
9. P1 computed per fold; paired 1-SE — **A** (§4).
10. P2a as Δ vs B0 — **A**.
11. Rank-only metrics for BF screens; CRISPRi descriptive — **A** (§7.1).
12. Applicability-domain distance — **A** (§7.3).
13. User download list now; Nusinow proteomics — **A** (§2; `docs/USER_DOWNLOAD_REQUEST.md`).

## Statistics (accept with revisions)
1. P1 within fold, then averaged — **A**. This is where it conflicted with prior-work #3, which asked for
   concatenated OOF predictions. The statistician is right: pooling fold-specific means gives B0 a spurious negative r.
   The published-comparable block may still use pooled values, flagged as such.
3. Fold-wise gene classes; frozen dev list committed before unblinding — **A**.
4. Normalised AP + AUROC; prevalence ≤ 0.5 headline; tie handling; constant-predictor unit test — **A**.
6. Paired 1-SE, complexity order, δ = 0.005 — **A**.
7. One named headline model — **A**.
11. Set A: line as unit, ranked in the test distribution; no per-gene r across RPE1 — **A**.
12. Operational definition of "context-specific survive" — **A** (§7.5 headline recall at fixed NEG FPR + nulls).
13. Affine-scaling validation by medians/SSMD, not AUROC — **A**.
14. Anchor agreement on deviations; evaluated cell type excluded; dev-only assertion — **A**.
2, 5, 8, 9, 10, 15, 16 (should) — **A**. Power statement added in short form.

## Prior work (revise)
1. Rung R re-runs (ridge, EN, RF, TrPLet-like) — **M**. RF is limited to a pre-declared random 500 variable genes,
   and EN to variable genes (compute).
3. Concatenated OOF P1 — **R** for selection (see statistics #1); kept in the published-comparable block.
4. Published-comparable block — **A**. The gene lists are used if obtainable from the papers' supplements.
5. Novelty claims corrected (Shi 2024 lineage-confounder null; EFFECT partial correlations) — **A**. LITERATURE.md
   to be amended; within-lineage r restored.
8. Extra local features (not-expressed hinge, deep deletion, paralog mutation, family size) — **A**.
2. RRR / matrix completion — **M**: diagnostic only (the algorithms review found low-rank Y hurt variable genes).
6. Library-stratified reporting — **M**: the NaN line group serves as library proxy; CRISPRScreenMap is not local.
   Lines screened in two libraries cannot be identified without ScreenMap, so this is requested from the user.
7. Expression transfer variants; E-only transfer model; kNN — **A** (§7.3, §5).
9. Proteomics arm — **A** (Nusinow).

## Data (revise)
1. No mean-fill; masked/pattern-wise fits; ≥50 observed held-out lines; report by NaN group — **A**.
2. GrowthPattern/ModelType/PrimaryOrMetastasis ablation-only — **A**.
3. E-only transfer model — **A**.
4. CN transform (+0.1, clip); no per-gene indicator; 1206 count; autosomes — **A**.
5. Fusion filter (self-fusions, non-coding partners, ≥5 lines) — **A**.
6. TPM on all genes, then subset; QuantSeq flagged — **A**.
7. Set A per background; verify Hs 729 — **A**.
8–10 (should) — **A**. Anchor shrinkage is explicit, and the correction is applied only if leave-one-anchor-out
   validation improves.

## Algorithms & compute (sound)
1. eigh + LOO, wide grid, grouped λ — **A**.
2. Leave-line-out lineage residuals; κ in grouped CV — **A**.
3. M2 on OOF M1 — **A**.
4. Factorised M4 with pre-check — **A**.
5. Nested stacking; averaging as the baseline — **A**.
6. NaN handling — **M**: pattern-wise exact fits instead of dropping >10%-NaN genes, so those genes keep
   predictions (data review #1).
7. Threading — **A**.
8. Numerics — **A**.
Should items — **A** (M1b, RRR diagnostic, M3 go/no-go, runtime table, λ re-tuned per ablation).

## Agent failure modes (good intent, needs code-level enforcement)
1–7 — **A** (§9).
8. Probe marked flawed; its bugs become tests — **A** (RESEARCH_TREE 1.3 marked [-]).
9–13 — **A**. The status block already sits at the top of RESEARCH_TREE.

## Biology (revise)
1. Rare-selective class — **A**.
2. Set A tiers by mechanism — **A**.
3. TP53 hypothesis split (learnable vs diagnostic) — **A**.
4. Isogenic TP53 counterfactual — **A** (H2).
5. Profile-swap null + nearest-lineage comparator — **A**.
6. Recall at fixed NEG FPR as headline — **A**.
7. EC leave-out for pluripotency — **A** (H3).
8. Post-mitotic divergence; fitness-only rule — **A**.
9. CRISPRi KO↔CRISPRi anchors — **M**: used if the data are obtainable; otherwise CRISPRi stays purely descriptive.
10. ENCODE-matched expression correction; rank variant; state matching — **A**.
11. Short-screen QC gating — **A**.
12–17 (should/could) — **A**. Cas9 delivery mode is recorded in the inventory (could).

## Not adopted from any review
None of the must items was rejected outright. The only conflict, P1 aggregation, was resolved in favour of the
statistics review.
