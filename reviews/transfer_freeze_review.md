# Independent review — transfer evaluation, pre-freeze (`freeze-transfer`)

Reviewer: independent (statistics + CRISPR screens). I did not build any of this.
Scope reviewed: DESIGN.md v3 §7 (and §4, §5 where §7 depends on them), RESEARCH_TREE.md D7–D12,
`configs/transfer_criteria.json`, `configs/harmonisation_frozen.json`, `configs/screen_plan.csv`,
`configs/control_split.json`, `configs/expression_harmonisation.json`, `docs/harmonisation_validation.md`,
`results/harmonisation/*` (summary tables only), `src/essml/harmonise.py`,
`scripts/transfer_criteria.py`, `scripts/harmonise_screens.py`.

## Verdict: **FIX FIRST** — do not tag `freeze-transfer` yet.

The harmonisation work is careful, well controlled and honestly documented: the blinding story
(controls-only reads, write-only cache, aggregate ceilings), the simulated-batch recovery
(corr(δ̂, δ) = 0.965), the wrong-line and within-set-permutation negative controls on step (iii), the
pre-stated and *honoured* Chronos adoption margin, and the per-screen reachability plan are all above
the standard I normally see. The problem is not the harmonisation; it is the **decision rule**. As
written, `configs/transfer_criteria.json` is not a pre-registration: it leaves at least nine choices
open that can be made after model predictions are in hand, its MDE is not the quantity DESIGN §7.6
specifies, its uncertainty model omits the dominant variance component, and its two headline gates are
asymmetric in a direction that makes "transfers" markedly easier to reach than "fails". Freezing it
now would freeze a rule whose outcome is still partly at the analyst's discretion.

Fix F1–F11 before the tag. F12–F18 should be fixed too, but could be recorded as declared deviations
if the PI prefers.

---

## A. Must-fix: the MDE and the uncertainty model

### F1. The MDE is not the quantity DESIGN §7.6 defines (severity: **critical**)

DESIGN §7.6: "MDE: 2 × the bootstrap SE of the replicate-ceiling **skill**", and "The minimum
detectable effect comes from the replicate ceiling **and the anchor noise**."

`scripts/transfer_criteria.py` computes `MDE = 2 × SE(auprc_all)` — twice the SE of the ceiling's
**absolute** AUPRC, not of the ceiling's *skill* (ceiling AUPRC − B0 AUPRC). B0 never enters the MDE
computation, and anchor noise never enters it at all. The verdict statistic, however, *is* a skill
(ΔAUPRC vs B0). So the MDE and the statistic it gates are on different footings: the MDE is the
sampling noise of a level, the statistic is a paired difference. There is no reason the two should be
comparable, and the paired difference will generally have a much smaller SE than the level.

**Fix.** Compute the ceiling's skill per directed pair — `AUPRC(ceiling predictor) − AUPRC(B0)` on the
same gene universe and label as the verdict statistic — and set `MDE = 2 × SE` of *that*, by the same
gene-block bootstrap. Add the anchor-noise term §7.6 requires (e.g. add in quadrature the bootstrap SE
of the same-pipeline anchor skill), or strike "and the anchor noise" from §7.6 as a recorded
deviation. Re-derive all eight MDE values.

### F2. The CI and the MDE ignore between-screen/between-study variance, which dominates (severity: **critical**)

The gene-block bootstrap resamples genes only. Both the ceiling SE and (per the file) the verdict CI
therefore describe gene sampling and nothing else. The numbers show how badly that understates the
real spread:

| background | ceiling (all genes) | gene-bootstrap SE | observed spread across study pairs |
|---|---|---|---|
| RPE1 | 0.661 | 0.0078 | `transfer_ceiling_pairs.csv` auprc_all 0.47–0.92; CEG/NEG-free 0.40–0.77 |
| hPSC | 0.516 | 0.0067 | 0.48–0.62 |

The between-study SD of the RPE1 ceiling is ≈0.10 — **13×** the gene-sampling SE used to set the MDE.
`docs/harmonisation_validation.md` §10 states this caveat and then concludes "'fails' is a strict but
reachable condition". It is reachable; that is the problem. "Fails" would be declared when T adds less
than 0.016 AUPRC over the gene mean — a margin six times smaller than the study-to-study variation of
the very benchmark that defines it. The same omission makes the "transfers" gate (lower CI > 0)
near-automatic: with ~19k genes and 11.9k blocks, any positive point estimate will clear zero.

This is compounded by a structural fact I did not see acknowledged anywhere: **every RPE1 screen uses
the same matched expression input** (`configs/screen_plan.csv`: "DepMap 26Q3 RPE1-ss sublines" for all
18 arms). Model T's input for the 18 RPE1 arms is therefore one expression vector plus a TP53 flag, so
T's *prediction* is effectively one of two vectors. All variation across the 18 "evaluable screens"
comes from the observed labels. Averaging over them and bootstrapping genes treats one prediction
measured against 18 correlated label sets as 18 independent screens.

**Fix (all three parts).**
1. Make the resampling hierarchical: resample **studies** (or arms within study) *and* gene blocks, so
   the CI carries the between-study component. With 7 RPE1 and 3 hPSC studies, prefer a
   leave-one-study-out / studentised interval over a percentile bootstrap, and say in the file which
   one is used.
2. Define the background skill as a **mean over study representatives**, not over arms (see F3).
3. State explicitly in the file that the prediction vector is shared across a background's screens, and
   that the CI is therefore a statement about label/gene sampling plus study choice, not about
   independent replication of the prediction.

### F3. The verdict population ≠ the ceiling/MDE population (severity: **high**)

The ceiling and MDE use **one representative per study, focused libraries excluded**: 7 RPE1 studies,
3 hPSC studies (`ceiling_study_representatives`). The verdict statistic uses
`backgrounds.RPE1.evaluable_screens` — **18 arms**, including:

- 9 Olivieri NT arms (6 TKOv3 + 3 TKOv2) — 50% of the mean, from one study, and these arms are
  near-replicates of each other (untreated arms of a drug-screen series);
- `brown2019_WT_PACE002_T18` **and** `_T30` — two timepoints of the same arm, counted twice;
- `bowden2020_WT_D19` / `_KO_D19` — a focused dual-gRNA DDR library, excluded from the ceiling as
  `focused`, but included in the skill mean on a gene universe orders of magnitude smaller.

hPSC is worse: `evaluable_screens` lists 4 screens, but `screen_level_notes.hPSC` says only Mair can
reach transfers/fails (Yilmaz: D12, no same-pipeline anchor; Ihry: no anchor). So the hPSC *background*
verdict would be computed from a mean that is 50–75% composed of screens the plan declares
verdict-incapable. Either the "inconclusive only" designation means something or it does not.

**Fix.** Define one population and use it for ceiling, MDE and verdict: one representative arm per
study (the existing most-negative-SSMD rule is fine and is control-based), focused/sub-genome libraries
excluded from the background mean, and **only verdict-capable screens** (same-pipeline anchor + ceiling)
contributing. Declare the resulting screen list explicitly in the file — it should be the same list the
MDE was computed on. State a minimum gene universe per contributing screen.

### F4. Gene blocks are narrower than DESIGN §7.6 specifies (severity: **medium-high**)

§7.6: "a gene block bootstrap (blocks = gene **families/complexes**)". The implementation uses paralog
families only: 11,934 blocks over 19,179 genes (≈62% singletons), max block 162. Complex membership —
the dominant source of correlated essentiality error (ribosome, proteasome, spliceosome, mediator) — is
not blocked at all, so the CIs are too narrow in exactly the regime that matters for core-essential
AUPRC.

**Fix.** Add complex co-membership (CORUM / the priors already in `cache/priors/`) to the block
definition, union with paralog components; or record the paralog-only choice as a declared deviation
from §7.6 with a stated justification. Re-derive the SEs and MDEs either way.

---

## B. Must-fix: loopholes — choices still open after predictions are seen

Each item below is a free parameter, an ambiguous definition or a missing tie-break. I list them as the
brief asks; any one of them permits a different verdict from the same predictions.

### F5. The ≥50%-of-anchor-skill gate is undefined *and* vacuous as applied (severity: **critical**)

Three distinct problems in one clause ("skill ≥ 50% of the same-pipeline anchor skill (same lab/library
anchors through the same harmonisation, D7 OOF)"):

**(a) The anchor-skill statistic is never defined.** Nowhere in `transfer_criteria.json`, DESIGN §7.6 or
`docs/harmonisation_validation.md` is it said what "anchor skill" is computed *on*: which label (the
anchor's own harmonised external effect < −0.5, or its DepMap effect?), which gene universe, which
label variant, and how several anchors are aggregated (mean? median? per-anchor minimum?). For Q3b, it
is not even stated whether the anchor recall uses the same 5%-FPR-on-NEG-evaluation construction.
This is the single largest loophole in the file: the denominator of the headline gate is at the
analyst's discretion after the fact.

**(b) "Same lab/library" is ambiguous, and the code resolves it to *lab*.**
`scripts/harmonise_screens.py` builds `lib_anchor` keyed on `e["group"]` (the correction group = lab),
not library. So `configs/screen_plan.csv` lists, for every TKOv3 RPE1 screen, the anchors
`hart2015_DLD1` + `hart2015_HCT116` (**TKOv1**, published guide-LFC) + `mair2019_HAP1_T18` (TKOv3,
counts). Library-level reading gives one anchor (HAP1); lab-level gives three, two of them from a
different library *and* a different data level. The two readings give very different denominators.
Note also that this is inconsistent with D12's own principle: D12 disqualified the Wang anchors for
Yilmaz because the anchors are published gene scores while the target is count-scored here — the same
mismatch (published guide LFC vs counts) is accepted for the TKOv1 anchors of the TKOv3 screens.

**(c) With the actual anchors, the gate is vacuous or undefined.** The only TKOv3 anchor is HAP1, whose
same-line deviation agreement with DepMap is ~0.09 (0.061 on its own, group mean 0.092–0.093; NN rank
34/1022; `docs/harmonisation_validation.md` §3 — "HAP1 agreement is poor at every timepoint"). Model T's
ΔAUPRC skill on HAP1 is therefore likely to be near zero and may well be negative. The rule then reads:
- anchor skill ≈ 0 → any skill > 0 clears "≥50% of ≈0": the gate does no work;
- anchor skill < 0 → "≥ 50% of a negative number" is satisfied by *negative* model skill: the gate
  inverts;
- and there is no rule at all for what happens when anchor skill is ≤ 0 or its CI covers 0.

Meanwhile the **"fails"** rule *does* require "same-pipeline anchors show positive skill (lower CI > 0)".
So a weak anchor blocks "fails" while leaving "transfers" wide open. That asymmetry is the opposite of
what a conservative pre-registration should do, and with HAP1 as the only TKOv3/hPSC anchor it is the
likely operating regime.

**Fix.** In the frozen file: (i) write the anchor-skill statistic out in full — label, gene universe,
label variant, aggregation across anchors, and the D7-matched model (single OOF model per anchor);
(ii) resolve "same lab/library" to **library *and* data level**, and list the anchor ids per screen
explicitly (for TKOv3/hPSC that is HAP1 alone — say so); (iii) require anchor skill to have a lower CI
> 0 for the ratio to be evaluated at all, otherwise the background is **inconclusive** regardless of
model skill; (iv) add an absolute floor so the gate cannot be cleared by a near-zero anchor (e.g. skill
≥ max(0.5 × anchor skill, a pre-stated absolute ΔAUPRC)); (v) make the anchor condition symmetric
between "transfers" and "fails".

### F6. No primary label / no primary MDE is declared (severity: **high**)

`labels` defines two labels (effect < −0.5, and the scale-free top-10% rank label). The MDE block gives
**three** Q3a values (`Q3a_dAUPRC`, `..._ctrlfree`, `..._rank10`) plus `Q3b_dRecall`. `Q3a.verdict_statistic`
names none of them, and `Q3a.primary_reported` is a *ratio* to the ceiling, which is not the verdict
statistic at all. DESIGN §7.4 lists all three as reported. So after seeing predictions, the analyst
picks which of three labels decides the verdict, with three different MDEs available.

Compounding this: Q3a per screen × 2 questions × 2 backgrounds × 3 labels × (T and the full model) is a
large multiplicity with no adjustment, while DESIGN §4 requires Holm adjustment for the pre-declared
ablation contrasts.

**Fix.** Name exactly one primary label and one primary MDE per question in the file (I recommend the
**CEG/NEG-free** variant as primary for Q3a — see F7), demote the other two to secondary with an
explicit "not verdict-bearing" note, and state that no verdict is taken on a secondary label. State
whether the verdict uses T only (DESIGN §7.4 says T is the primary evidence; keep the full model
secondary).

### F7. The headline ceiling and the Q3a verdict include genes whose values the affine step pinned (severity: **high**)

`auprc_all` (the ceiling reported in the file, 0.6613 / 0.5161, and the source of the headline MDE)
includes the CEGv2/NEGv1 **scaling** halves. Those are precisely the genes the affine map forces to
`-1.4677` and `+0.0072`. Their label (effect < −0.5) is therefore fixed by harmonisation, not observed,
and they are predicted well by any method including B0. This inflates the absolute AUPRC of the
ceiling, of T and of B0 alike, and so inflates the `AUPRC_T / ceiling` ratio that §7.4 and the file call
the primary reported evidence. The control-free variants exist (`auprc_ctrlfree`, 0.5436 / 0.4481) —
they should be the headline, not the alternative.

**Fix.** Make the control-free (or evaluation-half-only) variant primary for the ceiling, the MDE and
the Q3a verdict; report the all-gene variant as a secondary with a one-line note that the scaling
halves are pinned by construction.

### F8. The Q3b comparator B1 is undefined for set B (severity: **high**)

`Q3b.verdict_statistic` is "recall(T) − recall(B1)". For a non-transformed screen there is no lineage.
`src/essml/baselines.py` makes B1 fall back to the gene mean for an unseen lineage — i.e. **B1 ≡ B0**.
DESIGN §5 confirms the primary is no lineage ("A pre-committed lineage mapping is a secondary variant").
But the file also carries a *third* object, `nearest_lineage_comparator` ("B1 prediction for the closest
cancer lineage"), with no definition of "closest" (metric? harmonised expression space? k? which
lineage labels?). So "skill over B1" could mean skill over B0, skill over an unstated nearest-lineage
mean, or skill over the §5 pre-committed mapping — a choice worth a lot on Q3b, and one that
`configs/expression_harmonisation.json` (applicability-domain distances already computed: hESC/iPSC
0.34–0.45, nearest neighbours embryonal carcinoma) makes easy to exploit after the fact. DESIGN §5 also
warns that mapping hPSC to Germ Cell "would route B1 through the EC lines and confound H3".

**Fix.** State in the file: for set B, B1 ≡ B0 under the no-lineage primary, so the Q3b verdict
comparator **is** B0; the nearest-lineage comparator is secondary, with "closest" defined exactly
(metric, space, k) or deleted.

### F9. The expression-swap cap contradicts §7.5 and is underspecified (severity: **medium-high**)

DESIGN §7.5 states the profile-swap null as a *comparison*: "Matched must beat mismatched."
`transfer_criteria.json` converts it to an *absolute* test: "If the swapped-profile skill has a lower
95% CI bound > 0, the verdict … is capped at 'inconclusive'". These are different rules with different
outcomes. On core essentials the swapped profile will usually retain real skill over B0 (shared
essentiality is not profile-specific), so the absolute form plausibly caps everything at inconclusive —
or, under a different reading, nothing. Unspecified: which profiles are swapped, how many swaps, and
whether the swap statistic is the mean, median or maximum over swaps.

**Fix.** Restore the §7.5 comparison form (matched skill > swapped skill, with a stated margin and a
paired CI), enumerate the swap set per background, and fix the aggregation (I suggest the mean over all
declared mismatched profiles, with the maximum reported).

### F10. Tie-breaks and aggregation rules that are simply absent (severity: **medium-high**, each a loophole)

1. **Screen weighting** inside a background: "mean over the background's evaluable screens" — equal per
   arm, per study, or gene-count weighted? (F3 shows this changes the answer.)
2. **Gene universe** for the background mean: each screen has its own universe (TKOv1/v2/v3, Brunello,
   custom DDR, ~2.5k-gene Novartis chip). The mean mixes AUPRCs computed on different gene sets, and
   none of them is the 19,179/18,979-gene universe the MDE was computed on. No common universe or
   minimum-gene rule is declared.
3. **NaN handling**: what happens when a screen's statistic is NaN (no positives, too few NEG-evaluation
   genes)? `nanmedian` is used for the ceiling; the verdict path says nothing.
4. **Exact boundary conditions**: "skill ≤ 0" vs "> 0", "≥ 50%", "lower bound > 0" — all stated, good —
   but "upper 95% CI bound of skill < MDE" has no rule for equality, and no rule for the case where
   *both* the transfers and fails conditions are vacuously satisfiable because anchor skill is ~0.
5. **Overall verdict with 2 backgrounds**: "strict majority" of 2 means both must agree (the file says
   so). But if one background becomes non-evaluable at analysis time, "strict majority of 1" = that one
   background decides the whole question. State that explicitly or require ≥2 evaluable backgrounds.
6. **AP tie handling**: `ap()` in `scripts/transfer_criteria.py` and `auprc()` in `harmonise.py` break
   ties by array order (Entrez sort), unlike DESIGN §4's "sklearn's step-wise AP" with its constant
   predictor = π unit test. Rank-only ORCS screens and 2-guide median LFCs do produce ties.

**Fix.** Write each of 1–5 into the file as an explicit rule; for 6, use sklearn or average within tie
groups so §4's unit-test guarantee holds for the set-B path too.

### F11. The criteria name models that are not frozen (severity: **high**, process)

The file's statistics are defined in terms of "model T", "the full model", B0, B1, B2. But
`configs/finalists.json` **does not exist** (DESIGN §5: it freezes the headline model, model T and the
H3 no-EC variant at G4; `src/essml/guard.py` requires it for the test gate), and RESEARCH_TREE marks G2
(M1/M1b/R1/T) still in progress. So at `freeze-transfer` the criteria would reference models whose
feature sets and hyperparameters are not yet fixed, which is the one thing a pre-registration cannot
leave open: "skill of T" is not computable without further choices.

**Fix.** Either tag `freeze-transfer` only after `configs/finalists.json` is committed, with the exact
run ids / config hashes for T, the full model, the H3 no-EC variant and B0/B1/B2 recorded in
`transfer_criteria.json`; or split the tag into `freeze-transfer-criteria` (thresholds, which is what is
ready) and `freeze-transfer-models`, and say in the file that the former is void unless the latter
precedes any set-B prediction.

---

## C. Must-fix: D7 consistency, and H1/H2/H3

### F12. DESIGN §5 still contradicts D7 (severity: **medium-high**, spec consistency)

D7 is correct: with 5-fold CV each dev line is out-of-fold for exactly **one** model, so DESIGN §5's
"Anchors: the average of the **4 fold-models not trained on that anchor**" is impossible, and §5's
"Set A and set B: the average of the 5 dev fold-models" (an ensemble) is superseded by D7's
per-fold-model skill averaging. `transfer_criteria.json.decisions.D7` matches D7. But DESIGN v3 §5 — the
frozen spec — still carries the old text, and §7.6 is silent. A reader of the spec at freeze time gets
the wrong rule.

There is also a residual ambiguity in D7's own wording as transcribed into the file ("skill = mean over
the per-fold (5 dev CV) frozen models"): averaging the five models' *predictions* and computing skill
once is not the same as computing skill five times and averaging. D7's tree entry says the latter
("Set A/B skill is computed per fold model and averaged; the 5-model ensemble is reported as
secondary"), and that is the version matched to a single-OOF anchor — but the file's phrasing admits
both.

**Fix.** Amend DESIGN §5 "Frozen models for evaluation" (or add a §5 note recording D7 as superseding
it), and in `transfer_criteria.json` write: "skill is computed separately for each of the 5 fold models
and the five skill values are averaged; the 5-model prediction ensemble is secondary". Note also that
the CI must then account for the five fold models being trained on overlapping data — state that the
bootstrap is applied to the averaged-skill statistic.

### F13. H1 as written is a degenerate test, and it is confounded (severity: **high**)

Three problems:

1. **The permutation test has no sampling unit.** All RPE1 screens share one expression input
   (`screen_plan.csv`), so model T's predicted MDM2/PPM1D percentiles take **one** value for the WT
   screens and **one** for the null screens, differing only through each screen's library gene universe.
   "Screen-level permutation p < 0.05 across WT vs null screens" therefore permutes labels over points
   that are identical by construction: the p-value measures the library universes, not the hypothesis.
   With 4 WT and 5 null RPE1 arms the smallest attainable p is also driven entirely by how many arms
   each study happens to contribute (see F3).
2. **H1 is not independent of H2.** Since the WT/null contrast *is* the TP53 feature flip on the same
   expression input, H1 and H2 are two readings of a single quantity. They should not be reported as two
   pieces of evidence.
3. **Confounding.** The criterion pools "RPE1 WT screens, **hPSC**" as the WT group against RPE1-null.
   hPSC differs from RPE1 in cell type, expression input and applicability-domain distance
   (0.34–0.45 vs 0.18), so a positive result is attributable to cell type rather than TP53 status.

Also worth stating plainly: because MDM2/PPM1D-vs-TP53 is one of the strongest signals in the DepMap
training data, H1 is close to the "near-tautological" negative controls §7.5 already declares. It is a
sanity check on the TP53 feature, not evidence of transfer.

**Fix.** Restrict H1 to RPE1 WT vs RPE1 null (identical expression input), drop the permutation
p-value in favour of an honest statement that this is a single two-point counterfactual contrast (or
bootstrap the gene universe if a CI is wanted), report hPSC separately and descriptively, label H1 a
feature sanity check rather than independent evidence, and say so in the file.

### F14. H2 has no minimum-n, an ambiguous r, and an un-cancelled scale (severity: **high**)

1. **No minimum-n.** The gene set is the 15 p53/DDR genes filtered by |observed Δ| > 0.25. How many
   survive is unknown in advance and could be 1–3. "Sign agreement ≥ 70%" over 2 genes means both; a
   Pearson r over 2 genes is ±1 by construction and over 3 is meaningless. The criterion can be met or
   missed by one gene.
2. **Ambiguous r.** "sign agreement >= 70% over genes with |observed delta| > 0.25 AND Pearson r > 0.3"
   does not say whether r is computed on the filtered subset or on all 15 genes.
3. **The delta is not batch-cancelling as §7.5 claims.** Bowden's WT and KO arms are affine-scaled
   **per arm**, each to its own control medians. p53-competent cutting toxicity shifts the WT arm's
   whole distribution — exactly the effect H2 is trying to measure — and the per-arm affine map absorbs
   part of it before the delta is taken. (p53/DDR genes are excluded from the gene-wise step, correctly,
   but *not* from the affine step.) Bowden is also a focused dual-gRNA DDR library with no anchor, so
   its control medians rest on few genes.

**Fix.** Declare a minimum number of qualifying genes (I suggest ≥6, else H2 is inconclusive by
declaration); state the gene set r is computed on; compute the observed WT−null delta **also** on the
pre-affine LFC scale (or on a common scaling applied to both arms jointly) and pre-register which is
primary; state that the Bowden library is focused and record how its affine parameters were obtained.

### F15. H3 tests the model against itself, with a shrinkage-sensitive absolute threshold (severity: **medium-high**)

"predicted effect < −0.5 for >= half of the set with the 7 germ-cell/EC lines in training" uses only
*predicted* values; the observed hPSC effects are "reported descriptively". H3 can therefore "pass"
while the hPSC screens show no such dependency, and fail while the model ranks the pluripotency genes
at the very top of the screen. The absolute −0.5 bar is also a property of the model's shrinkage, not of
transfer: a well-calibrated ridge predicting rare selective dependencies will rarely clear −0.5 for any
gene. Two further gaps: the 8-gene set mixes genes that are reliably essential in hPSC KO screens
(POU5F1, LIN28A) with genes that generally are not (NANOG, ZFP42, DPPA4, PRDM14), so "≥ half of 8" may
be unachievable for reasons unrelated to transfer; and the "narrow EC transfer" fallback can be
triggered by a one-gene difference (4/8 with EC vs 3/8 without) with no stability check.

**Fix.** Add an observed-side condition (the gene must be depleted in the hPSC screens for the
prediction to count), add a rank-based primary (e.g. predicted in the screen's most-depleted *x*%, with
*x* fixed now) alongside the absolute threshold, pre-state which of the 8 genes are eligible on
literature grounds only (not on the observed set-B values), and pre-state a margin for the narrow-EC
call (e.g. a difference of ≥2 genes).

### F16. "Fails" is unreachable-in-spirit for hPSC, and the doc's own reasoning is inverted (severity: **medium**)

`docs/harmonisation_validation.md` §7 states: "hPSC agreement is low, so 'fails' will be near-unreachable
for hPSC (large MDE)." The computed MDE says the opposite: hPSC `Q3a_dAUPRC` = **0.0134** < RPE1's
**0.0157**, and hPSC `Q3b_dRecall` = 0.0306 < RPE1's 0.0372. A *lower* ceiling gave a *tighter*
gene-sampling SE and therefore an *easier* "fails". This is the clearest symptom of F1: the MDE as
implemented is unrelated to how much skill is achievable in that background. For hPSC the achievable
headroom is tiny (ceiling recall 0.118; AUPRC 0.516 against a prevalence of ~0.31 — see
`replicate_ceilings_study_level.csv`), yet the MDE is the smallest of the four.

**Fix.** Fixing F1 (ceiling *skill*, plus anchor noise) largely resolves this; then re-check the §7 and
§10 narrative against the recomputed numbers and correct the claim. Separately: with 3 hPSC studies and
6 directed pairs, state honestly in the file that the hPSC MDE rests on 6 dependent pairs.

---

## D. Items (4) and (5) of the brief

### (4) Did any harmonisation choice use model predictions, or non-control genes of non-transformed screens beyond the allowed replicate-ceiling computation?

**No model prediction was used anywhere.** I checked this carefully and found the claim sound:
`scripts/harmonise_screens.py` routes set-B reads through `guard.load_external_screen(..., controls_only=True)`
for every scaling/QC decision; full profiles go write-only to `cache/harmonised/` with
`harmonise_write_only` rows; the ceilings read those files through `guard._log("replicate_ceiling_aggregate", ...)`
and return aggregates. No model is imported in the harmonisation or criteria scripts. DepMap targets are
dev-only (`guard.load_targets`); every anchor is a dev line. Good.

**But three choices did use non-control genes of non-transformed screens beyond the ceiling computation:**

**F17 (severity: high) — D12's justification.** D12 has a sound *a priori* rationale (the Wang anchors
are published gene-level CS, Yilmaz is count-scored here, so the pipelines differ). But both
`docs/harmonisation_validation.md` §7 and `harmonisation_frozen.json → genewise.per_screen` also justify
it by the *outcome*: "applying the correction lowered Yilmaz's agreement with Mair and Ihry (hPSC
0.13 → 0.08)". That agreement is computed over all genes of three non-transformed screens. DESIGN §7.2
requires harmonisation choices to be "fixed on controls and anchors only". Selecting a harmonisation
variant on a cross-set-B all-gene agreement statistic is outside that, and it is not the
replicate-ceiling computation — it is a choice *made on* the ceiling. It matters twice: it sets Yilmaz's
harmonisation *and* its "no same-pipeline anchor → inconclusive only" status, which in turn shapes the
hPSC background.

*Fix.* Restate D12 on the pipeline-mismatch rationale alone (which stands on its own and is checkable
from the anchors' data level), delete the agreement-based justification from the config and the doc, or
record it explicitly in RESEARCH_TREE as a deviation from §7.2's "controls and anchors only" with the
reasoning that the ceiling statistic is observed-vs-observed. Note which it is; do not leave both
justifications standing side by side.

**F18 (severity: medium) — the Chronos joint fit, and D9's stated reason.**
(a) The TKOv3 Chronos refit was a **joint** fit over the HAP1 anchor *and* the set-B count screens
(Olivieri ×6, Brown PACE002/006, Mair ×2), sharing guide efficacy across them. The adoption statistic
(HAP1 dev r, +0.013 vs the +0.02 margin) therefore depends on non-control genes of non-transformed
screens. Chronos was **not** adopted, so the frozen pipeline is unaffected — but the decision statistic
was contaminated, and `harmonisation_frozen.json` keeps the route switchable via
`chronos_tkov3_decision.json`. If it is ever switched on for a sensitivity analysis, the set-B profiles
and the anchor become mutually dependent.
(b) D9 (expression harmonisation) is a recorded override of the rule's outcome — properly flagged in
`pi_flag`, and the rule was written before the comparison, which is good practice. But DESIGN §7.3
restricts candidate selection to "matched ENCODE/GEO profiles of DepMap **dev** lines only", and the
override's stated reasons cite a set-B profile: "Yilmaz hESC has 0.71M [protein-coding TPM sum]". The
Tzelepis HT-1080 argument is within-population (an anchor/dev line); the Yilmaz one is not.

*Fix.* (a) State in `harmonisation_frozen.json` that the Chronos decision statistic came from a joint
fit including set-B screens, and that re-enabling the route requires re-deciding on an anchor-only (or
anchor-plus-transformed) fit. (b) Re-justify D9 on dev-line evidence only (the HT-1080 case and the MAE
numbers suffice), or record the Yilmaz citation as a deviation from §7.3's selection population.

Two further items I checked and found **clean**: the Chronos cutting negative controls are the guides
targeting the NEGv1 **scaling** half (not the evaluation half) — correct; and the study representative
for the ceiling is chosen by most-negative **evaluation-half SSMD**, which is a control-gene choice and
is legitimate, though note it is an evaluation-half statistic used for a selection (see F19 below).

### (5) Is the evaluation-half / scaling-half control split respected everywhere?

**Almost.** The mechanics are right and consistently implemented:

| use | set used | correct? |
|---|---|---|
| affine scaling (`affine_params`, driver lines 269/434/632/671) | `scaling_set` (589) + `neg_scaling` (462) | yes |
| D10 affine targets (DepMap dev medians) | `scaling_set` + `neg_scaling` | yes |
| QC (`qc_metrics`, every call site) | `ceg_evaluation` (342) + `neg_evaluation` (462) | yes |
| Q3b 5% FPR threshold | `neg_evaluation` | yes |
| ceiling `recall_q3b` negatives | `neg_evaluation` (`is_negev`) | yes |
| ceiling `auprc_ctrlfree` | excludes all CEG ∪ NEG | yes |
| `scaling_set` construction | dev common essentials minus `ceg_evaluation` minus NEGv1 | yes (the 258-gene overlap with `ceg_scaling` is by design) |

Two exceptions:

- **F7 above** is the substantive one: the *headline* ceiling and MDE (`auprc_all`) include the scaling
  halves whose values the affine step pinned. That is a leak of the scaling half into the primary
  statistic, even though the split itself was honoured.
- **F19 (severity: low-medium).** The ceiling's study representative is selected by most-negative
  **evaluation-half** SSMD. Using the held-out half to *select* the arm that then defines the ceiling
  and the MDE couples the evaluation half to the benchmark. It is a small effect and defensible (it is a
  control-gene criterion, not a model criterion), but the cleaner choice is to select on the scaling-half
  SSMD and keep the evaluation half purely for QC. *Fix:* switch the selection statistic to the scaling
  half, or record the choice as deliberate.

A note, not a finding: `qc_min_controls = 10` applied to the *halves* is what makes the Tuano HMEC
libraries unassessable (2 NEG evaluation genes). That is a legitimate consequence of the split and is
declared (D-c), not a violation.

---

## E. Smaller points worth fixing while the file is open

- `_status` in both `transfer_criteria.json` and `harmonisation_frozen.json` still says "DRAFT";
  flip both at the tag and record the commit, or the frozen artefact will claim it is not frozen.
- `transfer_criteria.json.decisions` records D7, D10, D12 but not **D8** (the rare-selective 1% floor,
  which feeds `gene_classes.csv` → `frac_below` → the Q3b screen-specific definition), **D9** (quantile
  expression harmonisation, which determines every set-B model input) or **D11**. Add them, and add the
  `configs/expression_harmonisation.json` path to `inputs` — it is a binding input and is currently
  unlisted.
- `inputs` should carry content hashes (the project already does this via `configs/frozen_hashes.json`
  for the eval path); otherwise "frozen" criteria depend on four mutable files.
- O5 in the validation doc notes commit `5181e51` holds mid-run versions of the harmonisation files and
  the final tree is uncommitted. The tag must be made on a tree where the committed files are the ones
  the numbers came from; the G0 "pending" list (registry hashes, `freeze-eval` tag, clean-tree baseline
  re-run) should close first, since `guard` already requires a clean tree.
- D-e ("the blind write path and the ceilings run through the driver … the frozen guard has no API for
  them; suggest adding one at `freeze-transfer`") is the right call — add the API, so the set-B write
  path is enforced by the guard rather than by convention.
- `qc_pass` is evaluated at the `1_affine` step while the final profile may carry the gene-wise step.
  Harmless today (no set-B screen receives a gene-wise correction after D12), but state it.

---

## F. What I would need to see to sign off a freeze

1. MDE recomputed as 2 × SE of the ceiling **skill** (ΔAUPRC/Δrecall vs B0), with the anchor-noise term
   (F1), on the control-free label (F7).
2. A hierarchical (study + gene-block) resampling scheme, with blocks widened to families **and**
   complexes (F2, F4).
3. One declared population — the same list for ceiling, MDE and verdict (F3).
4. The anchor-skill statistic written out in full, "same lab/library" resolved to library + data level,
   an anchor-skill lower-CI > 0 precondition, an absolute floor, and symmetry between transfers and
   fails (F5).
5. One primary label and one primary MDE per question; secondaries marked non-verdict-bearing (F6).
6. B1-for-set-B stated as B0; nearest-lineage defined or deleted (F8). Swap null restored to §7.5's
   comparison form with the swap set enumerated (F9).
7. Items 1–6 of F10 written into the file; AP tie handling aligned with §4 (F10.6).
8. `configs/finalists.json` committed with run ids/config hashes for T, the full model, the no-EC
   variant and B0/B1/B2 — or the tag split in two (F11).
9. DESIGN §5 amended for D7, and D7's "average of skills, not of predictions" made explicit (F12).
10. H1 restricted and relabelled; H2 given a minimum-n, an unambiguous r and a pre-affine delta; H3
    given an observed-side condition and a rank-based primary (F13–F15).
11. D12 justified on pipeline grounds only (or the deviation recorded); the Chronos joint-fit dependency
    and D9's set-B citation recorded (F17, F18).
12. The §7/§10 narrative corrected where it contradicts the computed MDEs (F16).

Items 1–8 are the ones that decide whether the rule is a pre-registration or a menu. I would not tag
until they are done.
