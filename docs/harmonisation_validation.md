# Screen harmonisation: validation (G5a; DESIGN v3 §7.1, §7.2, §7.6, §9)

## How to reproduce

```
python scripts/chronos_tkov3.py --epochs 301          # TKOv3 Chronos refit (35 min)
python scripts/harmonise_screens.py --write-draft-config   # end-to-end (23 min)
```

- Module: `src/essml/harmonise.py`.
- Tests: `tests/test_harmonise.py` (16 toy tests with known answers); the full suite passes (106 tests).
- Inputs: `configs/harmonisation_screens.json` (arm registry), `configs/control_split.json`,
  `configs/gene_classes.csv`.
- Outputs: `results/harmonisation/*`, `configs/screen_plan.csv`, and the **draft** `configs/harmonisation_frozen.json`.
  The PI commits the draft at `freeze-transfer`.
- All numbers below come from the final run of 2026-10-26 (PI decisions D10 and D12 applied).

## Blinding

- **DepMap targets.** Read for dev lines only, through `guard.load_targets`. Every anchor is a dev line; the test
  lines P31/FUJ and SKM-1 were never read.
- **Non-transformed (set-B) screens:**
  - Scaling parameters and QC come from control genes only, via `guard.load_external_screen(..., controls_only=True)`.
  - Full harmonised profiles, including the Chronos TKOv3 set-B profiles, go to `cache/harmonised/` only. They are
    logged as `harmonise_write_only` rows in `results/test_access_log.csv`.
  - Replicate ceilings are aggregate data-vs-data statistics, logged as `replicate_ceiling_aggregate`.
- No model was trained or applied. Nothing gene-level from a non-transformed screen is in `results/`.

## 1. Frozen choices (draft)

| step | choice |
|---|---|
| Guide filter | Drop guides with mean T0/plasmid reads < 30. For RPM tables (Ting) the cut-off is < 5 RPM. Retained: 52,874/71,090 (Mair MEF) to 70,556/71,090 (Olivieri). |
| LFC | Per-sample RPM (within pool for the 2 Ihry pools), plus a 0.5-read pseudocount. LFC per end replicate vs the mean T0, then averaged over replicates. |
| Gene score | Median guide LFC over at least 2 guides (DESIGN default). Mean-of-guides scoring moved anchor agreement by ≤ 0.02. Chronos was tried for TKOv3 (§4) and not adopted. |
| **Affine (D10)** | Per arm: the scaling-set median maps to **−1.4677** and the NEGv1 scaling-half median to **+0.0072**. These are DepMap-matched targets (see note below). |
| Gene-wise | EB location correction per lab/library (see below). |
| QC | Computed on the evaluation halves (342 CEG / 462 NEG). Non-Avana libraries pass with SSMD < −1.0; Avana-type libraries with NNMD < −1.0. With fewer than 10 evaluation genes in either set, QC is "not assessable" (counted as a fail). |
| ORCS / BF / MAGeCK screens | Rank-only (percentile, low = essential). |
| CRISPRi screens | Descriptive only. |

**Affine targets (PI decision D10).** Each target is the median over the 1022 dev lines of that line's median gene
effect on the same gene set. Provenance is stored in `harmonisation_frozen.json → affine_targets`, and the values are
recomputed and checked on every run. DESIGN's literal −1/0 is reported as the alternative throughout.

**Gene-wise EB correction:**
- Model: d_ag = x_ag − y_ag = δ_g + e. The prior is δ_g ~ N(0, τ²), so estimates shrink toward *no* correction.
- τ² is the cross-anchor covariance of d, which needs at least 2 fit anchors.
- Noise is inflated by κ = 2. κ was fixed a priori, not tuned; κ = 1/2/4 is reported in
  `genewise_kappa_sensitivity.csv`.
- The scale term is used only with at least 8 fit anchors. Its sampling variance is a per-gene jackknife.
- The 15 p53/DDR genes are excluded from both fit and application: TP53, CDKN1A, TP53BP1, USP28, CHEK2, ATM, MDM2,
  MDM4, PPM1D, RRM2B, ZMAT3, TNFRSF10B, FDXR, SESN1, DDB2.
- A correction is applied to a library only if two conditions hold: the LOAO gain in same-line deviation correlation
  exceeds 0.01 on that library's own anchors, and at least half of those anchors improve.

The control split uses seed 20261026.

## 2. Why DepMap-matched targets (D10): both conventions

DepMap dev lines put our scaling set (dev common essentials minus the CEGv2 evaluation half) at a per-line median of
−1.468 (5–95%: −1.55 to −1.40) and the NEG half at +0.007. Mapping the scaling set to −1 therefore compresses
external screens 1.47× relative to DepMap. DESIGN's −1 was shorthand for "map to DepMap's scale". Our scaling set is
deeper than DepMap's own common-essential control list, which is why the literal value is wrong.

| metric | DepMap-matched (used) | DESIGN −1/0 (alternative) |
|---|---|---|
| Anchors, mean dev r after affine (raw: 0.197) | **0.200** | 0.164 |
| Anchors, NN recall@1 after affine (raw: 0.78) | **0.83** | 0.61 |
| RPE1 ceiling, dev r (study level) | 0.32 | 0.43 (inflated: shared −0.32·gene-mean term) |
| Simulation, no planted gene-wise distortion: RMSE change from step (iii) | 0.056 → 0.056 (identity) | 0.150 → 0.087 (absorbs the scale mismatch) |
| Gene-wise decisions | Yusa, Wang applied at library level (not to Yilmaz: D12); TKO not | same |

## 3. Anchors: same-line agreement with DepMap, before and after each step

Metrics:
- dev r = corr(x − DepMap dev gene mean, y_line − gene mean).
- NN = rank of the true line among the 1022 dev lines (1 = recalled).

Step 2 is always leave-one-anchor-out (LOAO). The LOAO-corrected anchor profiles are cached for the §7.6 anchor skill.
Mean dev r per group (`anchor_agreement.csv` has per-anchor values; all 23 anchors pass QC):

| group (anchors) | raw | affine | gene-wise LOAO | recalled at rank 1 (LOAO) |
|---|---|---|---|---|
| Moffat TKO (DLD-1, HCT 116: TKOv1; HAP1: TKOv3) | 0.097 | 0.092 | 0.093 | 1/3 (DLD-1). HCT 116 rank 11 (raw 4, affine 15); **HAP1 rank 34** (raw 71, affine 47) |
| Yusa v1 (MOLM-13, MV4-11, OCI-AML2/3, HT-29, HT-1080) | 0.212 | 0.222 | **0.258** | 5/6 |
| Wang 2015 library, published CS (10 AML lines, K-562, Raji) | 0.202 | 0.208 | **0.291** (with scale term: 0.298) | 10/12 |
| Brunello (A-375) | 0.325 | 0.297 | n/a (1 anchor) | 1/1 |
| Avana (HT-29, Doench) | 0.225 | 0.193 | n/a | 1/1 |

Gene-wise decisions (`genewise_decisions.csv`), as LOAO gain in dev r:

| library | gain | applied? |
|---|---|---|
| TKOv1 | +0.003 | no |
| TKOv2 | +0.001 | no |
| TKOv3 | −0.003 | no |
| Yusa v1 | +0.036 | **yes**, location; no set-B count screen uses it |
| Wang | +0.083 | yes at library level, location (scale term +0.090, not more than the location gain + 0.01). **Not applied to Yilmaz 2018 (PI decision D12)**, the only set-B screen in this library |
| Brunello, Avana | – | no (1 anchor) |
| Novartis, custom DDR, CRISPRi | – | no (no anchors) |

Controls on step (iii):
- **Negative control: wrong-line references.** Each fit anchor is paired with a random non-anchor dev line. The gain
  shrinks to +0.013 (Yusa), +0.027 (Wang) and −0.005 (TKO).
- **Mechanism check: within-set permutation.** Permuting references within the fit set reproduces the real gain
  (+0.034, +0.080). A location offset depends only on the fit-set mean of the references, so this is a gene-level
  library offset rather than a line match.

**Duration check** (`duration_check.csv`):
- HAP1 T6–T18 and HT-29 D7–D25: the affine step absorbs duration. The evaluation-CEG median is stable across
  timepoints.
- HT-29 agreement peaks at D13.
- HAP1 agreement is poor at every timepoint.

## 4. Chronos refit for TKOv3 (PI-authorised, 90-min time box)

**Setup:**
- Installed `crispr_chronos` 2.3.11 + TensorFlow 2.21.0 into the Python 3.13 environment by `pip install crispr_chronos`.
  This environment change is not in `configs/pip_freeze_initial.txt`.
- One **joint** fit (301 epochs, 35 min) over all TKOv3 count screens with a usable time axis:
  - the HAP1 anchor (T0 and T3–T18);
  - Olivieri NT ×6;
  - Brown PACE002 (T18/T30) and PACE006 (T25);
  - Mair H1 MEF (T4/T8) and laminin (T4/T9/T12).
- Brown Durocher "T3" was left out: its time unit is not stated.

**Why a joint fit:** Chronos shares guide efficacy across the screens of a library, and one screen identifies it
poorly.

**Inputs:**
- Each screen's own T0 serves as its pDNA reference; set-B screens have no plasmid reads.
- Reference reads < 30 mask that guide's late timepoints.
- Cutting negative controls are the guides targeting the NEGv1 scaling half.
- Gene effects get the same per-screen affine scaling.

**Adoption rule, fixed in advance by the PI:** adopt only if the HAP1 dev r gain is ≥ +0.02.

**Result** (`chronos_tkov3_comparison.csv`, `chronos_tkov3_decision.json`):
- HAP1 dev r: 0.061 (median LFC) → 0.074 (Chronos), a gain of **+0.013, below the margin**.
- HAP1 NN rank: 47 → 21.
- Evaluation-half SSMD improved for 10 of 11 TKOv3 screens. Examples: Mair laminin −1.32 → −2.13, MEF −1.08 → −1.31,
  Olivieri −1.67…−1.95 → −2.06…−2.20.

**Decision: Chronos not adopted; TKOv3 stays on median LFC.**

Caveat: HAP1 is the **single TKOv3 anchor**, so this is weak evidence either way. The SSMD gains suggest Chronos would
help QC. Under the pre-stated rule it does not qualify, and I did not change the rule after seeing the result. The
Chronos profiles stay cached for a possible sensitivity analysis. The route can be switched on only through
`chronos_tkov3_decision.json`.

## 5. Per-screen QC (set B; evaluation-half controls; `screen_qc.csv`)

After affine scaling, the 23 QC-passing KO count screens have:
- evaluation-CEG medians of −1.35 to −2.03 (median −1.51; the DepMap dev lines give −1.36);
- evaluation-NEG medians within ±0.13 of 0.

| screen (arm) | SSMD | pass? |
|---|---|---|
| Olivieri TKOv3 NT ×6 / TKOv2 NT ×3 | −1.67…−1.95 / −1.11…−1.62 | pass |
| Brown WT Durocher, PACE002 T18/T30, PACE006 T25; p53-null TKOv2 | −1.61, −2.11, −2.45, −2.10; −1.93 | pass |
| Haapaniemi Brunello WT / null | −0.42 / −1.23 | **WT fail** / null pass |
| Hart 2015 RPE1 | −1.81 | pass |
| Bowden DDR WT / KO (focused) | −1.88 / −2.13 | pass |
| Mair MEF T8 / laminin T12 | −1.08 / −1.32 | pass |
| Yilmaz D23 | −1.79 | pass |
| Ihry D18 | −1.29 | pass |
| Shifrut CD8+ T | −0.81 | **fail** |
| Ting CD4+ T (2.5k genes) | −1.55 | pass |
| CRISPRi: Tian 2021 / Tian 2019 D10, D28 / Dräger | −0.75 / −1.03, −1.25 / −0.52 | fail / pass / fail (descriptive) |
| ORCS rank-only: fetal NSC ×4, RPE1 799/1791/1792/1932/1933, MCF-10A ×2, HUDEP-2, organoid ×2 | −1.2…−2.3 | pass |
| ORCS: Drainas 3D 1934/1935, Carnevale 1905 | −0.52, −0.62, −0.79 | **fail** |
| ORCS: Tuano HMEC | – | not assessable (2 NEG genes) |

For reference, DepMap dev lines on the same halves have a median SSMD of −1.99.

## 6. Simulated-batch control (`simulated_batch.csv`; known answer)

**Setup:**
- Anchors are 3 or 11 DepMap dev lines; targets are 20 other dev lines.
- Planted distortion: x = 3·[c + λ_g(y − c) + δ_g] + 0.4 + ε, with δ_g ~ N(0, 0.25²), λ_g = exp N(0, 0.25²) in the
  scale scenario, and ε ~ N(0, 0.15²).

**Results (matched convention):**
- RMSE vs truth: 1.17 → 0.26 after affine → **0.089** after gene-wise with 3 anchors. With 11 anchors and the scale
  term: 0.102.
- corr(δ̂, planted δ) = **0.965** (3 anchors) and 0.936 (11 anchors).
- p53/DDR genes are untouched.
- With no planted gene-wise distortion, step (iii) is the identity (0.056 → 0.056).

## 7. Replicate ceilings (observed vs observed; matched scale; study level)

Arm pairs and both directions are summarised within each unordered study pair first
(`replicate_ceilings_study_level.csv`; pairs in `replicate_ceilings_study_pairs.csv`). DepMap RPE1-SS was not
available. Values are medians over study pairs, with ranges in brackets.

| background | dev r (all genes) | dev r (CEG/NEG-free) | AUPRC (< −0.5, CEG/NEG-free) | recall @ 5% FPR, screen-specific |
|---|---|---|---|---|
| RPE1: 7 studies, 21 pairs (Olivieri, Brown ×4 screens, Hart 2015, Haapaniemi null) | **0.32** (0.17–0.56) | 0.32 (0.16–0.53) | 0.53 (0.40–0.77) | 0.36 (0.14–0.69) |
| hPSC: 3 studies (Mair, Ihry, Yilmaz); after D12 every hPSC screen is affine-only | **0.13** (0.12–0.23) | 0.13 (0.10–0.21) | 0.45 (0.41–0.49) | 0.11 (0.10–0.13) |
| hPSC before D12 (Wang correction on Yilmaz), for the record | 0.08 (0.05–0.13) | 0.08 (0.05–0.13) | 0.42 (0.37–0.45) | 0.10 (0.07–0.11) |

**RPE1.** Agreement is moderate. It is best between the TKO-lab screens: Olivieri vs Brown p53-null 0.62 CEG/NEG-free
(DESIGN scale). Haapaniemi null and Hart 2015 agree least.

**hPSC.** Agreement is low, so "fails" will be near-unreachable for hPSC (large MDE).

**D12 (PI decision).** No gene-wise correction for Yilmaz 2018. DESIGN §7.2 defines corrections per lab/library
pipeline, with anchors processed the same way as the targets. The Wang anchors are published gene-level CS (mean
log2FC), while Yilmaz is scored here from counts (median guide LFC). The mismatch was visible: applying the correction
lowered Yilmaz's agreement with Mair and Ihry (hPSC 0.13 → 0.08). For the same reason, the Wang anchors are not a
same-pipeline anchor for Yilmaz, so its §7.6 verdicts are "inconclusive only". The override is recorded per screen in
`harmonisation_frozen.json → genewise.per_screen` (`gene_wise_applied: false`, reason, decision D12).

## 8. Screen plan (`configs/screen_plan.csv`, 94 rows)

**Reachable verdicts (Q3a / Q3b):**
- **RPE1 (14 QC-passing count arms):** transfers / fails / inconclusive.
- **hPSC:**
  - Mair (TKOv3; single, weak HAP1 anchor): transfers / fails / inconclusive.
  - Yilmaz: inconclusive only (D12: no same-pipeline anchor).
  - Ihry: inconclusive only (no anchor).
- **Primary T:**
  - Shifrut and Carnevale fail QC, so they are not evaluable.
  - Ting CD4 is inconclusive only (no anchor, no ceiling).
  - So the primary-T background has **no evaluable verdict**, which affects H4.
- **ORCS rank-only:** Q3a on the rank label is inconclusive only. Q3b is not evaluable.
- **Fetal NSC, MCF-10A, HUDEP-2, organoid:** no matched expression.

Exclusions are listed in `results/harmonisation/exclusions.csv`: marker sorts, the differentiation contrast, CRISPRa,
duplicates, unusable ORCS screens, and test or non-DepMap anchors.

## 9. Deviations from DESIGN v3 and open items

**Deviations:**
- **D10 (PI):** DepMap-matched affine targets (§2).
- **D12 (PI):** No gene-wise correction for Yilmaz 2018, and no same-pipeline anchor for it (§7).
- **D-a:** A lab-level TKO group. TKOv1/v2/v3 share only 10–33k guides, and each library has ≤ 2 anchors. Decisions
  are still made per library.
- **D-b:** Wang-library anchors are gene-level CS. No guide-level anchor data exist; this is resolved by D12 (§7).
- **D-c:** QC is "not assessable" with fewer than 10 evaluation controls.
- **D-d:** The Shifrut CFSE-dilution sort is treated as a proliferation readout (H4 names it). It fails QC anyway.
- **D-e:** The blind write path and the ceilings run through the driver, logged with `guard._log`. The frozen guard
  has no API for them; suggest adding one at `freeze-transfer`.

**Open items:**
- **O1:** Chronos for TKOv3 was not adopted (+0.013 < +0.02). It rests on a single anchor.
- **O3:** The MDE and the gene-family block bootstrap are drafted in `configs/transfer_criteria.json`
  (`scripts/transfer_criteria.py`; §10).
- **O4:** There is no K562 CRISPRi anchor.
- **O5:** Commit 5181e51 contains mid-run versions of the harmonisation files. The working tree is final and
  uncommitted.
- **O6:** The environment now includes TensorFlow and crispr_chronos.

## 10. Draft transfer criteria (`configs/transfer_criteria.json`; `scripts/transfer_criteria.py`, 18 min)

**Replicate-ceiling skill.** Computed from observed screens only, using one study representative per study (the arm
with the most negative evaluation-half SSMD). It is summarised as the median over directed cross-study pairs.

**Bootstrap.** Gene-block bootstrap with 1000 resamples:
- Blocks are paralog families: connected components of Ensembl paralog pairs with identity ≥ 30% in both directions.
  This gives about 11.9k blocks, the largest with 162 genes.
- Seed 20261127, which is `split_seed` + 101 (RPE1) and + 102 (hPSC).
- MDE = 2 × SE.

| background | Q3a ceiling AUPRC (all genes / CEG-NEG-free / rank-10%) | SE | **MDE Q3a ΔAUPRC** | Q3b ceiling recall@5%FPR | SE | **MDE Q3b Δrecall** |
|---|---|---|---|---|---|---|
| RPE1 (7 studies, 42 pairs) | 0.661 / 0.544 / 0.700 | 0.008 | **0.016** (CEG/NEG-free 0.020; rank 0.023) | 0.382 | 0.019 | **0.037** |
| hPSC (3 studies, 6 pairs) | 0.516 / 0.448 / 0.338 | 0.007 | **0.013** (0.013; 0.017) | 0.118 | 0.015 | **0.031** |

Caveat: the SE reflects gene sampling only. Between-study spread of the ceiling is much larger (RPE1 study-pair
AUPRC 0.40–0.77), so "fails" (upper CI < MDE) is a strict but reachable condition.

**Statistics:**
- **Q3a.** The reported statistic is AUPRC(T) relative to the ceiling. The verdict uses ΔAUPRC vs B0, with a
  gene-block CI.
- **Q3b.** Δrecall vs B1 (and vs B0) at 5% FPR on the NEGv1 evaluation half. Screen-specific genes have a dev
  dependency fraction ≤ 0.5.
- **D7.** Set-A/B skill is averaged over the per-fold models and compared against anchor skill from single OOF anchor
  predictions.

**Verdict rules:**
- **transfers:** lower CI of skill > 0, AND skill ≥ 50% of the same-pipeline anchor skill.
- **fails:** skill ≤ 0, AND upper CI < MDE, AND anchors show positive skill, AND QC passed.
- **inconclusive:** otherwise.
- **Expression-swap cap:** if swapped-profile skill has a lower CI > 0, the verdict is capped at inconclusive.
- **Overall:** a strict majority of evaluable backgrounds (RPE1 and hPSC, so both must agree).
- **Primary T:** inconclusive by declaration (QC failed).
- **HMEC:** not assessable.

**Hypotheses:**
- **H1:** MDM2 and PPM1D get a lower predicted percentile in TP53-WT than in TP53-null RPE1 inputs, for both genes,
  with a screen-level permutation p < 0.05.
- **H2:** WT−null delta over the 15 p53/DDR genes:
  - pass if sign agreement is ≥ 70% over genes with |observed Δ| > 0.25, AND r > 0.3;
  - primary pair: Bowden isogenic WT vs KO; secondary: Brown PACE002 WT vs Olivieri TKOv3 TP53−/−.
- **H3:** at least half of POU5F1, LIN28A, NANOG, SOX2, DPPA4, L1TD1, ZFP42 and PRDM14 have a predicted effect
  < −0.5 in hPSC with the EC lines in training. If this holds only with EC lines in training, the claim is
  "narrow EC transfer".
- **H4:** inconclusive by declaration.
