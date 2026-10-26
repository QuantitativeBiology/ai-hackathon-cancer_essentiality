# Round 2 review: external PI

Object: DESIGN.md v1 and reviews/round1_response.md.

## Verdict: minor revision

**(a) My round-1 items.** All 13 must/should items were adopted, and most are implemented well: T, the set A tiers, the profile-swap null, the core-essential test against observed labels, the RPE1 exclusion from anchors, and §12. Two are only partly done. The decision rule in §7.6 defines "fails" but not "transfers", and it names no statistic. The use of anchors as a comparator has a new circularity (change 2).

**(b) New problems.** These come from the anchor-correction loop, from the H2/H3 model variants, and from the unstated training set for frozen models.

**(c) Is Q3 decisive?** It is decisive for core essentials. For context-specific dependencies it is decisive only for screens that have a same-pipeline anchor. That is probably RPE1 (Hart 2015 TKOv1 with HCT116/DLD-1, and TKOv3 with HAP1) and possibly H1. Primary T-cell, HSPC and CRISPRi screens will almost certainly be "inconclusive" under the rule. That is acceptable, but it should be declared before the run, not found afterwards.

## Remaining changes

**1. §7.6, a complete rule (must).** For each question, name one primary statistic:
- **Q3a:** ΔAUPRC of T vs B0 against the observed labels.
- **Q3b:** Δ recall at 5% NEG FPR of T vs B1, using the *evaluation* half of NEGv1.

Then define all three verdicts.
- **Transfers:** the lower bound of the gene-bootstrap 95% CI of Δ is above 0, **and** Δ_nontransformed ≥ f × Δ_anchor, with f fixed in the json, for example 0.5.
- **Fails:** as written now, plus the condition that the ceiling-derived minimum detectable effect is smaller than Δ_anchor.
- **Inconclusive:** otherwise.

If the profile-swap control shows skill (matched − mismatched CI includes 0), the Q3b verdict for that screen is capped at "inconclusive". The json may set the numbers, but the DESIGN must fix the form of the rule.

**2. §7.2 step 3 and §7.6, anchor skill leaks the correction (must).** With 1–3 anchors per library, a gene-wise correction fitted with anchor X largely absorbs X's own deviations. That inflates Δ_anchor, which then makes "transfers" harder to reach and makes a "fails" verdict easier to justify. Compute every anchor-side quantity used in §7.6 from the **leave-that-anchor-out** corrected profile. Also, when a library has fewer than 8 anchors (Pacini 2021), restrict the correction to location only.

**3. §7.1, a pre-declared feasibility table (must).** Before G6, commit `configs/screen_feasibility.csv` with, per screen: library, lab, same-pipeline anchors (IDs and count), readout (counts/LFC/BF), matched expression (yes/proxy/none), and which verdicts are reachable ("fails" needs ≥1 same-pipeline anchor). This turns the inevitable "inconclusive" results for T-cell, HSPC and CRISPRi screens into pre-registered expectations instead of post-hoc excuses.

**4. §5 T / §9, which model and training set are frozen (must).** State that the models applied to set A and set B are the G4 finalists **refit on dev+test (1277 lines)** with frozen hyper-parameters, or else dev only, and choose one. Anchor predictions are dev-CV out-of-fold, so they come from a different fit than the set-B predictions. Either use the dev-only fit for both, or report anchor skill under both. Also count T and the H3 "no-EC" variant in `finalists.json`; the current limit is headline + 2. They must be frozen at G4, not trained after unblinding.

**5. §7.5 H2 (should).** In the TP53 counterfactual, flipping only the TP53 flag ignores the p53-target expression score, which carries most of the learned signal. Use arm-specific expression where it exists. Otherwise, report the flag-only delta as a lower bound and say that "fails" cannot be concluded from it.

**6. §3 set A (should).** Define the per-line metric that is ranked within the cancer-test distribution. It should be the deviation correlation on selective genes plus ΔAP vs B0, not the raw profile correlation, which is dominated by gene means. A1 is in effect one background, so add a sentence saying that set A cannot support a "fails" verdict alone.

**7. §5 T lineage mapping (should).** Commit the mapping from external sample to OncotreeLineage (or "none") in `transfer_criteria.json`. State which form is primary for T: mapped lineage or no lineage. Report the other as a sensitivity analysis.

**8. §7.5 CN-collateral test (could).** T sets CN to diploid, so this test must use the full-feature model on RPE1-SS lines with their real CN. Say so.

**No further changes needed:** §2 targets/NaN handling, §4 metrics and selection rule, §6, §8, §9 safeguards, §10 and the §12 ordering.
