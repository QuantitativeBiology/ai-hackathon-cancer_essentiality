# Round 1 review: external PI

Reviewer: external PI (computational biology). Object: DESIGN.md v0. Read: PROMPT.md, DESIGN.md, LITERATURE.md, docs/data_audit.md, docs/lit/*.md.

## 1. Verdict: major revision

The cancer half of the design is sound. It has patient-grouped CV, a locked test set, the lineage-mean baseline, per-gene metrics split by gene class, and a ladder of models. Most of my changes there are clarifications. The non-transformed half, which is half of the prompt, is not decisive as written. Its core-essential test does not use the external screen data at all. RPE1 sits in both the anchor set and the target set. Anchor lines are predicted in-sample. And there is no rule that turns the results into "transfers", "fails" or "inconclusive". The design also has no ordered work plan or go/no-go gates. The literature (prior_methods §5) asked for within-lineage r, and the design drops it.

## 2. Proposed changes

**1. §7.4, core-essential test (must).** As written, AUPRC of *predicted* effects for CEGv2 vs NEGv1 never uses the external screen. B0 gets the same score, so the test is uninformative. Replace it with two parts. (a) Observed transfer: per screen, the recall and precision of dev-defined common essentials among genes that are essential in the harmonised screen. Report this next to the held-out CEG/NEG AUROC (screen quality) and the cross-study ceiling. (b) Model transfer: AUPRC of model predictions against the *observed* labels (all genes, plus the common-essential subset), always paired with B0. Also list the common essentials that are lost or gained, for example in hPSC vs cancer, and report whether the model predicted each shift. Why: as written, "core essentials transfer" is true by construction.

**2. §7.2 step 3 and §7.4, anchor circularity and leakage (must).** (i) Remove RPE1, and every non-transformed line, from the anchors used to fit corrections that are then applied to non-transformed screens. Otherwise the harmonisation is fitted to make external RPE1 match DepMap RPE1. Use RPE1 only as the cross-study ceiling. (ii) Any anchor used as a positive control (HCT116, DLD-1, HT-29, A-375, K-562, AML lines) must be predicted by a model refit without that anchor or its PatientID mates (for example DLD-1/HCT-15, K-562/P2UR). Most anchors sit in dev or test, so as written their predictions are in-sample.

**3. §7.4, a pre-registered decision rule (must).** The primary contrast should be a within-batch difference-in-differences: model skill minus B0/B1 skill on the non-transformed screen, against the same quantity on transformed anchors screened in the same lab and library (Hart 2015 RPE1 vs HCT116/DLD-1; TKOv3 hPSC/RPE1 vs HAP1). Batch effects largely cancel in this contrast. Decide per cell background and per gene class:
- **Transfers:** skill over B0 has a gene-bootstrap CI above 0, and the drop from the anchor skill is less than a pre-set fraction, for example 50%.
- **Fails:** the CI includes 0 while the anchors in the same batch show skill and the cross-study ceiling is at least 2× the minimum detectable effect.
- **Inconclusive:** anything else.

Report the minimum detectable effect from the ceiling before applying the model.

**4. §7.4, a negative control for specificity (must).** Swap expression profiles between external screens (predict screen A from screen B's expression), and also predict from a centroid DepMap profile. If the context-specific correlation does not fall, the signal is shared batch or proliferation structure, not cell-specific prediction.

**5. §5/§6, a defined "transfer model" (must).** External cells are near-diploid with few mutations, and CN/M/F are often unprofiled. Pre-declare which features are available per external screen. Select, in dev CV, the best model that uses only those features (E, E-derived, lineage mapping). In addition, test the full model in dev CV with C/M/F set to neutral values on held-out lines. Pre-register the lineage mapping for each external sample (T cells → Lymphoid? hPSC → none → B0). The prompt says to profile "with the same feature types where data permits", so this choice must be explicit, not improvised later.

**6. New §12, ordered work plan and gates (must).**
1. Verification harness: toy, positive and negative controls must pass before any real metric.
2. B0/B1/B2 and M1 on dev. Gate: M1 beats B1 on P1, by a paired bootstrap over patients.
3. Started now, in parallel, because it is the long pole: external acquisition plus the anchor round-trip test. Gate: same-line recall at k=1 improves after harmonisation, otherwise use rank-only metrics for that library.
4. M2, then ablation.
5. M3, only if M2 leaves a gap on variable genes. M4 is time-boxed (for example 1 day) and dropped if it does not beat M2 in the first fold.
6. Freeze finalists, test once, refit on dev+test.
7. Set A, then set B.

Gate for set B: a category with no QC-passing screen that has matched expression is reported as "not testable", not filled with proxies.

**7. §3/§7, set A heterogeneity (must).** Split set A by transformation status using ModelTreatment and the literature: hTERT-only or near-normal (RPE1-SS, WPE1-NA22?, Hs 729) versus viral, SV40 or adenoma (C8166, HG-3, HA1E, BPH-1, organoid). Treat the 5 RPE1 sublines as one background for inference. Drop "per-gene correlation across ≥4 samples" as a test: across RPE1 sublines it measures aneuploidy noise, and across mixed lines it measures lineage. Keep it as descriptive only.

**8. §4, within-lineage per-gene r (should).** Add it as a secondary metric, computed on residuals after the training-fold lineage mean. It is the cleanest guard against a headline driven by lineage. Make "Δ P1 over B1 with a paired CI" the headline number, not raw P1.

**9. §4, how P1 is computed (should).** State that per-gene r is computed within each held-out fold and then averaged, with genes kept only if they vary in that fold. Base the 1-SE rule on *paired* fold differences, not on the SD of each model across folds.

**10. §4, P2a framing (should).** P2a is driven by common essentials, and B0 scores close to the ceiling on it. Always report it as Δ vs B0 and never alone. P2b vs B1 is the essentiality headline.

**11. §7.1/§7.4, readout strata (should).** Pre-declare rank-only metrics for screens with only BF or MAGeCK scores (Hart 2015 via ORCS, Yilmaz, Shifrut). Keep CRISPRi neurons and microglia as a descriptive stratum. Post-mitotic survival is not proliferation, so missed cell-cycle genes there are expected and are not evidence of model failure.

**12. §7.3, applicability domain (should).** Report each external profile's distance to its nearest dev lines (and its percentile against dev-to-dev distances) next to its performance. This separates "out of domain" from "fails".

**13. §2, missing modalities (should).** Send the user the exact download list now (data_resources.md: OmicsSignatures, harmonized_Sanger_MS_2022, RRBS methylation, metabolomics), because the ablation that the prompt requires depends on it. Pull Nusinow 2020 proteomics (open URL, about 375 lines) straight away as a partial substitute.

**14. §2, complexes (could).** CORUM is unavailable and CC BY-NC. Use EBI Complex Portal (curated, not text-mined) instead.

**15. §7, final model for external use (could).** State that the frozen finalist is refit on all 1277 cancer lines with the hyper-parameters fixed, and that no external data feed back into model choice.

**No changes needed:** the split design and locking, the leakage tests, the NaN policy, the exclusion of screen-proxy metadata, the runtime harness, the B1 shrinkage, and the M1 residual formulation.
