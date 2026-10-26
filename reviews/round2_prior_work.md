# Review round 2: use of prior work (DESIGN.md v1)

Reviewer: prior-work agent. Date: 2026-10-26.
Inputs: DESIGN.md v1, reviews/round1_response.md, LITERATURE.md (amended).

## Verdict

**Accept with minor revisions.** Most of my round-1 items are now in the text: rung R, the published-comparable
block, the corrected novelty claims in LITERATURE §1, within-lineage r, the extra M2 features, the proteomics arm
and model T.

The rejection of concatenated out-of-fold P1 (#3) is correct. Pooling fold-specific gene means gives B0 a spurious
negative r. Keeping pooled values only in the published-comparable block, flagged, is the right compromise.

Making RRR diagnostic-only (#2) and using the NaN group as a library proxy (#6) are acceptable, given the
algorithms evidence and that ScreenMap is not available locally.

Two accepted items are missing from the v1 text, and one new inconsistency has appeared.

## Remaining changes

**1. §7.3, reconcile "expression transfer variants — A" with the text. (must)** The candidates (i)–(iv) contain
neither Celligner nor Shi's GTEXDEPMAP recipe. GTEXDEPMAP used quantile normalisation + ComBat, without cPCA, for
non-malignant tissue (Shi 2024, PMC11358024). Two options:
- Add (v), quantile normalisation + ComBat with the external lab/study as batch. Candidate (iii) may be declared
  equivalent if it is a gene-wise location/scale fit. Also add (vi), Celligner-style cPCA + MNN (Warren 2021,
  PMC7782593).
- Or log a rejection in round1_response with a reason. A defensible one: cPCA targets tumour purity and stroma axes
  that cultured cells lack, and MNN is unstable with a handful of external samples.

Either way, cite GTEXDEPMAP in §7 as the nearest prior art. It was validated only qualitatively, with no
normal-cell screens.

**2. §7.4 and §7.5, add B2 (kNN in harmonised expression space) as a transfer comparator next to B0/B1 and the
nearest-lineage comparator. (must)** It was accepted in R1 (#7) but is absent. It is the natural similarity-transfer
baseline, and it costs nothing given the applicability-domain distance already computed in §7.3.

**3. §4/§5, put rung R on matched gene sets. (should)** R2 runs on *variable* genes and R3 on a random 500 variable
genes.
- State that every R-vs-M comparison uses per-fold P1 on the *identical* gene subset. For R3, that is M1/M2/headline
  restricted to the same 500 genes.
- Report the subset's B1 skill so the reader can see it is representative.
- In the published-comparable block, report R models on the intersection of their gene set with the DeepDEP and
  LEAP lists, and give the intersection sizes.

**4. §7.1, clarify the CRISPRi anchor choice. (should)**
- Using RPE1 (Replogle 2022) as a KO↔CRISPRi anchor conflicts with §7.2.3, which excludes RPE1 from anchors as
  circular. RPE1 is an A1/set-B evaluation background.
- Use K562 CRISPRi growth screens only (Horlbeck 2016). If RPE1 data are used, exclude every RPE1 evaluation from
  any CRISPRi-calibrated result.
- Use only fitness (growth) readouts. Perturb-seq transcriptional phenotypes are not fitness.

**5. §2 Prior knowledge: replace "CORUM if obtainable" with hu.MAP 3.0 and EBI Complex Portal. (could, repeated)**
Both are non-text-mined and openly licensed. CORUM is CC BY-NC and was unreachable.

No other change needed from the prior-work side.
