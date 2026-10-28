# Red-team review of REPORT.md (draft at commit 957edfb)

Reviewer: an independent agent, 2026-10-28, read-only. The harness blocked subagent file writes, so the PI transcribed
this from its hand-back. Line numbers refer to the draft at 957edfb. Each item's resolution is given in brackets; all
items were applied (D29). The replacement texts the reviewer proposed were applied with their source tags; they are
summarised here rather than repeated.

The reviewer checked REPORT.md against:
- PROMPT.md, RESEARCH_TREE.md (D0–D28), both independent reviews, reviews/implementation_check.md and
  results/scope_ledger.csv;
- the results JSONs (test, setA, setB RPE1/hPSC/summary, g2, g3 kernel/nonlinear/M5, mutation_dependencies,
  predictability, timing, registry);
- the figures, configs/transfer_criteria.json, the run configs and scripts/predict.py.

Its overall assessment: the tagged numbers are right. The problems were what the prose around them says,
interpretations the report left out, and statements about the shipped tool.

## Critical
1. **The shipped tool is not the evaluated headline, and the report's advice cannot be followed** (draft lines 36–39,
   339–340, 554–555, 719–720, 735–737, 833–835).
   - scripts/predict.py takes only `--expression`. For m5_avg it sets copy number, mutations, hotspots and fusions to
     diploid wild type and metadata (including lineage) to unknown.
   - The test score (0.319) and the NRAS recovery (0.631) come from M3 reading real hotspot, CN and lineage features.
     The expression-only configuration users run was never evaluated on cancer lines (implementation check R1-PI-5).
   - The draft advised "use m5_avg with genomic inputs", which no user can do.

   [Resolved: a post-hoc dev-CV diagnostic of the expression-only use (results/predict_expression_only_cv.json) and
   a pre-stated default rule (D29). The abstract, Results, Methods and Limitations were rewritten; the predict.py
   warnings were corrected.]

## Major
2. **Negative controls misreported** (203–209). Two of the three stricter controls missed a pre-stated expectation:
   perm_lineage P1 0.127 fell below its band (g3_kernel_summary controls/1, pass False), and controls/2 exceeded its
   bound. The headline/m1 attribution was also wrong. [Applied.]
3. **Set A's Q3b omitted** (48–49, 466–471, 530, 622–625). Q3b was positive in 5/5 RPE1 lines, every CI above 0,
   median 95th percentile of the cancer test lines. "Not detectable" was one-sided. [Applied.]
4. **"Transfers" does not single out the model** (44–49, 489–499, 515–516, 622). [Applied.]
   - B2 meets the same conditions.
   - T shows no skill on external screens of four cancer anchors.
   - The MDE is driven by between-study variance; with fixed pools it is 0.020, which the RPE1 interval would not
     lie below.
5. **The RPE1 result never tested expression harmonisation** (447–449). T's RPE1 input was DepMap's own expression
   (the mean of the six RPE1-ss profiles). [Applied.]
6. **Q3a mislabelled as core-essential transfer** (93–94, 455–456). The scaling halves of the core-essential lists are
   removed from its universe; the gene-mean comparator's AP is the more direct evidence. [Applied.]
7. **Task gaps hidden in the abstract and Discussion** (44–49, 552, 727–728). No iPSC or iPSC-derived screen was
   evaluated (the hPSC screen is H1 ESC), and no primary haematopoietic screen passed QC. [Applied.]
8. **Essentiality AUPRC barely interpreted** (27–49, 382–390, 607–612). P2a is close to its ceiling for the gene mean
   alone (0.905 vs 0.915), and selective nAP is modest (0.184). [Applied.]
9. **Published comparison selective** (289–293, 627–645). [Applied.]
   - Only the lowest of three r ≥ 0.2 counts was quoted.
   - TrPLet on DeepDEP's genes was omitted.
   - "Non-nested selection accounts for a substantial part of such gaps" was unsupported.
   - R4's original gap is an upper bound (D19).
10. **Paralog-loss claim goes beyond the data** (327–329). E-only AUPRC was 0.180 / 0.305 / 0.338 for ARID1B / CREBBP /
    STAG1; AUPRC is not prevalence-normalised. The D24 prior-exposure disclosure was missing. [Applied.]
11. **Missing disclosures.** [Applied.]
    - **(a) Predictability:** the tiers were seen for m1 first (D22b); the null covers the M1 path only (D22e); the
      pre-G6 anchor reads (D22f).
    - **(b) Review items implemented only in part:** R1-PI-5, R1-ST-2, R1-ST-8, TF-F16.
    - **(c) Registry provenance:** the headline CV row and its inputs have tree_clean = 0, so provenance rests on the
      D25b dry run.
    - **(d) Bootstrap bias:** the Q3b recall bootstrap bias reaches about −0.017.

## Minor
12. "Confidence" vs dev mean r in the Spearman statistic; give the approximate CI or drop the sentence. [Applied.]
13. Tier P1 is over all genes in the tier, not variable genes; also in the Figure 1B caption. [Applied.]
14. The G2 ablation baseline is the pre-D14 B1-residual m1_E (0.2997). [Applied.]
15. The abstract over-attributed signal to individual modalities; hotspots and fusions add nothing on top of E+C+M.
    [Applied.]
16. Feature descriptions: the headline M2 is paralogs only, and M3's pool includes the whole metadata block. [Applied;
    the PI also verified that M2 uses the gene's own CN and mutation status.]
17. The M4 row is a fold-0 value; add fold-0 M1. [Applied.]
18. Contradictions: there is no G5 gate file; the MS proteomics wording in Future work; "fitted exactly" vs merged
    patterns. [Applied.]
19. The quoted reviewer wording was incomplete ("top decile", "about one twentieth relative"). [Applied.]
20. The headline-on-RPE1 secondary lacked its CI and comparison with T and B2. [Applied.]
21. "Beats nearest neighbours on most genes" had no number. [Applied: +0.070 P1.]
22. Runtime: M3 used 6 worker processes; Figure 3 understates M2's cost. [Applied.]
23. The test refit trained on 1,022 lines against about 818 per CV fold. [Applied.]
24. Non-essential P1 magnitude, about a third of variable-gene P1 (D16). [Applied.]
25. "Same feature types where data permits" was not addressed. [Applied; whether DepMap genomics for the RPE1
    sublines could have been used is stated as not examined.]
26. Readability: gates undefined; "PI" and "external PI" are AI agents; undefined terms (TPM, TKO, ORCS, AUROC, HAP1,
    RPE1, D-G3-n); reference 14 uncited. [Applied.]

## Verdict (reviewer)
The report can be published once items 1–11 are fixed. All are wording or disclosure changes that need no new
analysis. Item 1 must be fixed in the text; evaluating the expression-only m5_avg in dev CV would be a new post-hoc
step and should be labelled as such if done. [It was done, and it is labelled as post hoc (D29).]
