# Round 3 — fresh external review (independent; DESIGN v2)

Reviewer: senior computational biologist, no prior exposure to v0/v1.

## Verdict
**Accept with revisions.** The design answers PROMPT.md's questions and the guard/freeze machinery is better than
most published work. Two problems remain: (a) a handful of thresholds are referenced but never set, plus four
genuine cross-section contradictions; (b) the plan has accumulated more arms than 12 cores and one project can
execute — the non-transformed half (the part PROMPT.md cares most about) is the part at risk if the cancer half
overruns. Changes below; cuts are items 8–14.

## Internal consistency

1. **§3 vs §4 — gene classes (MUST).** §3 recomputes classes per training fold; §4 claims the fixed
   `eval_genes.csv` means "fold differences never reflect which genes were scored". Untrue: *variable* membership
   (SD > 0.25) moves between folds, and P1 is defined on *variable* genes. Freeze one `gene_classes.csv` from all
   1022 dev lines (training-only, so not leakage) and use it for CV, test, set A and set B; report the per-fold
   recomputation as a robustness check only.
2. **§4 — `eval_genes` does not apply to A/B (MUST).** "≥50 observed test lines" is unsatisfiable for set A (15
   lines) and for external screens, yet §4 says "all metrics use this set". Add `eval_genes_external` = dev-evaluable
   genes ∩ genes measured in that screen/library, and say §4's set governs cancer dev/test only.
3. **§9 guard vs §7.2/§7.6 (MUST).** "Exactly one logged test access per finalist" blocks two required operations:
   harmonisation must read set-B *control-gene* scores **before** criteria are committed, and G6 applies frozen
   models to screens that arrive at different times. Rewrite: control/anchor-gene reads on set B are a logged
   pre-freeze exception (like `observed_counts`); after freeze, one access per (finalist × screen).
4. **§5 frozen models vs §7.6 anchor ratio (MUST).** External screens are predicted by the average of 5 fold-models;
   anchor skill comes from single out-of-fold models. An ensemble is strictly better, so "skill ≥ 50% of anchor
   skill" is biased toward *transfers*. Fix: predict each anchor with the average of the 4 fold-models not trained
   on it (and state the 4-vs-5 residual). Feasible — I verified HAP1, HCT 116, DLD-1, A-375, K-562, HT-29 are all
   `dev` (folds 2/3/4) in `configs/splits.csv`; MCF-7 and Hep G2 are `test`, as §7.2 says.
5. **Thresholds referenced, never set (MUST).** §7.1 QC "NNMD must pass" (no cutoff); §7.6 count rule (given only as
   "e.g. ≥2/3"); §7.6 MDE (no formula); §7.2.3 "strong shrinkage" (no prior) and "improves anchor agreement" (no
   margin); §7.5 H1–H4 say "thresholds in §7.6", but §7.6 covers only Q3a/Q3b. Put numbers in
   `configs/transfer_criteria.json`, or relabel H1/H3/H4 as descriptive (see item 9).
6. **§7.6 count rule must be per biological background (MUST).** Hart 2015 / Brown 2019 / Olivieri 2020 are three
   datasets of *one* genetic background; a 2/3-of-screens rule can be carried by RPE1 alone, and §7.2.5 already uses
   those same screens as the replicate ceiling. Count distinct backgrounds (RPE1, hPSC, primary T) and report
   per-background verdicts.
7. **Prioritise set A over set B (SHOULD).** Set A is the only batch-free transfer test and needs none of §7.2/§7.3.
   Make Q3's answer explicitly two-part (within-DepMap non-cancerous, then external harmonised) and require set A to
   be run and reported before set-B verdicts. This also de-risks items 11–12.

## Cuts / deferrals

8. **§5 — drop M4 (factorised NN) and the reduced-rank/low-rank diagnostics (MUST).** A day per go/no-go plus an MLP
   pre-check, against LITERATURE §1 saying linear models are hard to beat; move to REPORT.md "remaining ideas".
9. **§5/§7.5 — drop R3 (RF on 500 genes) and the P proteomics arm (375 lines → underpowered ablation) (SHOULD).**
   Keep R1, time-boxed R2, R4 with its non-nested optimism run. In §7.5 keep H2 (isogenic, batch-cancelling) and
   H1 as the only scored hypotheses; H1′/H3/H4 descriptive.
10. **§5 M1b — keep the 1/√p_b default weights; defer coordinate-descent block weights (SHOULD)** unless the §6
    ablation shows modality weighting matters.
11. **§7.3 — cut candidates (iv) rank features and (v) Celligner-style (SHOULD).** Celligner without tumour data is
    its own project; keep none / quantile / quantile+gene-wise (the GTEXDEPMAP recipe).
12. **§7.2.1 — default to median guide LFC; fit fresh Chronos only where an anchor shows it changes a verdict
    (SHOULD).** "A few hours each" × N screens is the largest unbudgeted cost in §10.
13. **§4 — keep the LEAP gene set only in the published-comparable block (COULD).**
14. **§9 — one red-team pass at G7 plus the independent reproduction of headline numbers; drop per-finding
    red-teams (COULD).**
15. **§10 (SHOULD).** Table omits ablation count × M1, the 2000-resample bootstrap, harmonisation and Chronos
    refits. Add per-gate wall-clock budgets and a total.

## Remaining misleading-headline risks

16. **§7.4 label circularity + a structurally null Δ (MUST).** Positives are "harmonised effect < −0.5" while the
    affine map forces DepMap common essentials to −1, so the CEG-evaluation-half labels are partly produced by the
    scaling; and ΔAUPRC vs B0 is near-zero by construction because B0 already ranks common essentials almost
    perfectly. Add a scale-free sensitivity label (bottom-q% of that screen's gene-score ranks, q set from its own
    CEG/NEG gap), and pre-declare that Q3a's evidence is absolute AUPRC vs the replicate ceiling plus the
    lost/gained-essential analysis — not Δ vs B0.
17. **§4 selection rule (SHOULD).** δ = 0.005 with "> 1 SE" is weak for naming a single headline model; require
    ≥ 2 SE, else the simpler model keeps the headline.
