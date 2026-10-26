# Review round 1 — biological validity (functional genomics / non-transformed screens)

Scope: §7, §4 labels, §2 derived features, §3 set A. **[local]** = computed here from 26Q3
`CRISPRGeneEffect.csv` + `Model.csv`.

## 1. Verdict

**Revise.** §7.2 harmonisation is strong — better specified than most published cross-screen work. The
*biological* half of §7 is not yet a test. Three substantive defects: (a) §4's gene classes structurally
exclude the lineage-restricted dependencies §7 exists to probe; (b) §3's "non-transformed set A" is mostly
virally transformed or oncoprotein-immortalised, so the headline claim would be mislabelled; (c) §7.4 has no
nulls or comparators and two hypotheses are near-tautological. All fixable in scope; two natural experiments
already in the local data make the fixes cheap.

## 2. Proposed changes

**1. §4 — add a rare-selective gene class. MUST.** *Rare-selective* = ≥3 lines with effect < −0.5 **and**
<5% of lines; use it in P2b and every §7.4 context test. **[local]** POU5F1 is < −0.5 in only 2/1361 lines
(NCCIT −0.88, NTERA2-cl.D1 −0.89, the two deepest in 26Q3); SALL4 min −1.35, LIN28B −1.21. The 5% floor
therefore labels pluripotency and most lineage-TF dependencies *not selective*, so P1 and P2b are blind to
the genes §7.4 pre-specifies. Check coverage first: NANOG is NaN in all 17 non-cancerous models **[local]**.

**2. §3/§7 — tier set A by transformation mechanism; stop calling all 15 "non-transformed". MUST.** Add a
`TransformationMechanism` column: T1 (hTERT, p53/Rb-competent) = five RPE1-SS sublines, plus HS729 only if
provenance confirms non-neoplastic origin; T2 (viral-oncoprotein immortalised) = HA1E (TERT + SV40 large T,
https://cellosaurus.org/CVCL_VU89) and the prostate panel (SV40-LT / HPV-E6E7 or tumourigenic); T3
(transformed) = C8166 (HTLV-1), HG3 (CLL), ARH-77 and IM-9 (EBV-lymphoblastoid, not myeloma —
https://cellosaurus.org/CVCL_1072, https://cellosaurus.org/CVCL_1305), adenoma organoid. Primary claim = T1;
T2/T3 are an "immortalised" secondary stratum. SV40 LT and HPV E6/E7 functionally inactivate p53 and Rb —
the exact axis §7.4 tests: **[local]** HA1E has no TP53-WT signature (TP53 −0.47, MDM2 −0.07), RPE1-SS48 has
TP53 +1.85, MDM2 −1.98. Pooling averages the effect away.

**3. §7.4 — split the TP53 hypothesis into learnable and unlearnable arms. MUST.** H1a (*learnable*): MDM2
and PPM1D dependency — DepMap has the TP53-WT contrast, so the model must predict these. H1b (*screen
diagnostic, not a model test*): positive effects for TP53, CDKN1A, TP53BP1, USP28, CHEK2, ATM (add the last
three; USP28–53BP1–p53 is the axis). Score by sign/rank plus an observed-vs-predicted slope, not magnitude:
**[local]** TP53 KO effect is +1.6 to +3.2 across RPE1-SS vs a pan-cancer mean of +0.38, so the suppressor
tail sits 4–8× outside the training range and any shrunk model underpredicts it. Extrapolation, not failure.

**4. §7.4 — use the isogenic TP53 pairs as a counterfactual. MUST.** For Brown 2019 (genome-wide RPE1
TP53-WT vs null, one lab/library) and Bowden 2020, compute the observed isogenic delta profile, then flip
only the TP53 feature in the model input and compare the predicted delta. Differencing within
lab/library/duration/background cancels every batch term §7.2 models, so this is the one §7.4 test whose
result cannot be blamed on harmonisation; it also checks TP53 status is actually *used*.

**5. §7.4 — add a profile-swap null and a lineage comparator to every context statistic. MUST.** Permute
which harmonised expression profile pairs with which screen to build a null for the
"pred − mean vs obs − mean" correlation, and report nearest-DepMap-lineage mean deviation as comparator
(epithelial/fibroblast for RPE1, germ-cell for hPSC, lymphoid for T cells). The gene-mean baseline gives
deviation ≡ 0, so §7.4 has neither null nor comparator and a positive correlation could come entirely from
lineage identity. §8's permutations cover training only, not transfer.

**6. §7.4 — define the headline statistic operationally. MUST.** *Screen-specific dependency* = harmonised
effect < −0.5 externally AND > −0.5 in ≥95% of DepMap lines; report recall at fixed FPR on held-out NEGv1,
model vs gene-mean vs lineage-mean vs cross-study ceiling. The question needs one pre-registered number;
correlation over *selective* genes answers a different one and, per item 1, omits the relevant genes.

**7. §7.4 — pluripotency: make a germ-cell-tumour leave-out the mechanism test. MUST.** **[local]** 26Q3
screens 7 germ-cell/embryonal-carcinoma lines (NCCIT, NTERA2-cl.D1, PA-1, 1411H, GCT27, SUSA, SCH). Train
with and without them; report whether hPSC pluripotency predictions (POU5F1, SOX2, SALL4, LIN28A/B, DNMT3B,
TDGF1, PRDM14) collapse. This converts "pluripotency regulators where measured" into a falsifiable answer to
*could a cancer-trained model ever predict them*: if yes, it is narrow transfer from two EC lines, not
general capability. It also bounds power — n≈2 dependent lines cannot support a per-gene model.

**8. §7.1/§7.4 — pre-specify the post-mitotic divergence and a readout rule. MUST.** For iPSC-derived
neurons, split genes into *proliferation-coupled* (replication, mitosis, spindle, nucleotide synthesis) and
*maintenance* (translation, proteostasis, mitochondrial, lysosomal), predict a priori that the first loses
effect and the second is retained, and evaluate separately. Include only proliferation/survival-dropout
arms; FACS-marker arms (Leng reactive state, Dräger phagocytosis, Shifrut CD25/IFNG sorts) are excluded or
reported apart. Chronos and DepMap model growth rate; in a non-dividing population the gene-effect→dropout
mapping changes qualitatively — the most predictable context effect available, currently unstated.

**9. §7.2 — give the CRISPRi stratum an anchor and a role. MUST.** Fit a KO→CRISPRi attenuation map on
CRISPRi screens of DepMap lines (K562, Horlbeck 2016 doi:10.7554/eLife.19760; RPE1 KRAB-dCas9
common-essential screen, Replogle 2022 https://pmc.ncbi.nlm.nih.gov/articles/9380471),
leave-one-anchor-out validated; attenuation is expression- and TSS-dependent, hence not affine. State the
positive role: CRISPRi makes no DSB, so this stratum separates cutting/p53 toxicity from genuine
non-transformed biology, and conversely H1b cannot be tested there. Without an anchor no CRISPRi number is
interpretable.

**10. §7.3 — fit, not merely validate, the expression correction; add a rank variant. MUST.** (i) Fit a
gene-wise location/scale map on matched DepMap/ENCODE lines (K562, HepG2, MCF-7, A-549, HCT116),
leave-one-line-out validated, mirroring §7.2 step 3. (ii) Pre-register a rank-transformed-expression variant
as the external-transfer configuration. (iii) Require expression from the **same activation state and
timepoint** as the screen (stimulated CD8+ T cells, CD34+ in expansion — not resting DICE/Blueprint), with
sensitivity across candidate references. Quantile normalisation is a within-sample monotone remap: it fixes
platform scale but not the composition shift between a cancer transcriptome and hPSC or stimulated primary
cells, and it distorts the lineage-marker genes that drive prediction. A rank model is invariant to monotone
distortion; reference choice is otherwise a hidden analyst degree of freedom.

**11. §7.2 — gate short screens on dynamic range before affine scaling. MUST.** Record NNMD and the raw
control gap |m_E − m_N|; below a pre-set threshold restrict that screen to rank metrics and AUPRC, and
report the scale factor applied. Shifrut's ~6-day T-cell screen permits few doublings, so common essentials
are shallowly depleted and dividing by a small gap amplifies noise into apparent effect sizes. §7.2 step 4's
claim that duration is "covered by affine scaling" is wrong in general: per Dempster 2019 assay length acts
gene-set-specifically (late dependencies deeper, early ones matched).

**12. §2 — replace the TP53 functional score with a p53/Rb-integrity composite. SHOULD.** TP53
damaging/hotspot + TP53 CN loss + MDM2/MDM4 amplification + CDKN2A + ATM/CHEK2/PPM1D + a viral-oncoprotein
flag (required by item 2); down-weight basal p53-target expression, a weak proxy for p53 competence in
unstressed cells and atypical in hPSC (low p21 despite high p53 protein), so it mis-ranks the very samples
§7.4 targets.

**13. §7.4 — demote the oncogene-addiction test; add a CN-collateral test. SHOULD.** "KRAS/BRAF/CTNNB1
predicted near 0" is near-tautological: the mutation feature is absent so the prediction is 0, the oncogene
is inactive so the observation is 0 — one input, no information. Move it to the negative controls; replace
with collateral/paralog-buffered dependencies requiring CN loss (SMARCA2 given SMARCA4 loss, ARID1B/ARID1A,
STAG1/STAG2, ME2/ME3, ENO1), predicted and observed null in near-diploid cells.

**14. §7.4 — drop "per-gene correlation across ≥4 set-A samples". SHOULD.** The five RPE1-SS lines differ by
one trisomy, so across-sample variance is largely screen noise and the correlation is uninterpretable.
Replace with a pooled analysis over the heterogeneous T1+T2 lines plus a descriptive
near-diploid-vs-trisomy contrast. **[local]** no clean aneuploidy direction exists (KIF18A: SS48 −1.34,
SS77 −1.18, SS6+7 −0.54, SS119+8 −1.60, SS51+7,22 −1.38), so do **not** pre-specify a directional one.

**15. §7.4/§8 — add the missing biological controls. SHOULD.** (i) *Positive, with DepMap support*:
GATA3/MYB/RUNX1/TAL1/IKZF1 (primary T), GATA2/MYB/GFI1B (CD34+ HSPC), SPI1/IRF8/CSF1R (iPSC-microglia) —
each has leukaemia-line support, so unlike pluripotency TFs these are fair tests. (ii) *Graded positive*:
genes unexpressed in the target type (TPM < 0.1) must be predicted and observed non-dependent; report FPR.
(iii) *Declared expected failures*: TCR-signalling genes (LCP2, ZAP70, LAT, CD3) are stimulation-dependent,
not cell-autonomous, and the enriched tail (CBLB, SOCS1, TCEB2, RASA2, TNFAIP3) lies outside DepMap's
dynamic range — report that tail separately with the model's predictive range, not as model failure.

**16. §2/§7 — state the missing-modality convention for primary cells. SHOULD.** "CN → 0 + indicator, no
damaging mutation" encodes *diploid, WT* — about right for primary T/HSPC, wrong for the RPE1-SS trisomies
and for hPSC with culture-acquired 20q11.21 or 12p gains (covering BCL2L1 and NANOG/GDF3). Capture reported
karyotype per screen; compare zero-fill against asserted-WT.

**17. §7.1/§11 — record Cas9 delivery as a §7.2 step-3 covariate. COULD.** Lentiviral / dox-inducible / RNP
/ Guide Swap. RNP gives a transient editing window, so mosaicism and in-frame repair attenuate effects with
no sustained selection for p53 loss — a different regime from lentiviral Cas9; the two should not share one
correction term.
