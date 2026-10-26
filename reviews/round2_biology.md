# Review round 2 — biological validity (DESIGN v1)

**[local]** = computed here from 26Q3 `CRISPRGeneEffect.csv`.

## 1. Verdict

**Accept with minor revisions.** All 17 round-1 items are in v1; §7.6, the §7.5 nulls, H2 and H3 are stronger
than I proposed. Six fixes remain, four in the new §7.5, where a definition and a gene list would inflate or
misfire the headline test.

**Self-correction:** my round-1 "POU5F1 < −0.5 in 2/1361 lines" was inferred, not counted. **[local]** the
true count is 18/1361 (1.32%); LIN28A 7, SALL4 11 — all still inside §4's *rare-selective* band, so item 1
stands. Correct the figure if quoted.

## 2. Remaining changes

**1. §7.5 — the screen-specific-dependency definition is too permissive. MUST.** v1 has "observed < −0.5 and
**not common-essential** in the cancer dev set", which excludes only genes essential in >95% of lines — so a
gene dependent in 60% of cancer lines qualifies. Restore the round-1 form: dependent in **≤5% of cancer dev
lines**. Otherwise the headline recall is dominated by near-core essentials that B0 already predicts — the
inflation the statistic exists to avoid.

**2. §7.5 H1 — drop USP7; move the effect-size clause to H1′. MUST.** **[local]** USP7 is < −0.5 in 974/1361
lines (71.6%): broadly essential, not TP53-WT-selective, so its transfer is guaranteed and would inflate H1.
Keep MDM2 (36.2%) and PPM1D (15.7%). Separately, H1's note "the effect size is 4–8× beyond the cancer range"
belongs to H1′: **[local]** MDM2 in RPE1-SS is −1.4 to −2.9, inside the cancer range (min −2.87).

**3. §7.2 step 3 — harmonisation can erase the hypothesis it enables. MUST.** All anchors are cancer lines
(HAP1 near-haploid, HCT116, DLD-1, A-375, K-562, HT-29, AML); none is p53-competent. A gene-wise
location/scale correction fitted on them and applied to RPE1 *TP53*-WT screens will read p53-competent
cutting toxicity as batch and shrink it. Exclude the p53/DDR gene set from the correction's fit **and**
application; report H1/H1′ before and after step 3.

**4. §7.5 — the CN-collateral test cannot detect the phenomenon as written. MUST.** RPE1-SS lines differ by
whole-chromosome *gains* (trisomy 7, 8, 7+22), whereas collateral/paralog-buffered lethality requires a
*loss*; and Chronos already CN-corrects, so residual dosage signal is near-zero by construction. Restate as
a null expectation across all non-transformed screens: SMARCA2, ARID1B, STAG1, ME2/ME3, ENO1 predicted ≈0
and observed ≈0 because the buffering paralog is intact. **[local]** well-posed — each is selective or
rare-selective (SMARCA2 3.4%, ARID1B 4.2%, STAG1 6.4%) and all sit at ≈0 in every RPE1-SS subline.

**5. §7.5 H4 — genes mis-assigned, and the dataset may not exist. MUST.** SPI1 belongs to the iPSC-microglia
CRISPRi stratum, which §7.1 makes descriptive-only, not to "CD34+/blood screens"; and the only genome-scale
CD34+ screen (Ting 2018) has an unresolved deposit, so H4 defaults to "inconclusive". Anchor H4 on the
obtainable stratum — Shifrut CD8+ T cells, with T-lineage genes (GATA3, MYB, RUNX1, IKZF1) — declare its
readout is proliferation under TCR stimulation, make the CD34+ arm (GATA2, MYB, GFI1B) conditional on data,
and require each H4 gene to be expressed in the matched profile; otherwise the model predicts ≈0 trivially.

**6. §5 model T — pre-commit the lineage mapping per screen. MUST.** "A fixed lineage mapping (or no
lineage)" is undefined and biologically loaded: mapping hPSC to Germ Cell routes B1 through the very EC lines
whose removal defines H3, confounding the EC leave-out. Commit the mapping in
`configs/transfer_criteria.json`, with **no lineage** as the primary variant.

**7. §2 — the p53/Rb composite was accepted but not updated. SHOULD.** It still reads "TP53 damaging or
hotspot status, plus a p53-target expression score". Add MDM2/MDM4 amplification (a main driver of MDM2
dependency alongside TP53-WT, so H1 needs it), TP53 CN loss and ATM/CHEK2/PPM1D; down-weight basal p53-target
expression. For Rb integrity use CDKN2A deletion/CN, not expression — deletion is the common lesion.

**8. §3 — assign tiers functionally, not from annotation. SHOULD.** A line is p53-competent if its observed
TP53 KO effect > +0.5 and MDM2 < −1. **[local]** all five RPE1-SS pass (TP53 +1.6 to +3.2); HA1E fails (TP53
−0.47, MDM2 −0.07), consistent with SV40 large T. This guards against annotation error and Cas9-driven p53
loss during screening. Wording: A1 calls all five RPE1-SS "aneuploid", but SS48/SS77 are the near-diploid
controls; and Hs 729 sits in a tier titled "virally transformed lymphoid/myeloid lines" — give it its own
unverified-provenance bucket.

**9. §7.1 — restore the neuron gene-set split. SHOULD.** v1 reduces it to "divergence is expected and
declared in advance". State the two sets (proliferation-coupled: replication, mitosis, spindle, nucleotide
synthesis; maintenance: translation, proteostasis, mitochondrial, lysosomal) and the direction — the only
a-priori-directional test available in a descriptive-only stratum.
