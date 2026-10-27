# G3 kernel track — report columns, controls, M1b tuning, M2 (D14 base), secondary splits

<!-- notes -->
## Summary and interpretation (written after the runs; numbers come from the tables below)

**Selection (DESIGN §4: paired ΔP1 > 0.005 and > 1 fold SE vs `m1_E_resB0`).** No kernel-track
candidate replaces the headline. **`m1_E_resB0` (run 20261026-222117, P1 variable 0.3081) stays
the current best.** The largest gains are M2 with paralog priors at +0.0034 (SE 0.0004, NB SE 0.0006;
5/5 folds) and tuned M1b on the full set at +0.0015. Both are consistent across folds but are below δ.

**1. Report columns (D16, D17).** Within-library and within-lineage P1 use the reviewer's
definitions and reproduce the review's numbers exactly for the frozen re-runs:
- m1_E_resB0: within-library 0.2790, within-lineage 0.2449;
- B2: 0.2256 / 0.1769;
- B1: 0.1709 / 0.

A third diagnostic column, P1 within lineage × library cells, was added after the lineage×library
control: 0.199 for the headline, against −0.007 for that control. That is line-level skill that
neither lineage nor library group can explain.

The Nadeau–Bengio SE is about 1.5× the fold SE, as the reviewer predicted. No conclusion changes:
every M1b/M2 gain is many NB SEs from 0 and far below δ.

**2. Negative controls (expectations were fixed in the configs before running).**
- *Global shuffle on the B0 path.* P1 0.0026, which passes (|P1| < 0.01).
  - My stated side-expectation ("λ goes to the upper grid edge") was **wrong**: λ went to the
    *lower* edge (100).
  - Reason: on the B0 path the closed-form LOO-r criterion tends to −1 as λ → ∞. The LOO
    prediction becomes the mean of the other lines, the same artefact as D3. So with no signal
    the criterion prefers the smallest λ. A toy run confirms it: pure noise gives LOO r
    −0.01 → −0.98 across the grid.
  - This is not leakage: outer-fold P1 is ≈ 0. But it means the LOO-r criterion is biased toward
    small λ for weak-signal genes on the headline path. See "decisions" below.
- *Expression deranged within lineage.*
  - Within-lineage P1 −0.004, which passes. No line-specific skill survives.
  - P1 0.1265 is *below* the expected band [B1 − 0.04, B1 + 0.01], so that check is a FAIL in the
    conservative direction. The ridge recovers less lineage information from donor profiles than
    the explicit shrunk lineage mean does.
- *Expression deranged within lineage × library group.*
  - P1 0.1755, within-library P1 0.133 ≤ B1's 0.171, and within-lineage×library P1 −0.007.
  - Within-lineage P1 is **0.071**, above the expected < 0.03, so that check is a FAIL.
  - So the library group (dev NaN-pattern group) explains a real share of within-lineage
    variation in gene effects, and features matched on lineage×library can carry it. It does not
    explain the headline: the headline keeps 0.199 within lineage×library cells.
  - It is still a caution. Up to ~0.07 of the headline's within-lineage P1 (0.245) could be
    library-group structure, which is screen-technical or confounded with it. This control's
    strata use the held-out lines' library groups (a CRISPR NaN mask) by design.
- The G2 control (global shuffle, B1 path) is unchanged: P1 = B1.
- The independent leakage review noted that the global shuffle maps about 217 dev lines onto
  test-line *feature* rows. Those are features, not targets, and the run is a null, so it cannot
  inflate anything. The stratified controls permute dev lines only.

**3. M1b block-weight tuning (D15).**
- Method: coordinate descent on the training-fold leave-patient-out criterion; largest NaN-pattern
  cluster, training-fold variable genes, grid θ ∈ {0, 1/64, 1/16, 1/4, 1} × the default 1/√p_b
  kernel, then ×½/×2 refinement; at most 25 eigh per outer fold.
- Tuned weights are stable across folds:
  - C 1/32 (E+C) or 1/64 (larger sets);
  - M 1/16;
  - H and F 1/64 or 0;
  - Meta 1/4;
  - composites 1/64–1/16.
- So the default weights over-weighted genomics by 16–64×.
- With tuned weights every set beats E-only, by +0.0004 (E+C) to +0.0015 (full). All are 5/5 folds
  with Holm p < 0.01 over the 6 declared contrasts. The C common-subset rule gives the same.
- **Genomics carries real but tiny signal (≤ 0.0015 P1), far below δ.** The G2 "every modality
  hurts" result came from the weights. The D15 condition for a negative claim is now met:
  "genomics does not add a selectable gain over expression".
- Caveat: sets with ≥ 3 free modalities exhausted the 25-eigh budget in one coordinate sweep, with
  no refinement, so their optimum is coarse.

**4. M2 on the D14 base.**
- Setup: B0-residual M1; local stage only for training-fold variable genes, others pass M1
  through; nested local-λ CV (review finding 15).
- Gains: none +0.0019, paralogs +0.0034, partners +0.0025. Paralogs vs none +0.0015 (Holm p 0.009);
  partners vs paralogs −0.0009 (n.s.).
- Pass-through fixed the G2 loss on all genes: P1-all is 0.1551 vs 0.1547 for the headline (G2
  M2: 0.141–0.147).
- The paralog gain is larger within lineage × library cells: +0.0104 (5/5). It is the one rung
  whose within-group gain exceeds its raw gain, i.e. it adds line-specific (genotype/paralog)
  information. It still does not clear δ on the primary metric.
- The local λ moved from 10 (none) to 31.6 (paralogs, partners) with nesting.

**5. EC leave-out and H3 no-EC variant (T trained without the dev germ-cell/EC lines).**
- Only **4 of the DESIGN §3 "7 germ-cell/EC lines" are dev lines**: NCCIT, NTERA-2 cl.D1, PA-1
  and SUSA. GCT27 is a test line; 1411H and SCH are not in the screened population.
  configs/transfer_criteria.json says "the 7 germ-cell/EC dev lines". The PI should correct the
  wording before G4.
- Removing them changes nothing elsewhere: ΔP1 ≤ 0.0002 for every model.
- Per-line skill of the EC lines barely moves.
- **The pluripotency genes are not predicted even WITH EC lines in training.**
  - Observed POU5F1 is −0.88 in NCCIT and NTERA-2. The predictions of m1_E_resB0 and T are
    −0.09 to −0.10 with EC lines in training and −0.04 to −0.07 without.
  - PRDM14 (observed −0.59 to −0.64) is predicted at about −0.15.
  - With only 2–4 EC lines, the cross-line ridge shrinks these rare dependencies almost to the
    gene mean.
  - For H3 (T must predict < −0.5 for ≥ 50% of the hPSC-essential pluripotency genes), this
    strongly suggests T will fail H3 whether or not the EC lines are included. The "narrow EC
    transfer" contrast will probably be uninformative.
  - This was observed on dev lines only. It is noted now so that nobody tunes for it after the
    fact.

**6. Leave-one-lineage-out (18 dev lineages with ≥ 20 lines; DESIGN §4 threshold).**
- m1_E_resB0 keeps skill on unseen lineages:
  - mean per-lineage P1 0.203, P2b 0.198;
  - T: 0.205 / 0.199;
  - B2: 0.120 / 0.142;
  - B1 = B0 by construction, so 0.
- M1 beats B2 in **18/18** lineages (mean +0.083).
- Weakest lineages: haematopoietic (Lymphoid 0.094, Myeloid 0.107), Head and Neck (0.166) and
  Bone (0.162).
- Strongest: Esophagus/Stomach, Pancreas and Uterus (0.26–0.28).
- For scale, the same model's within-lineage P1 under 5-fold CV is 0.245. That figure pools
  lineages and has same-lineage lines in training, so it is not a strict comparator. An unseen
  lineage costs roughly 0.04 P1, not the whole skill.
- This is the most transfer-relevant number on dev: set A and set B cells are an "unseen lineage"
  in the same sense.
- The first LOLO attempt failed after fitting all 28 folds: the frozen `evaluate` cannot score a
  fold in which no gene has ≥ 20 observed lines. Fix: only lineages with ≥ 20 dev lines are held
  out (`lolo_min_lines`). Smaller lineages stay in every training set, so each scored fit is the
  same fit as in full LOLO. Predictions are now saved before evaluation.
- LOLO fit times (sum over 18 folds, 6 threads, excluding the inline reference baselines): B1 ~630 s, B2 ~970 s, M1 ~1970 s, T ~2900 s.

**7. Runtimes (indicative only).** All runs used 6 threads on a machine shared with the nonlinear
track (LESSONS 14). Totals for 5 folds:
| run | time |
|---|---|
| M1-E | ~330 s |
| M1b tuned E+C | ~600 s |
| M1b tuned, full set | ~1580 s (tuning ≈ 200–250 s per fold) |
| M2, first prior variant | ~5060 s (30 inner/nested M1 fits per outer fold) |
| M2, later variants | ~900 s (shared inner-OOF cache) |

**8. Does anything look too good?** No.
- Every new gain is small, uniform across folds, and comes with matching within-library and
  within-lineage gains.
- The independent leakage review (fresh agent, did not write the code) of the M1b tuning, M2
  nesting and the run.py exclusion/LOLO code found no leakage into dev-CV results. Its fixes were
  applied (commit 959d96f):
  - saved fold models record the actual training IDs;
  - exact, versioned M2 inner-OOF cache key;
  - atomic cache writes;
  - a budget guard;
  - tests for exclusion and LOLO.
- The nonlinear track's M3 reads the same inner-OOF cache through the same `LocalStageM2.inner_oof`
  and exact key. That explains one cache hit per fold in M2 and is legitimate.
<!-- /notes -->

Generated by `scripts/g3_kernel_report.py` from `results/registry.csv` (latest ok run per name).
Dev lines only, 5 outer dev folds unless stated. P1 = mean per-gene Pearson r on variable
genes within each fold, then averaged (± fold SD). ΔP1 = paired fold-wise difference; SE =
SD/√5; **NB SE** = Nadeau–Bengio corrected SE (disclosure only; the selection rule uses the
fold SE). Within-library / within-lineage P1: truth and prediction centred per gene within
library group (configs/library_groups.csv) or lineage inside each held-out fold (definitions:
reviews/g2_review/recompute_p1.py `group_centre`; lineage partialling of neg_genes.py).
Within lineage×library: the same centring within lineage × library-group cells (added in G3
after the lineage×library control; diagnostic only).
P2b = normalised AP on selective genes with held-out prevalence ≤ 0.5. Fit s = sum over folds
(incl. inner tuning), measured with 6 threads on a machine shared with the nonlinear track:
indicative only (LESSONS 14).

## Results table (5-fold dev CV)

| run | P1 var ± SD | ΔP1 vs m1_E_resB0 (SE / NB SE; folds +) | ΔP1 vs B1 | within-library P1 | within-lineage P1 | within lineage×library P1 | P2b nAP | fit s | λ (other/var) |
|---|---|---|---|---|---|---|---|---|---|
| b1 | 0.1837 ± 0.0131 | -0.1244 (0.0048 / 0.0072; 0/5) | — | 0.1709 | 0.0000 | 0.0000 | 0.0867 | 105 | — |
| b2 | 0.2527 ± 0.0087 | -0.0554 (0.0023 / 0.0035; 0/5) | +0.0689 | 0.2256 | 0.1769 | 0.1271 | 0.1524 | 127 | — |
| r1_E | 0.3073 ± 0.0050 | -0.0008 (0.0008 / 0.0011; 1/5) | +0.1236 | 0.2770 | 0.2431 | 0.1953 | 0.1898 | 356 | {"0": [46415.888336127726]} |
| m1_E | 0.2997 ± 0.0067 | -0.0084 (0.0017 / 0.0025; 0/5) | +0.1159 | 0.2717 | 0.2367 | 0.1915 | 0.1845 | 370 | {"0": [21544.346900318822], "1": [21544.346900318822]} |
| t_E_status | 0.3082 ± 0.0040 | +0.0001 (0.0004 / 0.0006; 3/5) | +0.1245 | 0.2793 | 0.2455 | 0.2006 | 0.1914 | 378 | {"0": [46415.888336127726], "1": [21544.346900318822]} |
| m1_E_resB0 | 0.3081 ± 0.0047 | — | +0.1244 | 0.2790 | 0.2449 | 0.1988 | 0.1908 | 327 | {"0": [46415.888336127726], "1": [21544.346900318822]} |
| ctl_resB0_permuted | 0.0026 ± 0.0050 | -0.3055 (0.0031 / 0.0046; 0/5) | -0.1811 | 0.0031 | 0.0010 | -0.0005 | 0.0256 | 315 | {"0": [100.0], "1": [100.0, 215.44346900318845]} |
| ctl_resB0_perm_lineage | 0.1265 ± 0.0033 | -0.1816 (0.0016 / 0.0024; 0/5) | -0.0572 | 0.1198 | -0.0040 | -0.0016 | 0.0792 | 525 | {"0": [46415.888336127726, 100000.0], "1": [100000.0]} |
| ctl_resB0_perm_linlib | 0.1755 ± 0.0066 | -0.1326 (0.0031 / 0.0047; 0/5) | -0.0083 | 0.1327 | 0.0708 | -0.0065 | 0.1001 | 535 | {"0": [100000.0], "1": [100000.0]} |
| m1_E_permuted | 0.1838 ± 0.0131 | -0.1243 (0.0048 / 0.0072; 0/5) | +0.0000 | — | — | — | 0.1050 | 503 | {"0": [10000000.0], "1": [10000000.0]} |
| b1_noEC | 0.1836 ± 0.0129 | -0.1245 (0.0047 / 0.0070; 0/5) | -0.0002 | 0.1708 | 0.0000 | 0.0000 | 0.0863 | 186 | — |
| b2_noEC | 0.2526 ± 0.0087 | -0.0555 (0.0024 / 0.0036; 0/5) | +0.0689 | 0.2256 | 0.1772 | 0.1278 | 0.1525 | 213 | — |
| m1_E_resB0_noEC | 0.3080 ± 0.0046 | -0.0001 (0.0001 / 0.0001; 0/5) | +0.1242 | 0.2788 | 0.2448 | 0.1987 | 0.1907 | 508 | {"0": [46415.888336127726], "1": [21544.346900318822]} |
| m1bt_E_C | 0.3085 ± 0.0046 | +0.0004 (0.0001 / 0.0001; 5/5) | +0.1248 | 0.2795 | 0.2456 | 0.2001 | 0.1915 | 596 | {"0": [46415.888336127726], "1": [21544.346900318822]} |
| m1bt_E_C_M | 0.3090 ± 0.0046 | +0.0009 (0.0001 / 0.0001; 5/5) | +0.1252 | 0.2799 | 0.2461 | 0.2004 | 0.1918 | 1065 | {"0": [46415.888336127726], "1": [21544.346900318822]} |
| m1bt_E_C_M_H | 0.3089 ± 0.0045 | +0.0008 (0.0001 / 0.0002; 5/5) | +0.1252 | 0.2799 | 0.2460 | 0.2003 | 0.1923 | 1136 | {"0": [46415.888336127726], "1": [21544.346900318822]} |
| m1bt_E_C_M_H_F | 0.3089 ± 0.0045 | +0.0008 (0.0001 / 0.0002; 5/5) | +0.1251 | 0.2798 | 0.2458 | 0.1999 | 0.1920 | 1386 | {"0": [46415.888336127726], "1": [21544.346900318822]} |
| m1bt_E_C_M_H_F_Meta | 0.3091 ± 0.0048 | +0.0010 (0.0001 / 0.0001; 5/5) | +0.1253 | 0.2801 | 0.2457 | 0.1997 | 0.1920 | 1345 | {"0": [46415.888336127726], "1": [21544.346900318822]} |
| m1bt_full | 0.3096 ± 0.0045 | +0.0015 (0.0002 / 0.0003; 5/5) | +0.1259 | 0.2807 | 0.2464 | 0.2009 | 0.1928 | 1576 | {"0": [46415.888336127726], "1": [21544.346900318822]} |
| m2n_E_none | 0.3100 ± 0.0050 | +0.0019 (0.0002 / 0.0004; 5/5) | +0.1263 | 0.2819 | 0.2493 | 0.2057 | 0.1922 | 5063 | {"0": [46415.888336127726], "1": [21544.346900318822]} |
| m2n_E_paralogs | 0.3115 ± 0.0048 | +0.0034 (0.0004 / 0.0006; 5/5) | +0.1277 | 0.2835 | 0.2520 | 0.2092 | 0.1946 | 876 | {"0": [46415.888336127726], "1": [21544.346900318822]} |
| m2n_E_partners | 0.3106 ± 0.0042 | +0.0025 (0.0005 / 0.0007; 5/5) | +0.1268 | 0.2825 | 0.2513 | 0.2085 | 0.1945 | 934 | {"0": [46415.888336127726], "1": [21544.346900318822]} |
| t_E_status_noEC | 0.3081 ± 0.0040 | -0.0000 (0.0004 / 0.0006; 2/5) | +0.1243 | 0.2792 | 0.2454 | 0.2005 | 0.1912 | 377 | {"0": [46415.888336127726], "1": [21544.346900318822]} |

Run IDs: b1 `20261026-220921_b1_0e2a7a27`; b2 `20261026-221456_b2_289f5170`; r1_E `20261026-223037_r1_E_1df71977`; m1_E `20261026-225039_m1_E_e66eff5f`; t_E_status `20261026-224125_t_E_status_24b46e2a`; m1_E_resB0 `20261026-222117_m1_E_resB0_91e52f9d`; ctl_resB0_permuted `20261026-231349_ctl_resB0_permuted_28e789f7`; ctl_resB0_perm_lineage `20261026-232244_ctl_resB0_perm_lineage_e6fcbf0b`; ctl_resB0_perm_linlib `20261026-233916_ctl_resB0_perm_linlib_f3b0b9bb`; m1_E_permuted `20261026-200800_m1_E_permuted_bb3bec92`; b1_noEC `20261027-024853_b1_noEC_19d271a7`; b2_noEC `20261027-030055_b2_noEC_bd320770`; m1_E_resB0_noEC `20261027-023409_m1_E_resB0_noEC_1d8386a5`; m1bt_E_C `20261026-235133_m1bt_E_C_39d641a9`; m1bt_E_C_M `20261027-005238_m1bt_E_C_M_a09d90b1`; m1bt_E_C_M_H `20261027-011448_m1bt_E_C_M_H_3d5feb4e`; m1bt_E_C_M_H_F `20261027-013815_m1bt_E_C_M_H_F_0d699efd`; m1bt_E_C_M_H_F_Meta `20261027-020528_m1bt_E_C_M_H_F_Meta_1d50afff`; m1bt_full `20261027-002223_m1bt_full_48508032`; m2n_E_none `20261027-031420_m2n_E_none_8263b78a`; m2n_E_paralogs `20261027-044450_m2n_E_paralogs_62ce1575`; m2n_E_partners `20261027-050532_m2n_E_partners_ea101d11`; t_E_status_noEC `20261027-000415_t_E_status_noEC_fa0181c1`

## Negative controls vs expected outcomes (stated in the configs before running)

| control | P1 var | ΔP1 vs B1 | within-lineage | within-library | within lin×lib | λ | checks | pass |
|---|---|---|---|---|---|---|---|---|
| ctl_resB0_permuted | 0.0026 | -0.1811 | 0.0010 | 0.0031 | -0.0005 | {"0": [100.0], "1": [100.0, 215.44346900318845]} | abs(P1) < 0.01: ok | yes |
| ctl_resB0_perm_lineage | 0.1265 | -0.0572 | -0.0040 | 0.1198 | -0.0016 | {"0": [46415.888336127726, 100000.0], "1": [100000.0]} | B1-0.04 <= P1 <= B1+0.01: FAIL; abs(within-lineage P1) < 0.01: ok; headline - P1 > 0.08: ok | NO |
| ctl_resB0_perm_linlib | 0.1755 | -0.0083 | 0.0708 | 0.1327 | -0.0065 | {"0": [100000.0], "1": [100000.0]} | B1-0.04 <= P1 <= B1+0.03: ok; within-lineage P1 < 0.03: FAIL; within-library P1 <= B1 within-library + 0.01: ok; headline - P1 > 0.08: ok | NO |
| m1_E_permuted | 0.1838 | +0.0000 | — | — | — | {"0": [10000000.0], "1": [10000000.0]} | abs(P1 - B1) < 0.01: ok | yes |

## M1b block-weight tuning (B0 residual; training LOO coordinate descent)

| contrast | ΔP1 | SE | NB SE | folds + | p | p Holm | helps |
|---|---|---|---|---|---|---|---|
| m1bt_E_C vs m1_E_resB0 | +0.0004 | 0.0001 | 0.0001 | 5/5 | 0.00319 | 0.00575 | yes |
| m1bt_E_C_M vs m1_E_resB0 | +0.0009 | 0.0001 | 0.0001 | 5/5 | 0.000968 | 0.00484 | yes |
| m1bt_E_C_M_H vs m1_E_resB0 | +0.0008 | 0.0001 | 0.0002 | 5/5 | 0.0016 | 0.00575 | yes |
| m1bt_E_C_M_H_F vs m1_E_resB0 | +0.0008 | 0.0001 | 0.0002 | 5/5 | 0.0019 | 0.00575 | yes |
| m1bt_E_C_M_H_F_Meta vs m1_E_resB0 | +0.0010 | 0.0001 | 0.0001 | 5/5 | 0.000242 | 0.00145 | yes |
| m1bt_full vs m1_E_resB0 | +0.0015 | 0.0002 | 0.0003 | 5/5 | 0.00144 | 0.00575 | yes |
- C common-subset rule: m1bt_E_C vs m1_E_resB0 (C-covered lines only): ΔP1 +0.0005 ± 0.0001 (5/5 folds, p 0.00798)
- C common-subset rule: m1bt_E_C_M vs m1_E_resB0 (C-covered lines only): ΔP1 +0.0009 ± 0.0001 (5/5 folds, p 0.000885)
- C common-subset rule: m1bt_E_C_M_H vs m1_E_resB0 (C-covered lines only): ΔP1 +0.0009 ± 0.0001 (5/5 folds, p 0.00169)
- C common-subset rule: m1bt_E_C_M_H_F vs m1_E_resB0 (C-covered lines only): ΔP1 +0.0008 ± 0.0001 (5/5 folds, p 0.00105)
- C common-subset rule: m1bt_E_C_M_H_F_Meta vs m1_E_resB0 (C-covered lines only): ΔP1 +0.0010 ± 0.0001 (5/5 folds, p 0.000153)
- C common-subset rule: m1bt_full vs m1_E_resB0 (C-covered lines only): ΔP1 +0.0015 ± 0.0002 (5/5 folds, p 0.00238)

Tuned block weights per outer fold (multiplier θ_b of the default 1/√p_b kernel; E fixed at 1):

- m1bt_E_C fold 0: E 1, C 0.03125 (7 eigh)
- m1bt_E_C fold 1: E 1, C 0.03125 (7 eigh)
- m1bt_E_C fold 2: E 1, C 0.03125 (7 eigh)
- m1bt_E_C fold 3: E 1, C 0.03125 (7 eigh)
- m1bt_E_C fold 4: E 1, C 0.03125 (7 eigh)
- m1bt_E_C_M fold 0: E 1, C 0.03125, M 0.0625 (17 eigh)
- m1bt_E_C_M fold 1: E 1, C 0.03125, M 0.0625 (17 eigh)
- m1bt_E_C_M fold 2: E 1, C 0.03125, M 0.0625 (17 eigh)
- m1bt_E_C_M fold 3: E 1, C 0.03125, M 0.0625 (17 eigh)
- m1bt_E_C_M fold 4: E 1, C 0.03125, M 0.0625 (17 eigh)
- m1bt_E_C_M_H fold 0: E 1, C 0.03125, M 0.0625, H 0.015625 (25 eigh)
- m1bt_E_C_M_H fold 1: E 1, C 0.03125, M 0.0625, H 0.015625 (25 eigh)
- m1bt_E_C_M_H fold 2: E 1, C 0.03125, M 0.0625, H 0.015625 (25 eigh)
- m1bt_E_C_M_H fold 3: E 1, C 0.03125, M 0.0625, H 0.015625 (25 eigh)
- m1bt_E_C_M_H fold 4: E 1, C 0.03125, M 0.0625, H 0 (21 eigh)
- m1bt_E_C_M_H_F fold 0: E 1, C 0.015625, M 0.0625, H 0.015625, F 0.015625 (25 eigh)
- m1bt_E_C_M_H_F fold 1: E 1, C 0.015625, M 0.0625, H 0.015625, F 0 (25 eigh)
- m1bt_E_C_M_H_F fold 2: E 1, C 0.015625, M 0.0625, H 0.015625, F 0.015625 (25 eigh)
- m1bt_E_C_M_H_F fold 3: E 1, C 0.015625, M 0.0625, H 0.015625, F 0.015625 (25 eigh)
- m1bt_E_C_M_H_F fold 4: E 1, C 0.03125, M 0.0625, H 0, F 0 (25 eigh)
- m1bt_E_C_M_H_F_Meta fold 0: E 1, C 0.015625, M 0.0625, H 0.015625, F 0.015625, Meta 0.25 (25 eigh)
- m1bt_E_C_M_H_F_Meta fold 1: E 1, C 0.015625, M 0.0625, H 0.015625, F 0, Meta 0.25 (25 eigh)
- m1bt_E_C_M_H_F_Meta fold 2: E 1, C 0.015625, M 0.0625, H 0.015625, F 0.015625, Meta 0.25 (25 eigh)
- m1bt_E_C_M_H_F_Meta fold 3: E 1, C 0.015625, M 0.0625, H 0.015625, F 0.015625, Meta 0.25 (25 eigh)
- m1bt_E_C_M_H_F_Meta fold 4: E 1, C 0.015625, M 0.0625, H 0, F 0, Meta 0.25 (25 eigh)
- m1bt_full fold 0: E 1, C 0.015625, M 0.0625, H 0.015625, F 0.015625, Meta 0.25, composites 0.015625 (25 eigh)
- m1bt_full fold 1: E 1, C 0.015625, M 0.0625, H 0.015625, F 0, Meta 0.25, composites 0.0625 (25 eigh)
- m1bt_full fold 2: E 1, C 0.015625, M 0.0625, H 0.015625, F 0.015625, Meta 0.25, composites 0.015625 (25 eigh)
- m1bt_full fold 3: E 1, C 0.015625, M 0.0625, H 0.015625, F 0.015625, Meta 0.25, composites 0.0625 (25 eigh)
- m1bt_full fold 4: E 1, C 0.015625, M 0.0625, H 0, F 0, Meta 0.25, composites 0.015625 (25 eigh)

## M2 (D14 base, pass-through for non-variable genes, nested local-λ CV)

- M2 paralogs vs none: ΔP1 +0.0015 (SE 0.0003, NB 0.0004; Holm p 0.00939)
- M2 partners vs paralogs: ΔP1 -0.0009 (SE 0.0003, NB 0.0005; Holm p 0.0555)
- m2n_E_none local λ per fold: {"0": {"1": 10.0}, "1": {"1": 10.0}, "2": {"1": 10.0}, "3": {"1": 3.1622776601683795}, "4": {"1": 3.1622776601683795}}
- m2n_E_paralogs local λ per fold: {"0": {"1": 31.622776601683793}, "1": {"1": 31.622776601683793}, "2": {"1": 31.622776601683793}, "3": {"1": 31.622776601683793}, "4": {"1": 31.622776601683793}}
- m2n_E_partners local λ per fold: {"0": {"1": 31.622776601683793}, "1": {"1": 31.622776601683793}, "2": {"1": 31.622776601683793}, "3": {"1": 31.622776601683793}, "4": {"1": 31.622776601683793}}

## Selection rule vs m1_E_resB0 (DESIGN §4: ΔP1 > 0.005 and > 1 fold SE)

| candidate | ΔP1 | SE | NB SE | replaces |
|---|---|---|---|---|
| m1bt_E_C | +0.0004 | 0.0001 | 0.0001 | no |
| m1bt_E_C_M | +0.0009 | 0.0001 | 0.0001 | no |
| m1bt_E_C_M_H | +0.0008 | 0.0001 | 0.0002 | no |
| m1bt_E_C_M_H_F | +0.0008 | 0.0001 | 0.0002 | no |
| m1bt_E_C_M_H_F_Meta | +0.0010 | 0.0001 | 0.0001 | no |
| m1bt_full | +0.0015 | 0.0002 | 0.0003 | no |
| m2n_E_none | +0.0019 | 0.0002 | 0.0004 | no |
| m2n_E_paralogs | +0.0034 | 0.0004 | 0.0006 | no |
| m2n_E_partners | +0.0025 | 0.0005 | 0.0007 | no |

**Selected model after the G3 kernel track: `m1_E_resB0`.**

## Secondary split: leave-one-lineage-out (dev lines; lineages with ≥ 20 lines scored)

| lineage | n | b1_lolo P1 / P2b | b2_lolo P1 / P2b | m1_E_resB0_lolo P1 / P2b | t_E_status_lolo P1 / P2b |
|---|---|---|---|---|---|
| CNS/Brain | 102 | 0.0000 / 0.0000 | 0.1085 / 0.0904 | 0.1707 / 0.1322 | 0.1766 / 0.1384 |
| Lung | 102 | 0.0000 / 0.0000 | 0.1328 / 0.1057 | 0.2275 / 0.1616 | 0.2277 / 0.1603 |
| Lymphoid | 78 | 0.0000 / 0.0000 | 0.0108 / 0.0557 | 0.0940 / 0.0983 | 0.0975 / 0.0995 |
| Esophagus/Stomach | 72 | 0.0000 / 0.0000 | 0.2226 / 0.1584 | 0.2844 / 0.1997 | 0.2836 / 0.2008 |
| Bowel | 69 | 0.0000 / 0.0000 | 0.1175 / 0.1106 | 0.2325 / 0.1695 | 0.2355 / 0.1717 |
| Pancreas | 59 | 0.0000 / 0.0000 | 0.2027 / 0.1528 | 0.2824 / 0.2035 | 0.2827 / 0.2043 |
| Skin | 59 | 0.0000 / 0.0000 | 0.0827 / 0.1062 | 0.1937 / 0.1678 | 0.1980 / 0.1716 |
| Breast | 54 | 0.0000 / 0.0000 | 0.1770 / 0.1510 | 0.2520 / 0.1912 | 0.2558 / 0.1927 |
| Head and Neck | 51 | 0.0000 / 0.0000 | 0.0749 / 0.1142 | 0.1662 / 0.1853 | 0.1653 / 0.1834 |
| Soft Tissue | 42 | 0.0000 / 0.0000 | 0.1359 / 0.1613 | 0.1881 / 0.1953 | 0.1874 / 0.1915 |
| Ovary/Fallopian Tube | 41 | 0.0000 / 0.0000 | 0.1375 / 0.1740 | 0.2429 / 0.2387 | 0.2436 / 0.2384 |
| Myeloid | 38 | 0.0000 / 0.0000 | 0.0340 / 0.1087 | 0.1067 / 0.1607 | 0.1091 / 0.1631 |
| Bone | 36 | 0.0000 / 0.0000 | 0.0875 / 0.1386 | 0.1616 / 0.2067 | 0.1577 / 0.1998 |
| Uterus | 31 | 0.0000 / 0.0000 | 0.1643 / 0.1767 | 0.2624 / 0.2390 | 0.2648 / 0.2435 |
| Kidney | 28 | 0.0000 / 0.0000 | 0.0958 / 0.1774 | 0.1749 / 0.2272 | 0.1810 / 0.2327 |
| Peripheral Nervous System | 28 | 0.0000 / 0.0000 | 0.0971 / 0.1500 | 0.1693 / 0.2287 | 0.1737 / 0.2295 |
| Bladder/Urinary Tract | 25 | 0.0000 / 0.0000 | 0.1239 / 0.1945 | 0.1924 / 0.2571 | 0.1943 / 0.2587 |
| Biliary Tract | 24 | 0.0000 / 0.0000 | 0.1586 / 0.2306 | 0.2519 / 0.3030 | 0.2533 / 0.3083 |
| **mean (scored)** | | 0.0000 / 0.0000 | 0.1202 / 0.1421 | 0.2030 / 0.1981 | 0.2049 / 0.1993 |
- b1_lolo - b2_lolo: mean ΔP1 -0.1202 over 18 lineages, positive in 0
- m1_E_resB0_lolo - b2_lolo: mean ΔP1 +0.0828 over 18 lineages, positive in 18
- t_E_status_lolo - b2_lolo: mean ΔP1 +0.0846 over 18 lineages, positive in 18
- Run IDs: b1_lolo `20261027-105225_b1_lolo_7c8331b9`; b2_lolo `20261027-101234_b2_lolo_a08ae146`; m1_E_resB0_lolo `20261027-081651_m1_E_resB0_lolo_f17ddbd4`; t_E_status_lolo `20261027-091056_t_E_status_lolo_bcef307c`

## Secondary split: EC leave-out (4 dev germ-cell/EC lines removed from training)

| no-EC run | P1 var | P1 var with EC | Δ (SE) |
|---|---|---|---|
| b1_noEC | 0.1836 | 0.1837 | -0.0002 (0.0002) |
| b2_noEC | 0.2526 | 0.2527 | -0.0000 (0.0001) |
| m1_E_resB0_noEC | 0.3080 | 0.3081 | -0.0001 (0.0001) |
| t_E_status_noEC | 0.3081 | 0.3082 | -0.0001 (0.0001) |

Per EC line (all eval genes): per-line r and ΔAP vs B0, with EC lines in training → without.

| line | lineage | b1 | b2 | m1_E_resB0 | t_E_status |
|---|---|---|---|---|---|
| NCCIT (ACH-001578) | Testis | r 0.912→0.913; ΔAP +0.000→+0.000 | r 0.908→0.906; ΔAP -0.003→-0.003 | r 0.912→0.913; ΔAP +0.003→+0.003 | r 0.912→0.912; ΔAP +0.002→+0.003 |
| NTERA-2 cl.D1 (ACH-002291) | Testis | r 0.887→0.886; ΔAP +0.002→+0.000 | r 0.882→0.881; ΔAP -0.004→-0.006 | r 0.889→0.889; ΔAP +0.006→+0.006 | r 0.889→0.890; ΔAP +0.007→+0.007 |
| PA-1 [PA1] (ACH-001374) | Ovary/Fallopian Tube | r 0.882→0.883; ΔAP -0.005→-0.004 | r 0.884→0.885; ΔAP -0.003→-0.003 | r 0.885→0.885; ΔAP -0.002→-0.002 | r 0.884→0.885; ΔAP -0.003→-0.003 |
| SUSA (ACH-001668) | Testis | r 0.906→0.906; ΔAP +0.001→+0.000 | r 0.902→0.901; ΔAP -0.003→-0.005 | r 0.901→0.901; ΔAP -0.001→-0.002 | r 0.901→0.902; ΔAP -0.001→-0.001 |

H3 pluripotency genes (observed effect; predicted with EC → without EC):

| gene | line | observed | b1 | b2 | m1_E_resB0 | t_E_status |
|---|---|---|---|---|---|---|
| POU5F1 | NCCIT | -0.8822 | -0.04 → -0.01 | -0.05 → -0.01 | -0.09 → -0.04 | -0.08 → -0.04 |
| POU5F1 | NTERA-2 cl.D1 | -0.8870 | -0.04 → -0.01 | -0.08 → -0.04 | -0.09 → -0.05 | -0.10 → -0.07 |
| POU5F1 | PA-1 [PA1] | -0.3628 | -0.01 → -0.01 | -0.13 → -0.07 | -0.18 → -0.06 | -0.19 → -0.07 |
| POU5F1 | SUSA | -0.2760 | -0.10 → -0.02 | -0.12 → -0.01 | -0.27 → -0.07 | -0.28 → -0.08 |
| NANOG | NCCIT | 0.1270 | +0.22 → +0.22 | +0.22 → +0.20 | +0.16 → +0.16 | +0.15 → +0.15 |
| NANOG | NTERA-2 cl.D1 | — | +0.22 → +0.22 | +0.21 → +0.24 | +0.19 → +0.19 | +0.22 → +0.23 |
| NANOG | PA-1 [PA1] | — | +0.22 → +0.22 | +0.22 → +0.22 | +0.17 → +0.18 | +0.18 → +0.18 |
| NANOG | SUSA | — | +0.21 → +0.21 | +0.27 → +0.30 | +0.19 → +0.20 | +0.18 → +0.19 |
| SOX2 | NCCIT | -0.4294 | -0.17 → -0.15 | -0.13 → -0.10 | -0.35 → -0.28 | -0.37 → -0.29 |
| SOX2 | NTERA-2 cl.D1 | -0.1899 | -0.17 → -0.15 | -0.27 → -0.23 | -0.28 → -0.22 | -0.26 → -0.20 |
| SOX2 | PA-1 [PA1] | -0.5681 | -0.15 → -0.15 | -0.27 → -0.27 | -0.24 → -0.24 | -0.22 → -0.22 |
| SOX2 | SUSA | -0.3500 | -0.16 → -0.15 | -0.22 → -0.17 | -0.29 → -0.25 | -0.28 → -0.24 |
| LIN28A | NCCIT | -0.1358 | -0.07 → -0.10 | -0.07 → -0.10 | -0.06 → -0.11 | -0.06 → -0.10 |
| LIN28A | NTERA-2 cl.D1 | 0.2806 | -0.07 → -0.10 | -0.09 → -0.12 | -0.07 → -0.11 | -0.07 → -0.11 |
| LIN28A | PA-1 [PA1] | -0.0688 | -0.10 → -0.10 | -0.06 → -0.10 | -0.06 → -0.09 | -0.07 → -0.11 |
| LIN28A | SUSA | 0.2815 | -0.08 → -0.10 | -0.10 → -0.12 | -0.04 → -0.08 | -0.05 → -0.09 |
| DNMT3B | NCCIT | 0.0087 | +0.03 → +0.03 | +0.02 → +0.00 | +0.07 → +0.05 | +0.07 → +0.05 |
| DNMT3B | NTERA-2 cl.D1 | 0.2168 | +0.03 → +0.03 | +0.03 → -0.02 | +0.05 → +0.04 | +0.05 → +0.04 |
| DNMT3B | PA-1 [PA1] | 0.2115 | +0.04 → +0.04 | +0.02 → +0.05 | +0.06 → +0.04 | +0.06 → +0.05 |
| DNMT3B | SUSA | 0.0354 | +0.04 → +0.03 | +0.03 → +0.01 | +0.02 → +0.01 | +0.03 → +0.01 |
| PRDM14 | NCCIT | -0.6364 | -0.10 → -0.09 | -0.12 → -0.11 | -0.15 → -0.11 | -0.15 → -0.11 |
| PRDM14 | NTERA-2 cl.D1 | -0.5865 | -0.10 → -0.09 | -0.08 → -0.05 | -0.15 → -0.12 | -0.15 → -0.11 |
| PRDM14 | PA-1 [PA1] | -0.4364 | -0.07 → -0.07 | -0.07 → -0.03 | -0.12 → -0.05 | -0.12 → -0.05 |
| PRDM14 | SUSA | -0.1433 | -0.14 → -0.10 | -0.12 → -0.06 | -0.22 → -0.09 | -0.22 → -0.08 |
| ZIC2 | NCCIT | -0.1479 | -0.02 → -0.01 | -0.06 → -0.07 | -0.03 → -0.03 | -0.02 → -0.02 |
| ZIC2 | NTERA-2 cl.D1 | 0.0299 | -0.02 → -0.01 | -0.07 → -0.02 | -0.05 → -0.05 | -0.06 → -0.06 |
| ZIC2 | PA-1 [PA1] | 0.1923 | -0.02 → -0.02 | +0.04 → +0.03 | -0.06 → -0.05 | -0.07 → -0.06 |
| ZIC2 | SUSA | -0.1166 | -0.01 → -0.01 | -0.01 → -0.00 | +0.04 → +0.04 | +0.04 → +0.04 |
| L1TD1 | NCCIT | -0.4530 | +0.01 → +0.01 | +0.02 → +0.02 | +0.00 → -0.00 | +0.01 → -0.00 |
| L1TD1 | NTERA-2 cl.D1 | 0.0224 | +0.01 → +0.01 | +0.05 → +0.05 | -0.00 → -0.01 | +0.00 → -0.01 |
| L1TD1 | PA-1 [PA1] | 0.1168 | -0.02 → -0.01 | +0.00 → +0.02 | -0.05 → -0.02 | -0.05 → -0.03 |
| L1TD1 | SUSA | -0.0055 | -0.01 → +0.02 | -0.03 → -0.00 | -0.08 → -0.04 | -0.08 → -0.04 |
| DPPA4 | NCCIT | -0.0825 | +0.04 → +0.06 | +0.08 → +0.08 | -0.01 → +0.03 | -0.01 → +0.03 |
| DPPA4 | NTERA-2 cl.D1 | 0.0324 | +0.04 → +0.06 | +0.09 → +0.08 | +0.07 → +0.10 | +0.07 → +0.10 |
| DPPA4 | PA-1 [PA1] | 0.0905 | +0.08 → +0.08 | +0.06 → +0.11 | +0.04 → +0.06 | +0.04 → +0.06 |
| DPPA4 | SUSA | -0.2051 | +0.05 → +0.06 | +0.10 → +0.11 | +0.07 → +0.08 | +0.07 → +0.09 |

