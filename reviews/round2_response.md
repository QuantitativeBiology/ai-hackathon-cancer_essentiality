# Round 2 — PI response (DESIGN v1 → v2)

Verdicts: external PI minor; statistics accept after 3 fixes; prior work minor; data minor; algorithms no blocking
issues; agent failure modes targeted revisions; biology minor. All must items were adopted except where noted.

## Conflicts resolved
- **Screen-specific dependency definition (§7.5).**
  - Statistics asked for a dependency fraction ≤ 0.5; biology asked for ≤ 5%.
  - Decision: ≤ 0.5 is the headline and ≤ 5% the secondary.
  - Why: ≤ 0.5 removes the core-essential inflation both reviewers worried about while keeping enough positives
    per screen for a recall estimate; ≤ 5% is reported as the strict view.

## Modified or rejected
- **Biology "should": assign set-A tiers functionally (TP53 KO > +0.5, MDM2 < −1)** — R. This would condition the
  evaluation set on the very outcomes being predicted. Tiers stay mechanism-based (metadata); the functional p53
  check is reported descriptively.
- **External PI #4: refit frozen models on dev+test?** — decided as follows. The test set uses the model refit on all
  dev lines. Set A and set B use the average of the 5 dev fold-models, so they match the anchors' out-of-fold
  training size. The test set is never used for training.

## Adopted (by section)
- **§2:** CN transform log2((CN+0.1)/1.1) with 0 = diploid; instability on raw |CN−1| > 0.3 (data #2–3).
  Expanded p53/Rb composites (biology should).
- **§3:** A1 labelled a single-background case study (statistics #3).
- **§4:** fixed `configs/eval_genes.csv` with mask-only test counts (data #1); prevalences over observed lines.
- **§5 M1:** NaN-pattern merging, B1 for <50-row patterns, sub-kernels of one fold K, M1b weight tuning, R2 time box
  (algorithms #1–3).
- **§5 T:** no lineage as the primary variant (biology #6).
- **§5:** frozen-model protocol (external PI #4).
- **§7.1:** K562-only CRISPRi anchor (prior work should); neuron proliferation/maintenance split (biology should);
  committed screen plan (external PI #3).
- **§7.2:** location-only correction below 8 anchors; leave-anchor-out skill (external PI #2); p53/DDR genes excluded
  from correction (biology #3).
- **§7.3:** Celligner-style candidate; GTEXDEPMAP cited (prior work #1).
- **§7.4/7.5:** primary statistics named (external PI #1, statistics #1); B2 comparator (prior work #2); USP7 dropped
  from H1 (biology #2); CN null expectation (biology #4); H4 re-anchored on Shifrut T cells (biology #5).
- **§7.6:** "transfers" defined; MDE condition; block bootstrap; count rule; swap-null cap (statistics #2,
  external PI #1).
- **§9:** guard bound to finalist config hashes and a clean tree, plus set A; static CRISPR-read test; freeze tags;
  harmonisation freeze file; gate evidence files; registry start rows and extra columns; scope ledger;
  implementation check (agent failure #1–7).
- **Fixed in text:** LITERATURE.md Pacini anchor statement; untagged-note warning; stale probe numbers removed from
  the RESEARCH_TREE status (agent failure, "accepted but not done").
