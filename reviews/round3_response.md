# Round 3 — PI response (DESIGN v2 → v3, frozen)

Fresh reviewer verdict: accept with revisions.

## Must items
1. Frozen dev-wide gene classes — **M**. Used for *evaluation reporting* everywhere, so they are identical for all
   models. Fold-wise classes remain only where classes drive training decisions. This answers the statistics review's
   concern (selection on the outcome affecting model training) and the consistency concern together.
2. `eval_genes_external` — **A**.
3. Guard exceptions (control genes; per finalist × screen) — **A**.
4. Anchors predicted by the 4 non-training fold models — **A**.
5. Unset thresholds — **A**: committed at G5 in `transfer_criteria.json` from anchors/controls only, before any set-B
   prediction. Hypotheses without a committed criterion are descriptive only.
6. Verdicts per biological background — **A**.
7. Scale-free rank label; absolute AUPRC vs replicate ceiling as primary Q3a evidence — **A**.
8. Cut M4 — **R**. PROMPT.md asks for the best possible method, and M4 is the only rung that shares strength across
   genes. It is already gated by a cheap pre-check (MLP vs ridge) and time-boxed to 1 day.
   Reduced-rank/low-rank diagnostics — **A** (cut).

## Should items
- Median guide LFC as default; Chronos only where needed — **A**.
- Set A before set B — **A**.
- Defer M1b tuning, R3, Celligner — **A** (deferred, not cut).
- Drop the proteomics arm — **R**: it is cheap and already being built.
- 2-SE selection rule — **R**: the 1-SE + δ rule was agreed with the statistics reviewer in round 2, and changing it
  now would not be motivated by new information.
- Per-gate compute budgets — **A** (time boxes in §12 apply).

Design frozen at v3. Further changes require a RESEARCH_TREE decision note.
