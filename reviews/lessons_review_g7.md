# Fresh-agent review for the G7 pruning of LESSONS.md

Reviewer: a fresh agent with no stake in the project, 2026-10-28 (PROMPT.md: prune LESSONS at milestones with a
"what to stop doing / what worked" review). The harness blocked subagent file writes, so the PI transcribed this from
its hand-back. The PI adopted the proposed LESSONS.md unchanged (16 rules).

The reviewer read:
- PROMPT.md, LESSONS.md, DESIGN.md §9–12;
- RESEARCH_TREE.md: status summaries, D0–D28, the tree;
- reviews/implementation_check.md, results/scope_ledger.csv;
- reviews/g4_review/REPORT.md, reviews/g6_review/REPORT.md;
- reviews/transfer_freeze_review.md and its response;
- the full git log: 198 commits from 10-26 12:37 to 10-28 04:55, about 40 hours.

**Process note:** LESSONS was pruned only once before (at M2), and the fresh-agent question had not been asked
before G7. That is itself a lapse against PROMPT.md.

## 1. What the agents should stop doing

1. **Promising process safeguards in prose that no code enforces.** DESIGN §9 promised a scope ledger, an exclusion
   log, an implementation check at G0 and G4, a gate file for every gate, stored review prompts and a fixed status
   block. None was produced during the project. The ledger, exclusions and implementation check were rebuilt after
   the fact at G7 (c69aaf5, D28). There is no G5.json. The `_status` fields of the frozen files never flipped, and the
   registry has no training-ID hash. The ledger shows 8 arms with no recorded fate.
   - LESSONS rule 7 ("check accepted items at the next gate") was on the books from M2 and was never acted on.
     Rule 19 restates the same failure.
   - The safeguards that were code all ran and paid back: guard.py, gate.py and check_report_numbers.py.
2. **Accepting review items in bulk without an owner or a tracker.** 174 items were accepted. Of these, 31 are
   partially implemented and 9 are not implemented (implementation_check.md summary).
   - One response says an item was done when it was not: R1-BI-17, Cas9 delivery mode "recorded in the inventory".
   - The three design rounds were committed within 33 minutes (1a46561 14:19 → 4196826 14:52). Accepting was cheap
     and delivering was not.
3. **Freezing decision rules without dry-running every branch, and turning down the reviewer who said so.**
   - Before the freeze, the transfer reviewer said the 50%-of-ceiling condition was near-vacuous and that
     "transfers" was easier to reach than "fails" (transfer_freeze_review "Points for the PI" #1). The floor was
     declined (F5).
   - After the freeze, the dry run showed "fails" could not be reached (D21e). After G6, "transfers" had to be
     qualified in wording (D27; g6_review M1–M2).
   - Other late discoveries of the same kind:
     - POU5F1, H3's key gene, is absent from the only hPSC screen (g6 m6).
     - Only 4 of the "7 EC dev lines" exist (D20d).
     - The expression-harmonisation criteria tied on every candidate (D9 override).
4. **Building things that depend on the final model before the model was chosen.**
   - Predictability tiers were frozen on m1_E_resB0 (2410d36), then frozen again on m5_avg (4897899/6c5e004).
   - The characterisation ran three times, with 23 anchor reads each time (5b0ad1b, fa5dd40, e935165).
   - The first G4 dry run was deleted (4c9118c).
   - finalists.example.json first named m1 as headline (396f3ef).
   - G4, G6 and predict.py were re-integrated for m5_avg (593ead7–aed0b7e, about 10 h after the G4/G6 scripts).
5. **Letting parallel tracks share names and files.**
   - registry.csv got a corrupted row (rule 11).
   - D21/D22 were taken by the other track, so the G3 decision was renumbered to D23 (bf2c1ac).
   - The other track took over scripts/g3_report.py (85c600c).
   - The m5_avg integration was delayed "to avoid file collisions" (D23f).
6. **Setting control bounds by intuition, then "accepting them as interpreted".**
   - The m5_ctl_perm weight bound was 0.01 and the run gave 0.041 (D23b).
   - The within-lineage and lineage × library permutations fell outside their expected bands (D20b; 0.071 against
     an expected 0.03).
   - Each miss needed a post-hoc explanation. That is a credibility cost the project took on itself.
7. **Re-implementing the whole metric stack at every gate after metrics were frozen.**
   - After G0, every reproduction matched: G2 bit-identical, D24 all reproduced, G4 within 2.1e-8 over 298 values,
     G6 within 4.9e-5 over 288 values (rounding only).
   - The reviewer findings that changed decisions or wording all came from red-team questions: the G0 guard gaps,
     the G2 refutation of the line-offset idea (D16), and G6 M1–M3 (D27).
   - The fix is to rebalance, not to drop independent review.
8. **Appending to LESSONS without pruning or the fresh-agent question.**
   - LESSONS was pruned once, at M2. The header still says "last: M2". It grew from 9 to 19 rules across M3–M6.
   - PROMPT asks for the stop-doing / worked-well question at every milestone. A repo search finds no sign it was
     asked before now.
   - Several rules describe a single incident (rules 15 and 16).
9. **(Minor) Writing decision-log entries as mini-reports.** D21, D22 and D25–D28 run 20–45 lines each, and the tree
   is 642 lines. That works against the prompt's "keep your own context lean"; the detail belongs in docs/, with a
   link from the log.

Not on the stop list: the 56 one-line access-log commits (D21a). They are automatic and cost nothing. Both reviewers
rebuilt the unblinding order from them and the reflog.

## 2. What worked well (keep)

1. **The test split was locked on day 1** (3834bc9, before any model). The G4 reviewer verified it from a single
   commit.
2. **Blinding was enforced in code**: a guard tied to committed finalists, frozen hashes and tags, and an access log
   committed before each access. G4 and G6 each ran exactly once. No access was refused and nothing was re-run. Both
   reviewers checked the order of events from git alone (g4_review Task 2, g6_review §2).
3. **The selection rules were pre-declared and kept even when the outcome was unwelcome.** δ / 1-SE drove D14, M3 and
   M4 no-go (D19), the D24 "not met" call, and m5_avg chosen by the rule rather than by preference (D23f). Chronos was
   not adopted, following the margin set before the run (D11).
4. **Bad outcomes were forecast before any access**: H3 likely fails (D20e), "fails" unreachable (D21e). The G6
   reviewer could call the outcomes "forecast pre-access", which protected credibility.
5. **The independent review of the transfer criteria before the freeze** returned "FIX FIRST" (D13). It fixed the MDE
   quantity, the pseudoreplication and nine open choices. It was the largest single credibility gain.
6. **Dev-only dry runs, with fold 0 standing in for the test set, reproduced CV within 1.9e-9** before unblinding
   (f97baac, 2c1395e). The one-time evaluations then ran cleanly.
7. **Code reviews before unblinding found real bugs**: un-nested M3 early stopping (a38d6a3), training-ID recording
   (959d96f), and a silent skip of finalists (D22f). By contrast, the post-unblinding reviews found nothing critical.
8. **Cheap alternate views exposed design choices that cost accuracy**: the published ridge beside the ladder (D14),
   leave-one-out beside add-one (D15), and the plain average as a pre-declared M5 baseline (D23).
9. **Time boxes used probe-projected coverage rules.** The full and 500-gene R2 runs (108 h and 25 h projected), full
   R4 and M3 coverage were all decided by rule, with reasons in the ledger. Runtimes came from a quiet-machine refit
   (D25c).
10. **Number tagging plus check_report_numbers.py**: 356 tagged REPORT numbers pass --strict.
11. **One parallel review panel of distinct expertise** (rule 3) found problems that did not overlap.

## 3. What happens to each old rule

| Old rule | Fate | New # | Reason |
|---|---|---|---|
| 1 Data access | Keep, extended | 1 | Still general. "If obtainable" inputs lapsed silently (D28). |
| 2 Probes via shared eval code | Keep | 3 | General, and cheap to follow. |
| 3 Parallel first review panel | Merge with 4 | 4 | Two halves of one review-structure rule. |
| 4 Single fresh reviewer after panel | Merge with 3 | 4 | Same. |
| 5 Grep for stale copies | Keep, new evidence | 13 | It recurred (D20d, `_status`, D27 m5). |
| 6 Don't reopen decisions | Keep | 14 | General. |
| 7 Check accepted items at next gate | Merge into new 5 | 5 | It was never acted on; now tied to an owner and a check. |
| 8 Lock test split early | Keep | 2 | It worked; core rule. |
| 9 Stop adding process machinery | Rewrite, merge with 17 and 19 | 6 | Code-enforced machinery paid back. The real failure was prose-only promises. |
| 10 Independent metric cross-check | Rewrite | 8 | Keep the once-at-freeze check. Later full re-implementations found nothing. |
| 11 File lock for shared logs | Rewrite, broadened | 12 | It covers decision-ID and file collisions too. |
| 12 Leave-one-out with add-one | Merge with 13, 15 | 9 | One "second cheap view before a negative" rule. |
| 13 Published simple method as diagnostic | Merge with 12, 15 | 9 | Same. |
| 14 Quiet-machine timing | Keep, shortened | 15 | It worked (D25c). |
| 15 Scale scan of corrections | Merge with 12, 13 | 9 | Single incident; keep it as one example. |
| 16 Shared-permutation 1/√n bound | Rewrite, generalised | 10 | Too specific. The general failure is guessed control bounds (D20b, D23b). |
| 17 Gate checks as code | Merge with 9, 19 | 6 | This is the part that worked; it is now the rule. |
| 18 Reviewer reports as text | Keep | 16 | Harness fact; it happened twice. |
| 19 Ledger artefacts at design freeze | Merge with 9, 17 | 6 | A symptom of prose-only promises. |
| New | Add | 5, 7, 11 | Review-item triage; dry-run a rule's branches before the freeze; build model-dependent work only after selection. |

That leaves 16 rules: 8 kept (some shortened or extended), 3 rewritten, 3 merged groups and 3 new.
