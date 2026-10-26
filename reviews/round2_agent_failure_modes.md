# Review round 2: AI-agent failure modes (DESIGN v1, commit 48bd50d)

## Verdict

**Accept with targeted revisions.**

§9 and §12 adopt every round-1 must item. Several of them are still not enforceable as written.

**New risk 1: accepted items that were never done.** Round-1 #10 is marked accepted, but:
- LITERATURE.md:66 still says "~8 anchors already help";
- `harmonisation_labels.md` and `data_resources.md` still carry 1 and 0 verification tags;
- RESEARCH_TREE.md:16 still presents the flawed probe numbers as calibration. Only item 1.3 is flagged.

**New risk 2: plan size.** About 30 model, baseline and contrast arms plus 6+ hypotheses invite silently skipped arms, gates passed without evidence, and a REPORT covering only the arms that worked.

## Remaining concrete changes

**1. [MUST] §9 test guard: bind the guard to content, not existence.**
- `ESSML_UNBLIND=1` is settable by any agent; "finalists.json exists" is trivially met. The loader releases test, **set A** and set B targets only if all of these hold:
  - the calling run's config hash is listed in `finalists.json`;
  - `finalists.json` was committed in an earlier commit than HEAD;
  - the working tree is clean.
- Set A is currently missing from the guard. Add it.
- Add a static test failing if any file except `src/essml/io.py` reads `CRISPRGeneEffect` or the crispr cache (closes the `pd.read_csv` bypass).
- G7 checks that `test_access_log.csv` holds exactly one access per finalist plus baselines.

**2. [MUST] §9 frozen code: freeze against a git tag, and freeze transfer code too.**
- A hash an agent can update in the same commit is not a freeze. At G0, tag `freeze-eval`; a test asserts `git diff freeze-eval -- src/essml/metrics.py src/essml/evaluate.py` is empty unless a `D#` entry names the change.
- At G5, tag `freeze-transfer` the same way for transfer metrics and harmonisation code (§7.2–7.5).

**3. [MUST] §7.2 and §12: harmonisation and criteria as files the code checks.**
- §7.2 names no file. Add `configs/harmonisation_frozen.json`; set-B prediction code refuses to run unless it and `transfer_criteria.json` were committed before HEAD. Set A also requires `transfer_criteria.json`.
- G2–G3 model agents get no set-B observed scores; model changes motivated by set-B observations are labelled post hoc.

**4. [MUST] §12 gates: evidence files, checked by the next gate.**
- Each gate writes `results/gates/G<n>.json`: commit, evidence run_ids, evaluated pass condition, non-producing reviewer.
- G4 unblinding refuses without G0–G3 files; G6 refuses without G5.
- Replace vague pass conditions ("Choose features") with the §4 selection-rule output, written into the gate file.

**5. [MUST] §9 registry: complete the schema.**
- Write a row at *start* (`started`), so crashed and killed runs stay visible.
- Add columns: dirty-tree flag (dirty runs cannot be cited), splits.csv hash, roles accessed, training-ID hash.
- Flag the deliberately non-nested R4 run `leaky=true`; code excludes leaky runs from selection.
- REPORT states the number of dev configurations evaluated.

**6. [MUST] Scope ledger (new §12 column, or `results/plan_ledger.csv`).**
- One row per DESIGN arm (B0–B2, R1–R4, M1–M5, T, ablation contrasts, H1–H4, set-A tiers): `done / dropped (reason) / not reached`.
- REPORT.md includes the script-generated ledger; dropped arms are reported, never omitted.
- Mark the arms PROMPT.md requires; optional arms start only after the core arms reach G4.

**7. [MUST] Round-response audit.**
- For each A/M item in `round1_response.md`, record the implementing file and section; a non-producing agent checks the list before v2 is frozen.
- Fix the three gaps in the verdict now.

**8. [SHOULD]**
- §9 stopping: use the selection rule's δ **and** 1 SE; ideas are logged in the tree before they run.
- §9 reproduction: list the headline numbers explicitly in G4 and G6; no `[!]` mark before reproduction passes.
