# Post-hoc check: DepMap screen maps vs the library proxy (2026-10-28, after G7)

**Status:** post hoc, after the report was final (G7). Nothing was re-run, and no headline, metric or selection
decision changes. Script: `scripts/library_map_check.py` → `results/library_map_check.json`.

**Inputs:**
- The user supplied `CRISPRScreenMap.csv`, `ScreenSequenceMap.csv` and the 26Q3 release README on 2026-10-28
  (DATA.md).
- Each model was mapped to the libraries of the screens DepMap combined into its gene effects (CRISPRScreenMap →
  ScreenSequenceMap.Library).
- That mapping was cross-tabulated with the library proxy used throughout: `configs/library_groups.csv`, the dev
  lines' gene-effect NaN-pattern groups. Groups 0–3 are the four most frequent patterns and group 4 is the rest.
- No gene-effect values were read. Test lines enter only as aggregate counts.

## Findings
1. **The proxy matches the libraries closely.**
   - Each of groups 0–3 is exactly one library combination: group 0 is Avana only (621 dev lines), group 1 is
     Avana + KY (154), group 2 is Humagne-CD (143) and group 3 is KY (47).
   - Group 4 (57 lines) is mixed:
     - 31 lines in other multi-library combinations;
     - 17 Avana-only and 9 Humagne-CD-only lines whose NaN pattern differs from their library's main pattern.
   - Adjusted Rand index 0.949; 96.3% of dev lines are in their group's majority library.
2. **Test lines (aggregate) agree.**
   - G4's proxy counts were 167 / 32 / 36 / 9 / 11 for groups 0–4.
   - The real counts are 168 Avana, 34 Avana + KY, 41 Humagne-CD, 9 KY and 3 in other combinations.
   - Per-line test groups were not saved, so no per-line check was made.
3. **Library goes with screen length.** Median screen length is 21 days for Avana, 24 for Humagne-CD, 14 for KY and
   17.5 for Avana + KY. The "library" strata therefore also separate screen durations.
4. **Screens combined into the dev lines' gene effects:**
   - none is drug-treated;
   - 9 dev lines have an engineered screen;
   - 892 dev lines are 2D standard screens, and the rest include 90 organoid screens (77 3D, 13 2D) and other
     formats.
   - No reported result conditions on these.
5. **Set A (DepMap's RPE1-derived sublines):**
   - all six are single Avana screens of 21 days, flagged as engineered;
   - ACH-002466's screen fails QC and is not in the gene effects, which is consistent with the five set-A lines
     reported.
6. **Release identity.** The README's first line is "DepMap Public 26Q3". It holds no checksums.

## What this changes in what was reported
- **Within-library P1** (dev and test), the **lineage × library permutation control** and the **"49 NaN
  patterns, four major line groups"** statement:
  - For about 94% of dev lines, the proxy groups are the real library combinations, so the reported numbers stand.
  - Group 4 mixes libraries, so centring within it removes less of the library offset for those lines.
  - The size of that effect was not measured (no re-run), but it touches about 6% of lines.
- **"Library" means library combination.** For example, group 1 lines carry both an Avana and a KY screen. It is
  also confounded with screen length, so "library group" in REPORT should read as screen technology plus
  duration.
- **Predictability, "lineage means agree across DepMap library groups"** (REPORT §4):
  - The quoted comparison, groups 0 vs 1, is Avana-only vs Avana + KY lines. These share the Avana library, so it
    is not a comparison of independent libraries.
  - The comparison that shares no library, group 0 (Avana) vs group 2 (Humagne-CD, a Cas12 library), was already
    computed (`results/predictability/characterisation.json`, `lineage_mean_crosslib_02`).
  - It shows the same ordering: median r 0.324 high tier, 0.102 medium, 0.011 low.
  - The claim holds; REPORT now names the libraries and quotes both comparisons.
- **Transfer analysis:** unaffected. "Library-matched anchors" refers to the external screens' libraries (TKOv3,
  Brunello, …), not to DepMap's.
- **Not possible with these files:** re-processing DepMap screens from raw reads, which would also need the read
  counts and guide maps.
