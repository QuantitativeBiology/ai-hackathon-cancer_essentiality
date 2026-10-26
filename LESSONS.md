# LESSONS — process rules for running this project

Each rule is one line, with its evidence in brackets. Pruned at each milestone (last: M2, 2026-10-26).

1. Check data access (APIs, bot walls, licences) before planning around a source. [DepMap portal API returns a Cloudflare challenge; found in the first 5 minutes, not after the design.]
2. Even "calibration" probes go through the shared evaluation code, or are not run. [The fold-0 probe tuned λ on validation and imputed with validation lines; caught by the R1 agent-failure reviewer.]
3. Run the first review round with parallel reviewers of distinct expertise; they find non-overlapping problems. [R1: library-structured NaNs (data), the pooled-r artefact for B0 (statistics), the circular RPE1 anchor (PI); no overlap.]
4. After one full panel, use a single fresh reviewer, plus targeted reviewers only for new risks. [R2's 7 reviewers returned mostly "minor"; R3's one fresh reviewer found the unset thresholds and the inconsistent gene-class definitions.]
5. When a fact changes, grep for and fix every copy in the same commit. [R2 had to patch stale probe numbers and a wrong Pacini claim.]
6. Don't reopen a logged decision without new evidence; check the decision log first. [The M4 and 1-SE disputes recurred after agreement.]
7. Check accepted review items against code or text at the next gate, not at the end. [Round 2 found three "accepted but not done" items.]
8. Lock the test split and commit it before any modelling. [A 104-line patient-sharing structure made the split non-trivial; locking it early kept every later choice blind to test lines.]
9. Stop adding process machinery once the design is frozen; the next unit of value is a result. [DESIGN.md grew from 199 lines to ~440 before any model ran.]
10. Cross-check headline metrics with an independent agent writing its own metric code before freezing evaluation. [The red-team G0 cross-check matched B1 to 6 decimals and separately found 3 guard gaps that tests had missed.]
