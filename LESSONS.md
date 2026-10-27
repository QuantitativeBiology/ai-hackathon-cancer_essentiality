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
11. Shared append-only logs need a file lock as soon as runs go parallel. [registry.csv got a corrupted row when G2 switched to parallel workers.]
12. Always run leave-one-out ablations alongside add-one; a weighting default can make every addition look harmful. [G2: every add-one modality hurt, while leave-one-out showed C +0.009 and E +0.021, Holm p < 0.001.]
13. Run the published simple method as a diagnostic beside each pre-registered choice; it exposes design choices that cost accuracy. [The pre-registered B1-residual target cost 0.008 P1 against raw ridge (R1), 5/5 folds.]
14. Time only on a quiet machine; parallel agents make runtimes noisy, so re-time finalists alone. [G2 fit times varied about 1.5× with 6 vs 12 threads and shared jobs.]
15. Before calling a corrective stage no-go, scan its scale (base + c·correction) on the out-of-fold predictions; a mis-scaled real signal looks like a small gain. [M2/M3 at full scale gave +0.0034/+0.0018; their plain average, i.e. 1/3-shrunk corrections, gave +0.0094.]
16. A permuted-input control shares one permutation across all genes, so its chance level scales with 1/√(n lines), not 1/√(n cells); set bounds accordingly. [m5_ctl_perm weight 0.041 against a 0.01 bound; other seeds 0.035–0.037.]
