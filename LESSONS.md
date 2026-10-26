# LESSONS — process rules for running this project

Each rule is one line, with its evidence in brackets. Pruned at each milestone.

1. Check data access (APIs, bot walls, licences) before planning around a source. [DepMap portal API returns a Cloudflare challenge; found in the first 5 minutes, not after the design.]
2. Even "calibration" probes go through the shared evaluation code, or are labelled unreportable at once. [The fold-0 probe tuned λ on validation and imputed with validation lines; caught by the R1 agent-failure reviewer.]
3. Run the reviewers in parallel with distinct expertise; they find non-overlapping problems. [R1: the data reviewer found library-structured NaNs, the statistician found the pooled-r artefact for B0, the PI reviewer found the circular RPE1 anchor; none overlapped.]
