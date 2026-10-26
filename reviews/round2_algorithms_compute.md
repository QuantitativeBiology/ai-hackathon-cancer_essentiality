# Round 2 review: algorithms and compute

## Verdict
Round 1 items 1-5, 7 and 8 are correctly absorbed. The pattern-wise NaN fit (item 6, modified) is statistically better than dropping genes and is affordable if implemented as below. There are no blocking problems, but four details need fixing.

Cost estimates, from my round-1 timings (eigh of an 814x814 kernel 1.4 s; 15-lambda LOO grid for all genes about 7 s):
- **Pattern-wise M1.** One eigh per pattern is 49 x ~1.4 s, about 70 s per fit. The LOO grid is paid once in total, because each gene belongs to one pattern. Pattern sub-kernels must be **submatrices of the single fold-level K**. Do not re-standardise features per pattern: that would change the kernel and break the equivalence unit test. Nested fits for M2/M5 (5 inner folds) multiply this to about 6 min per outer fold, so about 30 min for the whole CV. Acceptable.
- **Block-LOO by patient.** Patient groups need (K+lambda I)^-1 per lambda, about 0.05 s each. 49 patterns x 15 lambdas is about 40 s. Fine.
- **Memory.** Under 1 GB for the eigendecompositions. No concern.

## Remaining changes

1. **Pattern merging (must, section 5 M1).** Many of the 49 patterns probably differ in few lines or few genes. Merge any pattern whose row set differs from a larger pattern by <=1% of training lines. For merged genes, drop the extra rows by exact downdating or simply refit on the union row set. Record the number of merged patterns and genes. This cuts eigh calls by perhaps half and removes tiny-pattern noise. Patterns with fewer than 50 observed training lines get B1 only. State this rule explicitly, because the held-out >=50 rule in section 2 does not cover training.
2. **Unpenalised lineage (must).** An unpenalised lineage block does not fit the plain `eigh(K)` trick. The leave-line-out B1 residual already removes the lineage mean, so drop the "lineage unpenalised" wording. Otherwise an implementer will add lineage dummies to X and penalise them.
3. **M1b weight tuning (must).** Each weight vector needs a new eigh. Tune with a coarse grid on the training-fold LOO criterion (not an outer loop), on the **dominant fully observed pattern** and variable genes only, by coordinate descent over about 4 blocks x 3 levels, at most ~25 eigh, about 1 min. Fix the weights per outer fold, then reuse them for all patterns. Do not tune them nested inside the ablation. Use block weights 1/sqrt(p_b) as the default, to be compared against the tuned weights. The tuned run must also be reported at the default weights, so that the tuning gain is visible. Ablation then re-tunes lambda only (not weights), as written.
4. **Rung R cost (should, section 5).**
   - R1 is cheap (seconds).
   - R2 EN: per-gene `enet_path` with `precompute=True`, a fixed alpha-path, and warm starts. Variable genes number about 1.8k, times 5 folds, divided over 12 processes, roughly 1-2 h. Set `max_iter` and a tolerance, and run it with 1 BLAS thread per worker. Report non-convergence counts.
   - R3 RF: 500 genes x 5 folds x about 2 s on 1 thread, about 10 min on 12 workers.
   - R4 nested top-k is minutes.
   - Time-box R2 to 2 h. If it exceeds that, subsample to 500 variable genes, as for R3, and say so.
5. **Stacking inner folds (should, M5).** Compute inner-fold OOF predictions once per outer fold and share them between M2, M3, M4 and M5 (as already planned in round 1). The pattern-wise M1 inner fits are the dominant cost. Cache by (fold, pattern) key.
6. **Unit test (should, section 8).** Add a toy test where lines with NaN have wildly different values: the pattern fit must equal the per-gene ridge on the observed rows, including the LOO residuals, with intercept and centring recomputed on the observed rows only. Centring must use observed-row means for each pattern; a common failure is centring with the full-fold mean.

No further change needed to M2, M3, M4, M5, the threading plan, or the numerics list.
