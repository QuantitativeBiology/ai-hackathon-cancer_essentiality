# Round 1 review: algorithms and compute

## 1. Verdict
The ladder is sound and CPU-feasible. B0, B1, B2, M1, M2, M3, M5 are the right families. M4 is the highest-risk rung. Its architecture is under-specified and its 2 h/CV-run estimate is optimistic. The main gap is that the ladder lacks the cheap strong alternatives. These are reduced-rank and low-rank-target variants, other kernels, and per-gene or per-block λ. My timing tests (dev fold 0, 814 train / ~208 val, 15.8k standardised expression features, residual from the gene mean only, not from B1) show these give almost nothing over the plain linear kernel. That is itself a finding worth pre-registering rather than assuming.

Timing experiments (dev only, 12 threads, float64):

| Experiment | Result |
|---|---|
| Kernel (814x814) + `eigh` | 1.4 s. One eigendecomposition then serves any λ grid and any target set. |
| Shared λ sweep | Mean r (all genes / top-10% variable): λ=3e3 0.133/0.296; 1e4 0.146/0.314; 3e4 0.154/0.318; 1e5 0.152/0.303. |
| Per-gene λ via closed-form LOO, 15-point grid | 0.152/0.314. No gain over the best shared λ. 6222 genes sit at the grid edge (λ≈3e5). Cost 7 s for 15 λ values on all 15.8k genes. |
| RBF (bandwidth g·median sq. distance, ridge 0.1–1) | Best 0.157/0.317, essentially identical to linear. |
| Low-rank target, truncated SVD at rank 20/50/100 | 0.156/0.276, 0.157/0.298, 0.157/0.308. Slight gain on all genes, loss on variable genes. |

## 2. Proposed changes

1. **§5 M1, priority must.** Do not use a "λ per gene by inner CV".
   - Do the `eigh` of K once per fold. Compute closed-form LOO residuals, e_i = (r_i − ŷ_i)/(1 − H_ii), for a log grid of λ, each at O(n²) per λ for all genes at once. LOO gives the per-gene choice at almost no extra cost.
   - Extend the grid to at least 1e2–1e7. The optimum was at the edge in my test.
   - Treat shrinking λ to a few groups of genes (variable vs the rest) as the default. Per-gene λ is noisy and gave no gain.
   - Note: LOO within the training fold is not patient-grouped. Use grouped folds (leave-patient-out, via block-wise LOO) if lines of the same patient exist in dev. At minimum, collapse near-duplicates before LOO.
2. **§5 M1, priority must.** Fix the target definition.
   - Fit the ridge on R = Y − B1(lineage), with B1 computed inside the training fold. The line's own lineage mean must be built in a leave-line-out way when computing the training residual. Otherwise training residuals are too small, because the line is in its own group mean. In small lineages this is a real optimism and a mismatch between training and test.
   - Choose the shrinkage κ by the same grouped CV. The probe (`scripts/probe/probe_dev_fold0.py`) did not use a leave-line-out lineage mean. Do not reuse its numbers as a ceiling.
3. **§5, priority should.** Add **M1b**: RBF/correlation kernel and a sum-of-kernels variant (one kernel per modality: E, C, M/H/F, meta), with kernel weights fixed by inner CV. This replaces concatenating modalities into a single standardised matrix. Rationale: modality blocks have very different dimensions (15.8k expression vs sparse mutations). Under plain concatenation, expression swamps the others and the ablation understates them. Use a block weight w_b on each block's standardised features, scaled by 1/√p_b. Expect ≈ equal performance to linear on expression alone (my test), so keep it to one extra rung, not a search.
4. **§5, priority should.** Add **M1c** reduced-rank regression as a *diagnostic*, not a rung. Truncated SVD on the ridge fitted values, or low-rank Y, was neutral to negative on variable genes in my test. If RRR is reported, report P1 on variable genes, where it loses, rather than all-genes mean r, where it wins slightly.
5. **§5 M2, priority must.** "≈18k small per-gene ridges, minutes" is fine, but design it as a stacked second stage with out-of-fold M1 predictions (nested, see item 9) rather than fitting on in-sample M1 residuals. In-sample M1 residuals shrink toward zero and make local-feature coefficients biased. Implementation: put the local features and the M1 out-of-fold prediction in one per-gene ridge. Use batched normal equations (p_local ≈ 10–30 features) with `np.linalg.solve` on stacked (G,p,p) arrays. This is <1 min, no Python loop needed.
6. **§5 M3, priority should.** Cheaper and safer LightGBM:
   - Set `n_jobs=1` per model and parallelise over genes with joblib (loky), or give LightGBM threads but run genes serially. Never both.
   - Preselect features inside each outer-train fold: ~200 features per gene (own-gene E/C/M, the paralog features, the top-correlated expression genes, 20 PCs). Use 100–300 trees, depth ≤4, learning rate 0.05, `min_data_in_leaf` ≥ 20, early stopping on an inner split.
   - Train on the M2 residual, not on Y.
   - Estimate: 2k genes × 5 folds × ~1 s ≈ 3 h single-thread, ≈ 15–20 min on 12 processes. Fine. Report the gain over M2 on the same genes. Prior experience suggests small (<0.01 mean r) gains. Set a go/no-go: drop M3 if the paired gene-level gain is not positive in ≥4/5 folds.
7. **§5 M4, priority must.** Specify the architecture and make it CPU-feasible. A (cell × gene) pair model with 1022×18k = 18M pairs is wasteful on CPU. Prefer a factorised model:
   - Cell encoder: MLP from a *reduced* input (top 5k variable genes of E, or 256 PCs of E plus C/M summaries) to d=64.
   - Gene side: a free per-gene embedding plus a fixed prior from training-fold Y (SVD loadings of the gene axis, rank 64) and from paralog/complex features.
   - Output: y_hat[c,g] = b_g + <u_c, v_g> + local-feature term. This is the matrix-factorisation form. One forward pass over a cell mini-batch (batch 64) gives all 18k outputs by one matmul: 64×64 · 64×18k.
   - Cost: an epoch is 800 cells × (5k×512 + 512×64 + 64×18k) ≈ 0.7 GFLOP per cell fwd, about 3 TFLOP/epoch fwd+bwd. This is ~1–3 min/epoch on 12 cores, 50 epochs ≈ 1–2 h per fold. Reduce with PCA input (256-d): ≈ 5 s/epoch. Use the PCA input as the default.
   - Loss: masked MSE (a NaN mask, no filling). Weight decay and dropout on the input. Early stop on an inner split.
   - M4's only possible edge over M1 is nonlinearity. Pre-register a cheap check first: does an MLP on the top-256 PCs of E predicting the top-100 PCs of R beat the kernel ridge? If not, drop M4.
8. **§5 B2, priority could.** kNN on correlation is a cheap rung (milliseconds). Also use it as an input feature for stacking (neighbour-mean residual). Keep.
9. **§5 M5, priority must.** Stacking without leakage.
   - Level-1 predictions for stacking must be out-of-fold over the *training* set (nested: outer test fold, inner K-fold grouped by patient).
   - Weights: non-negative per-gene-class weights (variable vs other), not per gene (noisy). Plain averaging of the top 2–3 is the baseline for the stacker to beat.
   - Cost: stacking multiplies fitting by (inner folds); M1/M2 are cheap, M3/M4 are the issue. Reuse the *same* inner-fold OOF predictions for λ selection, M2, and the stack.
10. **§2 / §10 NaN targets, priority must.** Filling NaN with the training-fold gene mean in Y is acceptable only for the unmasked closed-form ridge. It shrinks apparent variance and biases r. The sensitivity check (dropping genes with >10% NaN) should also be the *default for the fit*. For the stacker and M4, use masks. Evaluation keeps observed pairs only (already stated). Compute per-gene r using only lines with observed targets; skip genes with <20 observed held-out lines.
11. **§10 Pitfalls, priority must.**
    - **Standardisation inside folds:** the probe does it correctly (`xm, xs` from train). Keep zero-variance filters (`std>0.5`) in the fold as well. Cap standardised features at ±5–10 so rare mutations/fusions do not dominate. Use float64 for K and `eigh` (float32 `eigh` was not tested, but conditioning at λ=1e2 is poor).
    - **Numerical stability:** use `eigh` or Cholesky on K + λI. Never invert. Add jitter. Avoid `np.linalg.solve` repeated per λ at 0.1 s each when `eigh` once costs 1.4 s. Mean-centre X and Y with training means; the intercept is then unpenalised. Do not penalise the lineage-indicator block.
    - **Threading on Windows:** set `OMP_NUM_THREADS`, `MKL_NUM_THREADS`, `OPENBLAS_NUM_THREADS` *before* importing numpy. For joblib/LightGBM stages set them to 1 and use 12 workers. For the BLAS stages (kernel, `eigh`) use 12 threads in a single process. Do not mix: use `threadpoolctl.threadpool_limits` as a context manager per stage. In torch, call `torch.set_num_threads(12)` and `set_num_interop_threads(1)`. Use `if __name__ == "__main__"` guards (Windows spawn). Large arrays passed to loky workers are memmapped, so pass fold indices instead of copies.
    - **Metric ties:** constant predictions give NaN r; `nan_to_num` → 0 is consistent with the design.
12. **§10 Memory/runtime table, priority should.** Replace "seconds / minutes" by the measured/estimated numbers below, logged by the harness.
    - X (1022 × 15.8k) float64 = 130 MB; K = 8 MB; Y (1022×18.4k) float64 = 150 MB. Even with nested copies (×10) it is <4 GB. M4 peak with PCA input <3 GB. No memory risk anywhere except external-data preprocessing.
    - Runtime per fold (12 threads): B0/B1 <1 s; B2 <1 s; M1 ≈ 1.4 s `eigh` + ~7 s for a 15-λ LOO grid + ~1 s predict; M1b similar plus kernel build (cdist ≈ 3 s); M2 <1 min; M3 15–20 min with 12 processes; M4 5–15 min with PCA input; M5 inner-fold multiplier ×5 for M1/M2 (≈ 5 min), M3/M4 reuse.
    - Full ablation (~20 configs × 5 folds) on M1 costs ≈ 20 × 5 × 15 s ≈ 25 min. Ablation on M3/M4 should be skipped; do it on M1/M2 only, as the design states.
13. **§6 Ablation, priority should.** Ablation with a kernel sum needs the block-weight scheme of item 3, otherwise "leave-one-out" changes the total feature count and hence the effective λ. Re-tune λ per configuration (cheap with `eigh` reuse, ~10 s).
14. **§4 Selection, priority could.** The 1-SE rule with 5 folds is fine. Because LOO-λ costs nothing, report the learning curve (train size 25/50/100%) once, to show whether more data would help. It helps interpret the ceiling.

No change needed to B0, B1 (apart from item 2), B2, the evaluation metrics, or the split design from an algorithmic standpoint.
