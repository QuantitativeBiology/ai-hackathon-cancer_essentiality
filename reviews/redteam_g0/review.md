# Red-team review: G0 evaluation infrastructure and G1 baselines

Reviewer: independent agent (did not write the code). Date: 2026-10-26. Scope: DESIGN v3 §4, §8 and §9 checked against
`src/essml/{metrics,evaluate,genesets,guard,baselines,modelstore,features,priors,io}.py`, `scripts/run.py`,
`scripts/make_gene_classes.py`, `scripts/make_eval_genes.py`, `configs/{eval_genes,gene_classes}.csv`, `tests/`,
`results/registry.csv` and `results/runs/*`. I changed nothing in src/ or scripts/. The test suite passes: 78 passed.

## Bottom line

- **No critical bug was found in the metric numbers.**
- My independent recomputation of the B1 run agrees with the registry to every digit stored, fold by fold (see the
  cross-check section below).
- Two parts are **correct**:
  - the metric definitions: axis, NaN masking, per-fold averaging, AP sign, π, normalised AP;
  - the CV protocol: inner tuning, lineage means and B2 scaling all use training lines only.
- The major findings are about **guard bypasses** and **provenance**, not wrong numbers.

## Cross-check

Script: `reviews/redteam_g0/crosscheck.py`; output: `reviews/redteam_g0/crosscheck_output.json`.
- **Own code only.** The script does not import essml.metrics, evaluate or baselines, nor run.py. Its Pearson,
  tie-grouped step-wise AP, normalised AP and B1 are written from scratch.
- **Targets: dev only.** They come from `guard.load_targets(dev_ids)`.
- **Inner splits.** sklearn's StratifiedGroupKFold is used only to reproduce the inner-split indices; the scoring is
  my own code.
- **κ re-tuned independently.** It was re-tuned by my own inner CV, with the variable class recomputed on each inner
  training part.

| quantity | registry (`20261026-151647_b1_0e2a7a27`) | own metrics on saved OOF pred | own B1 + own κ tuning |
|---|---|---|---|
| P1 variable (mean of 5 folds) | 0.18374 | 0.183737 | 0.183737 |
| P1 variable SD | 0.01309 | 0.013091 | 0.013091 |
| P2a (class all) | 0.90627 | 0.906274 | 0.906274 |
| P2a SD | 0.00367 | 0.003670 | 0.003670 |
| P2b nAP, selective, prevalence ≤ 0.5 (extra) | 0.08675 | 0.086748 | 0.086748 |

- **Per-fold values.** P1 variable is 0.181592, 0.197083, 0.183532, 0.193108, 0.163368. These equal `per_fold.csv`
  and `summary.json`.
- **Gene count.** There are 2,118 variable ∩ evaluable genes in every fold.
- **Chosen κ.** My inner CV picked κ = {10, 10, 20, 10, 10}, the same as the run.
- **Predictions.** The largest absolute difference between the saved OOF predictions and my B1 is 2.4e-7, which is
  float32 storage.

Sanity checks:
- `configs/gene_classes.csv` *variable*: 0 mismatches against my dev-wide recomputation (SD with ddof=1 over observed
  lines, > 0.25).
- `configs/eval_genes.csv`: 0 genes fail the ≥40-per-dev-fold rule. There are also 0 genes that pass that rule but are
  missing from the list, so the 18,434 − 17,915 genes excluded are all explained by the test rule (≥50 observed test
  lines). I did not check the test half of the rule: that would need a mask call on test lines.
- Rare-selective is ignored, as instructed. The runs used **8,195** rare-selective evaluable genes (the old
  definition). The committed table now gives **2,991**.

## Findings

### Major

**1. The guard can be bypassed through its own parameters and private helpers, and the static test does not catch
it.**
- **Where:**
  - `src/essml/guard.py:182` (`load_targets(..., splits_path=, crispr_path=)`);
  - `guard.py:157` (`_read_rows`);
  - `guard.py:48` (`CRISPR_PARQUET`);
  - `tests/test_no_direct_crispr_reads.py:14`.
- **Failure scenario:** model code calls one of these and gets test or set-A targets with no `ESSML_UNBLIND`, no
  finalist and no log row:
  - `guard.load_targets(test_ids, splits_path=tmp_csv_marking_all_dev)` returns full test targets;
  - `guard._read_rows(test_ids, guard.CRISPR_PARQUET)` does the same;
  - `pd.read_parquet(guard.CRISPR_PARQUET)` does the same.
- **Why the static test misses it:** its regex only matches `CRISPRGeneEffect`, `crispr\.parquet` or
  `load_matrix('crispr')`. `CRISPR_PARQUET` has an underscore, so it does not match. Scripts already call private
  guard functions (`harmonise_screens.py` calls `guard._log`).
- **Fix:**
  - Make roles come only from the committed `configs/splits.csv`: drop the `splits_path` and `crispr_path` parameters
    from the public API, and inject them for tests with monkeypatch.
  - Extend the static test to flag `CRISPR_PARQUET`, `_read_rows` and `splits_path=` outside guard.py and tests/.
  - Better still, have `_read_rows` check its caller (or re-check roles) itself.

**2. `load_control_targets` releases control-gene values for the cancer TEST set without unblinding. This goes beyond
the §9 exception.**
- **Where:** `guard.py:241-256`.
- **What the spec allows:** DESIGN §9 allows control-gene reads "of external screens for harmonisation and QC". This
  function accepts any role, including `test`.
- **What it releases:** the control set is CEGv2 ∪ NEGv1 ∪ the 848 dev common essentials, about 2,000 genes.
  - Common essentials are most of the positives in per-line P2a.
  - CEG/NEG values give each test line its location and scale.
- **Failure scenario:** a calibration or harmonisation step calls `load_control_targets(test_ids)`, perhaps for "QC",
  and its output feeds a model or a threshold. That model's test P2a is then optimistically biased. Only a
  `control_only` log row records it.
- **Fix:** refuse role `test` in `load_control_targets` unless unblinded. Allow only `noncancerous` (set A) before
  the transfer freeze, and only if the design really intends that for set A.

**3. Committed results files hold CRISPR-derived per-line statistics for every line, test lines included.**
- **Where:**
  - `results/data_audit/crispr_line_stats.csv`: 1,361 lines, with per-line median and `frac_below_thr`. That column
    is exactly each line's P2a prevalence π.
  - `results/data_audit/crispr_gene_stats.csv`: gene mean, SD and fraction below −0.5 over all lines, test lines
    included.
- **Origin:** both were written by `scripts/audit_data.py` with `unguarded=True` before the split existed.
- **Failure scenario:** later code reads these "audit tables" and gets test-line target information with no guard
  call and no log row. Examples: a per-line essentiality offset, or a gene filter based on all-line SD. The static
  test cannot see this, because the files are CSVs under results/.
- **Fix:**
  - Remove the test, set-A and unassigned rows from `crispr_line_stats.csv` and rewrite `crispr_gene_stats.csv` from
    dev lines only, or move both to cache/.
  - Add `crispr_line_stats|crispr_gene_stats` to the static-test pattern.

**4. Provenance: the registry rows cannot be reproduced from the commits they cite, and the freeze was re-done after
the runs.**
- **Where:** `scripts/run.py:241-251` and `:341`; `results/registry.csv`; `configs/frozen_hashes.json`;
  `results/gates/G0.json`.
- **Dirty trees and wrong commits:** all four registry rows have `dirty_tree=1`.
  - B1 and B2 cite commit 721c1f5. Between that commit and the code actually used, `run.py` changed by 92 lines,
    `baselines.py` by 326 and `guard.py` by 234.
- **No tag:** there is no `freeze-eval` tag, which §9 requires.
- **Freeze rewritten:** `frozen_hashes.json` was rewritten at 16:14, after the B runs (15:09–16:11).
  - genesets.py's hash went from `7910…` (in G0.json) to `0d34…`.
  - metrics.py, evaluate.py and guard.py kept their G0 hashes, which is why the cross-check matches.
  - Baselines were not re-run, although the post-freeze protocol requires it.
- **Class tables not recorded:** the registry has no hash of `eval_genes.csv` or `gene_classes.csv`. That is why the
  stale rare-selective numbers cannot be detected from the registry.
- **Failure scenario:** REPORT.md cites a run id, and `check_report_numbers.py` passes. But the run's numbers came
  from uncommitted code and an older class table that nobody can reconstruct.
- **Fix:**
  - Tag `freeze-eval`.
  - Re-run B0, B1 and B2 on a clean tree at the tagged commit.
  - Add to the registry columns the sha256 of `eval_genes.csv`, `gene_classes.csv` and the four frozen modules (and
    the D8 table version).
  - Make `run.py` refuse a non-`leaky`, non-debug run when the tree is dirty, apart from the access log (see minor
    finding 1).

**5. Full set-B profiles are written to cache/ and can be read later without any guard. `controls_only` also loads
the whole table first.**
- **Where:**
  - `scripts/harmonise_screens.py`: `blind_write` writes `cache/harmonised/<id>.parquet`; `replicate_ceilings` reads
    it back with `pd.read_parquet`.
  - `guard.py:271` (`df = reader(path)` before the columns are subset).
- **Status:** this is a design-sanctioned write. But once written, the files can be read by any later script with no
  unblinding, and nothing checks for that. The static test only knows the DepMap file names.
- **Failure scenario:** a later "transfer QC" or model-selection script loads `cache/harmonised/*.parquet`. Set-B
  outcomes then shape the model or the criteria before `transfer_criteria.json` is committed.
- **Fix:**
  - Route reads of `cache/harmonised/` through a guard function, gated by the same conditions as
    `load_external_screen`.
  - Add `cache/harmonised` to the static-test pattern, with harmonise_screens.py on the allow-list.
  - In `controls_only`, give the reader the gene list, so that non-control genes are never materialised in a
    caller's closure.

### Minor

**1. The dirty-tree flag carries almost no information.**
- `results/test_access_log.csv` is tracked, so every logged `mask_only` or `control_only` access dirties the tree.
  This is why all registry rows have `dirty_tree=1`.
- `unblinding_check` will also refuse ("dirty tree") until the log is committed. A run that makes a mask call and
  then requests test targets will always fail.
- `--untracked-files=no` (`guard.py:117`, `run.py:247`) means an untracked script still counts as a clean tree.
- **Fix:** leave the access log out of the dirty check and count untracked `.py` files under src/ and scripts/.

**2. The finalist config hash is self-declared.**
- `config_hash` is a value the caller passes in, and it covers only the config JSON. It does not cover the code or
  `seeds.json` (`run.py:238`).
- Any script can pass a finalist hash it copied from the file.
- **Fix:** have `guard.load_targets` compute the hash itself from a config path, together with the code-tree hash.

**3. `splits_sha256` is fragile to line endings.**
- It hashes raw bytes (`run.py:254`).
- `configs/splits.csv`, `eval_genes.csv` and `gene_classes.csv` are stored with CRLF line endings in git (`i/crlf`):
  pandas wrote them on Windows.
- Regenerating them on Linux, or a future `.gitattributes`/autocrlf renormalisation, changes the hash without any
  change in content.
- **Fix:** use `sha256_lf` (from freeze_eval.py), or hash the parsed table.

**4. `FoldTransformer.fit` is not covered by the run sentinel.**
- `check_fit_ids` (`features.py:391-402`) blocks only test and noncancerous IDs.
  - It knows nothing about the held-out CV fold.
  - It passes silently if splits.csv is missing.
- Fits go through neither `run.guarded_fit` nor FIT_LOG. A future M1/M2 model that calls
  `FoldTransformer.fit(raw, dev_ids)` would standardise using the held-out fold, and no test would notice.
- **Fix:** take a `forbidden` ID set (the current Sentinel) and record each fit in FIT_LOG. Raise if splits.csv is
  missing.

**5. The F (fusion) absent flag is never set.**
- `fusion_matrices` builds its index as E lines ∪ lines with any fusion row (`features.py:297`).
- A line that was never fusion-profiled therefore looks the same as "no fusions", with absent = 0.
- **Fix:** take the index from the lines that were actually profiled.

**6. Partners outside the expression matrix are dropped from the mean.**
- `priors._partner_mean` (`priors.py:317`) leaves such partners out of the denominator.
- The docstring (`priors.py:24`) says a missing partner counts as expression 0.
- **Fix:** make the code and the docstring agree; choose one before M3 uses this feature.

**7. Training statistics are fitted on a different set of rows than they are applied to.**
- In `FoldTransformer._fit_num`, mean, SD and `min_count` use only the training rows present in the block.
- At transform time, absent training rows are zero-filled.
- So the standardised training columns are not exactly mean 0 / SD 1, and `min_count` ignores the absent lines. There
  is no leakage; this is a documentation and consistency point.

**8. Some pre-registered secondary outputs are not yet produced.**
- §4 lists within-lineage per-gene r (lineages with ≥20 held-out lines), results by NaN/library group, and the
  published-comparable block. `evaluate.py` produces none of them.
- No existing number is wrong. These must be added before the report.

**9. D8 timing.**
- The rare-selective 1% floor was adopted after B0, B1 and B2 rare-selective P1 values were already in the registry.
  The B2 run finished at 16:11; the freeze is dated 16:14.
- RESEARCH_TREE says it was decided "before any model comparison". The change is model-neutral, since the class is
  reporting-only, but the wording should say what had already been seen.

**10. Dev-wide classes include the held-out lines in CV.**
- This is sanctioned by §3. The variable and selective lists are chosen partly from the held-out fold's own values,
  which slightly favours models that predict extreme held-out values.
- It is the same for every model; note it in the report.
- Optionally, report fold-wise classes as a sensitivity check.

**11. The registry has no write lock.**
- `write_registry` does read-modify-replace with no lock. Two concurrent runs can drop each other's rows.

## Verified correct (no action needed)

**Metrics**
- Pearson is taken over lines for each gene: column-wise on [lines × genes].
- NaN truth entries are masked from both the mean and the deviations.
- Constant prediction or truth gives r = 0 by an exact test. Fewer than `min_obs` observations gives NaN.
- NaN predictions on observed entries raise an error.
- Spearman uses masked average ranks.

**Evaluation protocol**
- Every metric is computed within each fold and then averaged over folds; there is no pooling. The SD uses ddof=1.
- The fixed evaluable list is intersected with the class and with "≥ min_obs held-out lines".

**AP and normalised AP**
- Score = −pred and label = truth < −0.5, so more negative means more essential.
- π is the held-out observed prevalence. nAP = (AP − π)/(1 − π).
- Genes need ≥3 positives and ≥3 negatives. The headline subset is prevalence ≤ 0.5 on held-out lines.
- B0 gives nAP = 0 and AUROC = 0.5, as expected.

**Gene classes**
- Observed lines only; SD with ddof=1; thresholds as in §4, with D8 for rare-selective.

**Inner CV and baselines**
- Inner CV only sees `Y_tr`, with variable genes taken from the inner training part.
- The sentinel and FIT_LOG cover every fit in run.py.
- Tests confirm that perturbing held-out targets changes neither predictions nor the chosen κ.
- B1 has no leave-one-out issue: held-out lines are never in the lineage sums.
- B2 standardisation and top-variance selection use training lines.

**Guard and storage**
- `modelstore` reloads only dev targets.
- The CRLF-normalised hashing of the frozen modules is robust; `test_frozen` passes.
