# Review round 1: AI-agent failure modes

Reviewer: agent-failure-modes reviewer. Object: DESIGN.md v0 (commit 1a46561). Read: PROMPT.md, DESIGN.md, RESEARCH_TREE.md, LESSONS.md, DATA.md, LITERATURE.md, docs/lit/*.md, scripts/make_splits.py, scripts/probe/probe_dev_fold0.py, configs/.

## 1. Verdict

**Sound in intent, but the safeguards are not yet enforced.** The scientific controls are good: a locked patient-grouped test split, leakage assertions, positive and negative controls, separate CEG/NEG halves for scaling and evaluation, an anchor-based ceiling, and a registry. Almost every safeguard, though, is a sentence rather than a mechanism. Nothing in code stops an agent from reading test rows, editing the metric after seeing results, tuning harmonisation until transfer looks good, or quoting a number that no registry row backs.

There is already concrete evidence of the risk. The dev-fold probe (`scripts/probe/probe_dev_fold0.py`):
- calls `Y.fillna(Y.mean())` on all dev lines before splitting, so the validation-fold targets are imputed and validation lines feed the fill values;
- filters genes on NaN rate over train plus validation;
- picks the best of 4 λ on the validation fold itself;
- uses "top-10% SD" where P1 uses SD > 0.25.

Its numbers (0.15/0.31) sit in RESEARCH_TREE.md with no committed output and no review. Minor for a probe, but it is the "runs but computes something subtly different" failure, and it broke DESIGN §2 ("NaNs are never imputed in evaluation") in the first script that computed a metric.

## 2. Reference spot-check (WebFetch/WebSearch, 2026-10-26)

| Claim | Source checked | Result |
|---|---|---|
| GSE156284 = Mair 2019 H1 hPSC CRISPRseq | NCBI E-utilities esummary | **Verified.** Title matches; PMID 30970261. Deposited 2020-08, after the paper, which is fine. |
| GSE128210 = Brown 2019 RPE1 TP53 WT/null | NCBI E-utilities esummary | **Verified.** Title matches; PMID 31464370. |
| Chang & Zhang 2023: ridge ρ 0.276 vs DeepDEP 0.137 | bioRxiv full text | **Verified.** The correlation type for ρ is not stated explicitly. |
| Pacini 2021: 908 lines, 17,486 genes, 168 overlap; "8 anchors already help", 30–40 give r > 0.995 | PMC7955067 | **Verified with a nuance.** 8 is the *smallest subset tested* that still beat no correction, not a recommended minimum. LITERATURE.md §4 and DESIGN §7.2 lean on "~8 anchors suffice". |
| Dempster 2019: nAUC 0.53→0.91; 21 vs 14 d; late medians −0.402/−0.269 | PMC6925302 | **Verified.** |
| LEAP, *Sci Rep* 2026, doi:10.1038/s41598-026-56381-0 | doi.org | **The DOI resolves** to a nature.com article page (paywalled). Content not read. |

Estimated hallucination rate: **0 of 6 fabricated**; 1 of 6 mildly over-generalised.

Process gap: LITERATURE.md says the source notes carry `[opened]/[search]/[unverified]` tags. In fact `harmonisation_labels.md` has 1 tag and `data_resources.md` has 0, so their claims carry no verification status. `nontransformed_screens.md` honestly flags 18 `[unverified]` items, including several data locations the external-screen assembly depends on (Ihry 2019, Ting 2018, Shifrut 2018 GEO).

## 3. Proposed changes

**1. [MUST] §8: frozen evaluation module with hash check.**
- Before any model beyond B0/B1 is run, implement `src/essml/metrics.py` (P1, P2a/b, secondary metrics, observed-entry masking, bootstrap) and `src/essml/geneclasses.py`.
- Unit-test both against the toy cases in §8. Add a regression test that NaN targets are never scored and that fill values use training rows only.
- Record their sha256 in `configs/frozen_hashes.json`. A pytest test fails if a hash changes.
- Changing the module after freeze needs a decision-log entry and a re-run of every affected registry row, *baselines included*.
- Write the per-fold gene-class lists to `results/geneclasses/` once, at freeze, and have all models read them from there.
- *Rationale:* stops "silently changing metrics", which PROMPT.md forbids. Cost: an hour.

**2. [MUST] §3: test-set access guard and log in code.**
- All data access goes through `load_split(role)`. The default loaders drop `test` and `noncancerous` rows.
- Test rows are returned only when `ESSML_UNLOCK=test` is set *and* `configs/finalists.json`, committed earlier, lists the calling config hash.
- Every unlock appends to `results/test_access_log.csv` (UTC time, commit, config hash, caller).
- The final report states the number of test accesses, which should be one per finalist plus baselines.
- Probes and ad-hoc scripts must use the same loader.
- *Rationale:* "touched once" is currently an honour rule. Agents in long sessions forget it, or "just check" after a context reset.

**3. [MUST] §7: non-transformed sets are held-out sets too.**
- Freeze the model (weights, features, preprocessing; hash in `configs/finalists.json`) before computing predictions on set A or set B.
- Any model change after seeing transfer results is labelled post hoc and reported separately.
- *Rationale:* the 15 DepMap non-cancerous lines and the external screens are the second headline test, and v0 gives them no access rule.

**4. [MUST] §7.2: harmonisation is fixed blind to the model.**
- Choose every harmonisation option (step 1 scoring, whether to apply step 3, whether to add the step-4 kinetics term, quantile settings in §7.3) using only control genes, anchors and cross-study agreement.
- Commit the choices to `configs/harmonisation_frozen.json` *before* model predictions on external screens exist. A reviewer checks the commit order.
- *Rationale:* PROMPT.md's "harmonise before concluding failure" invites iterating harmonisation until transfer looks good. This is the design's largest forking-paths risk.

**5. [MUST] §7.4 and §4: numeric pre-registration of success and failure.**
- For each pre-specified hypothesis and each Q3 question, state the metric, the expected sign, and the threshold that counts as "transfers", "fails" or "inconclusive". Example: core essentials transfer if the model's held-out CEGv2/NEGv1 AUROC is ≥ 0.9 and within 0.05 of the cross-study ceiling.
- Make the failure-attribution rule numeric too ("exceeds the ceiling" by how much, with which CI).
- Commit all of this before external data is scored. The statistics reviewer can set the values.
- *Rationale:* prevents post-hoc reinterpretation and overclaiming.

**6. [MUST] §9: registry by construction, and numbers traceable to it.**
- Every run goes through one harness, `scripts/run.py`. It appends a row *at start* (`status=started`) and updates it at the end (`ok/failed/killed`).
- Row fields: run_id, commit, dirty flag, config hash, splits.csv hash, seeds, roles accessed, a hash of the training row IDs, metrics path, runtime and peak RSS, and the producing agent.
- Runs from a dirty tree cannot be cited.
- Every number in RESEARCH_TREE.md and REPORT.md cites a run_id. `scripts/check_report_numbers.py` compares the cited numbers against the registry.
- REPORT.md states how many dev configurations were evaluated.
- *Rationale:* addresses losing track of runs, fabricated or unverified numbers, and selective reporting. Failed runs become visible automatically.

**7. [MUST] §8: independent reproduction by reimplementation, with blinded prompts.**
- An agent that did not produce the headline numbers recomputes them (test P1/P2 for finalists and baselines; transfer AUROC/AUPRC). It reimplements the metrics from the §4 text in a separate script, reading only the committed prediction files and targets. Re-running the same scripts would reproduce the same bug.
- Tolerance is 1e-6.
- Reviewer prompts (stored in `docs/review_prompts/`) give the question, the files and the spec. They never give the PI's expected value, the probe numbers or the desired conclusion.
- For every `[!]` finding, a separate red-team agent is asked for the strongest reason it is wrong (leakage, confounding, selection).
- *Rationale:* subagents tend to confirm the PI's hypothesis, and re-running shared code reproduces shared bugs.

**8. [MUST] RESEARCH_TREE and §8: quarantine the probe numbers.**
- Mark 0.15/0.31 as `[probe, flawed: validation-fold imputation and λ selected on the validation fold]`.
- Do not use them as expectations or success thresholds.
- Add the probe's defects as named regression tests in change 1.

**9. [SHOULD] New §12: stopping rules.**
- Each idea enters the tree with its expected gain and cost.
- A model rung or branch is abandoned (`[-]`) when its dev P1 gain is < 1 SE.
- Stop exploring when 3 consecutive ideas give < 1 SE, or at the time budget. Only then write `finalists.json` and unlock the test set.
- *Rationale:* PROMPT.md's "keep going" with no criterion leads to unbounded scope drift and implicit dev-set overfitting.

**10. [SHOULD] LITERATURE.md and §7.1: verification policy.**
- Any claim that justifies a design choice, or appears in REPORT.md, must be `[opened]` with a quote and a location.
- `[unverified]` claims are leads only and never appear in REPORT.md.
- Re-tag `harmonisation_labels.md` and `data_resources.md`.
- An external screen enters the inventory only after a script downloads it and records URL, accession, checksum and licence in DATA.md.
- Correct the Pacini "~8 anchors" wording.

**11. [SHOULD] §4: post-freeze bug protocol.**
- If a metric, label or split bug is found after freeze or after test access: fix it, log a decision entry, re-run *all* affected arms including the baselines, and report both versions.
- Never re-run only the model.

**12. [SHOULD] §2 and §6: exclusion accounting.**
- Every filter writes counts and IDs to `results/<run>/exclusions.csv`. Filters include the >10% NaN genes, missing modalities, fusion partners >20%, the 1214-line subset and harmonisation drops.
- Paired comparisons assert identical gene and line sets.
- *Rationale:* silent denominators are a common source of fake gains.

**13. [SHOULD] RESEARCH_TREE.md: context-loss header.** A fixed status block lists what is frozen (with hashes), the test access count, open defects and pending reviews. The PI trusts a subagent summary only when it lists run_ids and output paths.

**14. [COULD] §8: per-run leakage fingerprint.** Log the training ModelIDs per run and check them against the test and non-cancerous IDs.

**15. [COULD] §7.2: fix the CEGv2/NEGv1 halves now.** Commit them, seeded, to `configs/control_splits.json` before any external screen is processed.
