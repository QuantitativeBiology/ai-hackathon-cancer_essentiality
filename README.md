(Human-written text starts)

This is an exercise in assessing the bioinformatics capabilities of agentic AI. It is part of a hackathon run by the Beltrao Lab at ETH Zurich, to which we (INESC-ID & IST – Universidade de Lisboa) contribute, in which Claude was used to run bioinformatics projects (https://github.com/evocellnet/ai-hackathon/). All results in this project were generated entirely by AI, with minimal human intervention. This repository is not intended to be read as a traditional publication or used as a validated bioinformatics tool.

The goal of this project was to predict CRISPR-Cas9 gene essentiality in human cancer cell lines from their molecular profiles (DepMap Public 26Q3), and to test whether the model generalises to non-transformed human cells. It was executed by Claude Opus 5.5, with delegation to Opus 5.5 and Sonnet 5.5 subagents.

None of the results have yet been verified by a human. We plan to verify the agent's methodology and benchmarking, and will update this section accordingly.

The rest of this repository is entirely AI-generated. It includes the code implementing the predictive model and an AI-written report in the style of a journal article, describing the final method, the path taken to reach it, and benchmarking against other methods.

(Human-written text ends)
---

# Predicting CRISPR gene essentiality from molecular profiles (DepMap 26Q3)

The project predicts per-gene CRISPR knockout effects (Chronos gene effect) for a cell line that has molecular profiles but has never been screened. It then tests whether the model transfers to non-transformed human cells: hTERT-RPE1, hPSC and others, using published screens.

## Status

**Complete (2026-10-28).** The final report is [REPORT.md](REPORT.md). The cancer test set (G4) and the non-transformed screens (G6) were each opened once; their numbers appear only in the two sections directly below. All other numbers are from development cross-validation (CV). Full record: [RESEARCH_TREE.md](RESEARCH_TREE.md).

**Post-hoc additions after the final report (2026-10-28).** Neither touched the test set or the external screens,
and no reported number changed.
- **Screen maps:** DepMap's screen maps confirm the library proxy used for library-stratified results
  ([docs/library_map_check.md](docs/library_map_check.md)).
- **Genomics input:** `predict.py` now accepts user-supplied genomics, and with it reproduces the evaluated headline
  ([docs/predict_genomics.md](docs/predict_genomics.md)).

| Milestone | State | Outcome |
|---|---|---|
| M0 Setup | done | Prompt, data checksums, external resource downloader |
| M1 Literature | done | [LITERATURE.md](LITERATURE.md) and four surveys in `docs/lit/` |
| M2 Design | done | [DESIGN.md](DESIGN.md) v3 frozen after 3 review rounds (`reviews/`) |
| G0/G1 Evaluation code and baselines | done | Frozen metrics, test-set guard, run registry; tag `freeze-eval` |
| G5 Transfer data | done | 16 non-transformed screens acquired and harmonised (QC failures excluded); transfer criteria pre-registered (tag `freeze-transfer`) |
| G2 First models and modality ablation | done | Ridge-type models beat both baselines; [docs/g2_results.md](docs/g2_results.md) |
| G3 Model ladder | done | Headline: plain average of ridge, its paralog per-gene stage and LightGBM (M5); [docs/g3_m5_results.md](docs/g3_m5_results.md) |
| Mutation-driven dependencies | done | Expression alone misses RAS/RAF hotspot addiction; the headline recovers most of it through LightGBM ([docs/mutation_dependencies.md](docs/mutation_dependencies.md)) |
| Per-gene predictability | done | Frozen tiers from the headline model: 1,632 high, 6,198 medium and 10,085 low genes, plus 519 not evaluable (1,054 / 921 / 143 / 41 of the variable genes). Every prediction ships with a confidence score and its tier ([docs/predictability.md](docs/predictability.md)) |
| G4 Freeze finalists, one-time test evaluation | done | Finalists frozen (tag `freeze-finalists`); one access each on 255 test lines; headline test P1 0.319 (see below); every number independently reproduced ([reviews/g4_review/REPORT.md](reviews/g4_review/REPORT.md)) |
| G6 Apply to non-transformed cells | done | Essential-gene ranking (Q3a) "transfers" by the pre-registered rule in RPE1 and hPSC, but the gain over the gene mean is small; screen-specific dependencies (Q3b) are inconclusive; independently reproduced ([reviews/g6_review/REPORT.md](reviews/g6_review/REPORT.md)) |
| G7 REPORT.md | done | Every decimal number is tied to a results file (`scripts/check_report_numbers.py --strict`); red-team review applied ([reviews/report_redteam_g7.md](reviews/report_redteam_g7.md)); LESSONS pruned after a fresh-agent review |
| M8 Post hoc (after G7) | done | Screen maps registered; the library proxy matches the real library combinations (adjusted Rand index 0.949), so no number changes. `predict.py` takes DepMap-format genomics and reproduces the headline's dev-CV predictions (max difference 2.4e-7) |

## Test set (G4, opened once, independently reproduced)

255 held-out cancer lines, grouped by patient and stratified by lineage. Every model was refit on all 1022 dev lines and accessed the test set once. Brackets give 95% patient-bootstrap CIs.

| Model | P1 variable | ΔP1 vs B1 | ΔP1 vs B2 | P2b nAP |
|---|---|---|---|---|
| B1 lineage mean | 0.166 [0.151, 0.181] | — | −0.083 | 0.075 |
| B2 expression kNN | 0.249 [0.231, 0.266] | +0.083 | — | 0.140 |
| M1 expression ridge | 0.310 [0.295, 0.324] | +0.144 | +0.061 | 0.180 |
| T transfer model (expression + TP53/RB1) | 0.311 [0.296, 0.325] | +0.145 | +0.062 | 0.180 |
| **M5 average (headline)** | **0.319 [0.304, 0.333]** | **+0.153 [0.138, 0.170]** | **+0.070 [0.062, 0.078]** | **0.184** |

Paired with M1, the headline gains +0.009 [0.009, 0.010] on the test set, the same as in CV. The gain is real but small, and it is confined to variable genes. Headline P1 by frozen per-gene tier: high 0.383, medium 0.196, low 0.106. Test refit time on 4 threads: 29 min for the headline, 1.4 min for M1.

## Transfer to non-transformed cells (G6, opened once, independently reproduced)

The transfer model T (expression plus TP53/RB1 status) was applied to published CRISPR screens of non-transformed cells. Expression and screens were harmonised to DepMap with a procedure frozen before any screen was compared with a prediction.

| Background | Question | Pre-registered verdict | Skill over the comparator [95% CI] |
|---|---|---|---|
| RPE1 (4 studies) | Q3a: ranking of the screen's top-10% depleted genes (ΔAP vs gene mean) | transfers | +0.032 [0.022, 0.040]; replicate ceiling 0.031 |
| hPSC (1 study) | Q3a | transfers | +0.009 [0.004, 0.014] |
| RPE1 | Q3b: recall of screen-specific dependencies (vs lineage mean) | inconclusive | +0.029 [−0.003, 0.061] |
| hPSC | Q3b | inconclusive | +0.008 [−0.029, 0.038] |

How to read these verdicts:
- **"Transfers" is weaker than the word suggests.** The gain over the DepMap gene mean is real and specific to the cell's own expression profile (it disappears with another cell type's expression), but it is small: about 5% relative AP. Every CI lies below the pre-registered minimum detectable effect. A "fails" verdict was unreachable, because no transformed anchor screen could certify the pipeline's sensitivity; this was known before the data were opened.
- **hPSC rests on a single screen.**
- **A simple baseline meets the same conditions.** Expression nearest neighbours (B2) also beat the gene mean, so "transfers" does not show that T beats generic expression similarity.
- **The RPE1 input was DepMap's own expression** of its RPE1 sublines, so the RPE1 result tests transfer across screens, not across expression batches.
- **Q3b:** context-specific dependencies are not detectable in the published screens. In DepMap's own screens of RPE1-derived lines they were (5 of 5 lines).
- **Non-cancerous DepMap lines** (5 RPE1-derived lines) are predicted as well as cancer test lines (median 93rd percentile).
- **Embryonal-carcinoma-specific dependencies** (NANOG, SOX2, LIN28A, PRDM14) in hPSC are not predicted. This was forecast before access.
- **Not evaluated:** no iPSC, iPSC-derived or primary haematopoietic screen. The iPSC and neuron screens are CRISPRi and had no anchor screen to calibrate them; the primary T-cell screens failed QC. The one hPSC screen is in H1 embryonic stem cells.

## Results so far (dev CV only, not the test set)

**Data and protocol:**
- 1022 development cancer lines, 5 folds grouped by patient and stratified by lineage.
- P1 is the mean per-gene Pearson r across held-out lines, on variable genes (SD > 0.25). It is computed within each fold and then averaged.
- P2b is the per-gene normalised average precision for essentiality (effect < −0.5).
- Fit time is summed over the 5 folds at 6 threads. The machine was shared with other jobs, so treat it as indicative.

| Model | P1 variable ± fold SD | ΔP1 vs ridge (M1) | P2b nAP | Fit time (5 folds) |
|---|---|---|---|---|
| B0 per-gene mean | 0 (by construction) | −0.308 | 0 | <1 s |
| B1 per-gene lineage mean (shrunk) | 0.184 ± 0.013 | −0.124 | 0.087 | ~2 min |
| B2 expression kNN | 0.253 ± 0.009 | −0.055 | 0.152 | ~2 min |
| R1 published-style ridge (raw target) | 0.307 ± 0.005 | −0.001 | 0.190 | ~6 min |
| M1 kernel ridge on expression, gene-mean residual | 0.308 ± 0.005 | — | 0.191 | ~5.5 min |
| M1b + copy number, mutations etc. (tuned weights) | 0.310 | +0.0015 | 0.193 | ~26 min |
| M2 + per-gene paralog stage | 0.3115 | +0.0034 | 0.195 | ~15 min |
| M3 LightGBM on M1 residuals | 0.310 | +0.0018 | 0.193 | ~35 min |
| M5 learned stacking of M1, M2, M3 | 0.318 ± 0.004 | +0.0099 | 0.196 | ~3.7 h |
| **M5 plain average of M1, M2, M3 (headline)** | **0.3175 ± 0.005** | **+0.0094** | **0.196** | **sum of its three inputs, ~1 h** |

**How the best model was chosen:** a more complex model replaces a simpler one only if it gains more than 0.005 in P1 and more than 1 SE. M2, M3 and the genomics weights each fall short. Their plain average with M1 clears it (+0.0094; better in all 5 folds), and learned stacking adds nothing over plain averaging. The simpler expression-only M1 is kept as a second finalist.

**Robustness:**
- Within library group, P1 is 0.289 for the headline M5 average and 0.279 for M1.
- Within lineage × library cells, it is 0.211 for M5 and 0.199 for M1.
- On held-out whole lineages, M1 scores 0.203 and beats B2 in 18 of 18 lineages. M5 was not run on this split.

**Other published methods,** run on the same genes as the best model:
- an elastic net is 0.030 below it;
- a TrPLet-like feature-selection method is 0.038 below it on its original feature grid, and about 0.02 below with a larger grid.

## Repository layout

| Path | Contents |
|---|---|
| `src/essml/` | Library code: data I/O, features, prior knowledge (paralogs, partners), baselines, frozen metrics and evaluation, test-set guard, screen and expression harmonisation, transfer evaluation |
| `src/essml/models/` | Kernel ridge (M1/M1b, R1, model T), per-gene stage (M2), LightGBM (M3), neural pre-check (M4), published methods (R2/R4), stacking (M5) |
| `scripts/` | Entry points: splits, features, `run.py` (common launcher plus registry), reports, harmonisation, gates |
| `configs/` | Locked split, seeds, evaluable genes and gene classes, run configs, frozen hashes, harmonisation and transfer criteria |
| `results/` | Run registry, test-access log, per-run per-gene metrics, gate evidence, harmonisation outputs |
| `docs/` | Data audit, literature surveys, per-gate result write-ups |
| `reviews/` | Design review rounds and independent result reviews |
| `tests/` | Verification suite: toy data, positive and negative controls, leakage and blinding checks |

## How to reproduce

1. **Environment** (Python 3.13, CPU only):
   ```
   pip install -r requirements.txt
   python -m pytest -q
   ```
2. **DepMap data:** download the DepMap Public 26Q3 files listed in [DATA.md](DATA.md) from the DepMap portal into `data/`.
   - The portal blocks scripted downloads, so this step is manual.
   - Verify them with `cd data && sha256sum -c ../configs/depmap_local_sha256.txt`.
   - Parquet caches under `cache/` are built automatically on first use.
3. **External resources:**
   - `python scripts/download_external.py` fetches paralogs, PPI, gene sets and the Hart lists.
   - `python scripts/fetch_external_screens.py` fetches the published non-transformed screens and expression data.
   - Sources, versions, checksums and licences are listed in [DATA.md](DATA.md).
4. **Features:** `python scripts/build_features.py`. The split (`configs/splits.csv`) and the gene lists are committed; `scripts/make_splits.py` regenerates the split.
5. **Runs:** `python scripts/run.py --config configs/runs/m1_E_resB0.json --threads 6`.
   - Every run is logged in `results/registry.csv` with config hash, commit and data hashes.
   - Reports: `scripts/g2_report.py`, `scripts/g3_kernel_report.py` and `scripts/g3_report.py`.
6. **Predict for new lines:** `python scripts/predict.py --expression my_lines.csv --out predictions.csv`.
   - The output gives each gene's predicted effect, confidence, tier and an applicability-domain distance.
   - For non-DepMap profiles, add `--harmonise quantile`.
   - It serves the headline m5_avg by default (`--model m1_E_resB0` for the expression-only model).
   - **With genomics:** add `--cn`, `--damaging`, `--hotspot` and `--metadata`. These take DepMap-format tables,
     and DepMap release files are accepted as they are. With all four, the script computes the evaluated headline,
     including its RAS/RAF hotspot recovery. Partial inputs and calls from other pipelines were not evaluated
     ([docs/predict_genomics.md](docs/predict_genomics.md)).
   - **Without them,** m5_avg's genomic inputs are diploid wild type (DESIGN §2); a warning says so.
   - This expression-only use was evaluated post hoc in dev CV only: P1 0.3135, against 0.3081 for the expression-only ridge and 0.3175 for the headline with genomics. It has no test-set number, and it does not recover RAS/RAF hotspot dependencies ([docs/predict_expression_only.md](docs/predict_expression_only.md)).
7. **Test set:** `src/essml/guard.py` refuses test and external-screen targets unless the run's config hash is listed in a committed `configs/finalists.json` and `ESSML_UNBLIND=1` is set. Each access is logged in `results/test_access_log.csv`.

## Data and licensing

- **No DepMap data is in this repository.** `data/` and `cache/` are git-ignored, and `tests/test_no_data_tracked.py` checks that no raw data, binary matrices or files over 5 MB are tracked.
- Only derived summaries are committed: per-gene metrics, gene lists, and line metadata for the split.
- DepMap Public data are CC BY 4.0.
- External resources keep their own licences, listed in [DATA.md](DATA.md). Some published supplementary tables may not be redistributed, so they are downloaded by script and never committed.
- No CRISPR/RNAi data from evaluation lines is used as model input. Text-mining resources such as STRING are not used.

## Key documents

| Document | Contents |
|---|---|
| [PROMPT.md](PROMPT.md) | The task as given, verbatim |
| [DESIGN.md](DESIGN.md) | Pre-registered design: splits, metrics, models, transfer criteria |
| [RESEARCH_TREE.md](RESEARCH_TREE.md) | Progress tree, status summaries and decision log |
| [LESSONS.md](LESSONS.md) | Process lessons for agent-run research |
| [LITERATURE.md](LITERATURE.md) | Literature synthesis |
| [DATA.md](DATA.md) | Data sources, versions, checksums and licences |
| [REPORT.md](REPORT.md) | Final report, journal-style |
