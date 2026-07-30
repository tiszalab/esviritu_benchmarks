# Analysis Audit Handoff

## Instructions for the next agent

Work from this document on the HPC checkout. First verify every finding against the current code, committed input tables, raw tool outputs, and installed tool versions. Then propose a staged implementation plan for user approval.

Do **not** overwrite committed statistics or figures until:

1. the relevant behavior is covered by a focused test;
2. any tool-output unit assumptions have been checked against raw outputs and tool documentation;
3. the corrected analysis completes successfully; and
4. old and new results have been compared and explained.

This is an audit, not a mandate to accept every proposed remedy unchanged. Distinguish confirmed implementation defects from methodological choices that require the user's scientific judgment.

## Repository state at handoff

- Audit base commit: `f35b39d` (`main`, message: `update analyses`)
- Handoff branch: `analysis-audit-handoff`
- Repository: `tiszalab/esviritu_benchmark_rplots`
- Local review date: 2026-07-30
- No analysis code was changed during the audit.
- The local review could inspect code and committed outputs but could not rerun the analysis because the complete execution environment and raw benchmark outputs are on the HPC.

## HPC Locations
repo: `/data/tisza/analyses/mjt_projects/esviritu_benchmarks/esviritu_benchmark_rplots`
main tool outputs: 
  - `/data/tisza/analyses/mjt_projects/esviritu_benchmarks/benchmarks_mock/results`
  - `/data/tisza/analyses/mjt_projects/esviritu_benchmarks/benchmarks_negative1/results`
  - `/data/tisza/analyses/mjt_projects/esviritu_benchmarks/benchmarks_positive1/results`

Note: ignore scripts outside of the repo.

## Recommended workflow

1. Inventory the HPC data and raw output files used by each notebook/script.
2. Record exact versions of EsViritu, Kraken2/k2, Centrifuger, Metabuli, Sylph, R, Python, databases, and taxonomy files.
3. Add small synthetic fixtures and tests before modifying parser or metric behavior.
4. Resolve the authoritative mock-community scoring design with the user.
5. Fix confirmed parsing and unit defects.
6. Refit statistical models and regenerate outputs into a separate comparison directory.
7. Compare old versus corrected tables and figures.
8. Update documentation, dependency locks, and provenance manifests.

---

# Priority 0: decisions required before implementation

## A. Choose one authoritative mock-community scoring pipeline

There are currently two materially different implementations.

### Python implementation

- File: `scripts/mock_benchmark/metrics.py`, especially lines 20-28 and 140-278.
- Evaluated ranks: family, genus, species.
- ANI tiers:
  - 100/98 -> species
  - 95/90 -> genus
  - 85/80 -> family
- Implements RPM and within-viral relative Bray-Curtis.
- Excludes Sylph from abundance/Bray-Curtis because its abundance is coverage-based.
- Produces `rank_metrics.tsv`, `ani_aware_detection.tsv`, and `per_genome_detection.tsv`.

### R implementation

- File: `notebooks/benchmark_stats.Rmd`, especially lines 565-657.
- Evaluated ranks: family, genus, species, subspecies.
- ANI caps:
  - 100/98 -> subspecies
  - 95/90 -> species
  - 85/80 -> genus
- Independently computes `stats/profile_vs_truth_metrics.csv` from `unified_long.tsv` and ground-truth files.
- The committed mock models and figures use this R-derived table.
- `paths$ani` declares `ani_aware_detection.tsv` at lines 23-30, but that table is not used later in the notebook.

### Required decision

Define and document:

- the biologically justified resolution threshold at each ANI;
- whether subspecies is a valid comparable rank for every tool;
- whether the primary abundance endpoint is fixed-denominator RPM, within-viral composition, or both;
- how coverage-based Sylph abundance should be handled; and
- which generated table is the single source of truth for downstream models and figures.

Do not maintain two silently divergent scoring implementations.

## B. Define the estimand for parameter selection

The headline settings are manually specified in `notebooks/benchmark_stats.Rmd`, lines 42-54. Determine whether these were:

- vendor defaults;
- prespecified recommended settings;
- selected using these benchmark outcomes; or
- selected to optimize a sensitivity/specificity tradeoff.

If selected after inspecting these data, label the headline comparison as tuned/exploratory and add a prespecified/default-setting comparison. Do not select parameters and estimate confirmatory performance on the same observations without disclosure.

---

# Priority 1: confirmed or high-confidence correctness issues

## 1. Redesign the R divergence-capped truth calculation

### Location

`notebooks/benchmark_stats.Rmd`, lines 624-650.

### Current behavior

Truth relative abundances are calculated before rank-specific eligibility filtering. In capped mode, truth genomes that are not eligible at the current rank are removed with:

```r
filter(d, rank_order[cap_rank] >= rank_order[[rk]])
```

Predictions are neither folded to the permitted ancestor nor filtered by a matched eligibility rule.

### Consequences

- The retained capped truth vector can sum to less than one.
- A correct finer-rank prediction for an excluded divergent genome can be counted as a false positive.
- The calculation deletes truth taxa rather than implementing a common taxonomic cap.
- The resulting strict-versus-capped figures do not cleanly estimate performance at realistic recoverable resolution.

### Expected correction

Design a transformation that applies a compatible rank representation to both truth and predictions. Likely approaches include:

1. roll each truth genome and its corresponding prediction to the genome's permitted ancestor; or
2. define an explicit rank-eligible evaluation universe, then filter and renormalize both sides consistently.

Add tests for:

- a low-ANI genome predicted at species level;
- a low-ANI genome predicted only at genus/family level;
- mixed ANI tiers in one community;
- multiple genomes sharing the same capped clade;
- truth and prediction vectors summing as intended after transformation.

## 2. Fix Sylph positive-control detection

### Location

`scripts/positive/analyze_positive_controls.py`, lines 244-276.

### Current behavior

Any nonempty Sylph profile sets:

```python
detected_reads = 1
correct_family_reads = 1
```

The parser does not verify that the expected accession, species, or family was reported.

### Consequence

Any unrelated Sylph hit is treated as a true detection. Sylph logistic-regression and LD50 results are therefore not trustworthy until expected-taxon matching is implemented.

### Expected correction

Use the Sylph accession/contig identifier to map every reported hit onto the same expected taxonomy used for other tools. Define detection at a documented rank. Preserve unrelated hits as false positives rather than true detections.

Add fixtures containing:

- an empty profile;
- the expected accession;
- a different accession in the expected family;
- an unrelated viral accession; and
- malformed output.

## 3. Remove Sylph from read-based negative-control FP percentages

### Location

`scripts/negative/analyze_negative_controls.py`, lines 268-301 and 324-327.

### Current behavior

For Sylph:

```python
fp_reads = dt.height
```

The number of profile rows is divided by input reads and reported as `fp_pct`.

### Consequence

This mixes taxon/reference hits with reads. It is not a read false-positive percentage and is not comparable to read counts from other tools.

### Expected correction

Report Sylph with unit-appropriate endpoints, such as:

- number of false reference hits;
- number of false species/families; and
- optional coverage/abundance burden in Sylph's native units.

Use `NA`, not zero, for read-based endpoints that Sylph cannot provide.

## 4. Audit and reconcile read-pair versus individual-read units

### Relevant locations

- `scripts/mock_benchmark/metrics.py`, lines 30-43
- `scripts/positive/generate_mock_reads.sh`, lines 1-36
- `scripts/positive/analyze_positive_controls.py`, especially lines 62-83, 104-135, 155-185, 205-237, and 327-334
- `scripts/negative/analyze_negative_controls.py`, especially lines 29-43 and 324-327

### Internal assumption

`metrics.py` states that EsViritu reports individual reads while Kraken2, Centrifuger, and Metabuli report paired fragments and require a multiplier of two.

The positive generator uses `wgsim -N` and labels depth as read pairs. The positive script divides tool-reported correct-family counts directly by that pair count. The negative script uses EsViritu readstats as the denominator for all tools without applying the documented multiplier.

### Required HPC validation

For one paired sample per tool, manually reconcile:

- FASTQ records in R1 and R2;
- root/classified counts in native tool output;
- counts in parsed output;
- EsViritu readstats;
- the intended denominator; and
- whether each native report counts pairs, mates, or classified sequences.

Record the evidence in a machine-readable test fixture or provenance note. Do not assume every paired tool has identical semantics solely because `--paired` was used.

### Impact boundary

Binary positive-control detection (`correct_family_reads > 0`) is not affected by a factor-of-two error. Sensitivity/read-recovery plots and negative-control read percentages are affected.

## 5. Stop converting parser failures into biological zeros

### Locations

- `scripts/positive/analyze_positive_controls.py`: repeated `except Exception: pass`
- `scripts/negative/analyze_negative_controls.py`: repeated `except Exception: pass`

### Consequence

- Positive parsing failures become nondetections.
- Negative parsing failures can become apparent perfect specificity.
- Missing files, empty valid outputs, malformed outputs, and failed tool runs are not distinguishable.

### Expected correction

Represent run state explicitly, for example:

- `ok`
- `valid_empty`
- `missing_output`
- `parse_error`
- `tool_failed`

Log file path and exception details. Exclude failed runs from biological denominators and report their count separately. Tests should assert that malformed files raise or produce an explicit failure state rather than zeros.

## 6. Redefine ANI-aware precision/F1 using consistent counting units

### Location

`scripts/mock_benchmark/metrics.py`, lines 202-276.

### Current behavior

- TP/FN are counted per ground-truth genome.
- FP is counted per unique predicted taxon.
- Multiple genomes sharing one expected clade can each become TP from one predicted clade.
- Predictions at ranks with an empty expected set are skipped because of `if not exp_set: continue`.

### Consequence

`TP / (TP + FP_aware)` mixes genome and taxon units and is not conventional precision.

### Expected correction

Choose one of these explicit designs:

1. score unique resolvable clades for TP, FP, and FN; or
2. report genome-level recall separately from taxon-level false-positive burden.

Do not combine genome TPs and unique-taxon FPs in one precision/F1 statistic. Add shared-clade and empty-expected-rank tests.

---

# Priority 2: statistical reanalysis

## 7. Replace or justify the mock-community ordinary linear model

### Location

`notebooks/benchmark_stats.Rmd`, lines 661-671 and 846-872.

### Current behavior

The files are named `mock_lmm_*`, but the model is ordinary `lm`:

```r
metric ~ base_tool * rank + sample
```

Each sample contributes measurements for every tool and rank. A sample fixed effect adjusts block means but does not model residual correlation among repeated tool/rank observations. F1 and Bray-Curtis are bounded.

### Required work

Inspect the number and construction of independent mock communities. Then choose a method matching the actual design, for example:

- sample-blocked permutation tests;
- a mixed/repeated-measures model with suitable random structure;
- cluster-robust standard errors; or
- rank-specific paired analyses with multiplicity control.

Check residuals and sensitivity to bounded-outcome modeling. Rename outputs so they accurately describe the fitted model.

## 8. Remove or correct unpaired, unadjusted plot-level Wilcoxon tests

### Location

`notebooks/benchmark_stats.Rmd`, lines 901-1029.

### Current behavior

`stat_compare_means()` runs four unpaired Wilcoxon comparisons in each facet with no visible multiplicity correction. Across three metrics, four ranks, two truth modes, and four tool comparisons, this can produce up to 96 nominal tests.

### Expected correction

Prefer annotations derived from the authoritative paired model. If direct Wilcoxon tests remain:

- pair by sample explicitly;
- use `paired = TRUE`;
- define the complete family of tests;
- adjust across that family; and
- report adjusted numerical p-values or intervals rather than only stars.

## 9. Clarify the positive-control model and contrasts

### Location

`notebooks/benchmark_stats.Rmd`, lines 130-270.

### Issues to resolve

1. The prose says accession is a random effect, but it is included as a fixed block.
2. `brglmFit` produced an `AS_mixed` bias-reduced estimator in the rendered output; verify whether calling this model "Firth logistic" is technically accurate for the installed `brglm2` version and set the desired estimator explicitly.
3. `emmeans(m, ~ base_tool, at = list(lr_c = ... c(10, 100) ...))` saves one estimate per tool, not separate 10- and 100-read estimates.
4. The model assumes tool effects vary with depth but not with genome/ANI, despite divergence being a central scientific question.

### Expected correction

- State whether inference is conditional on the seven selected genomes or intended to generalize.
- If estimates are wanted at 10 and 100 read pairs, retain depth in the EMM grid/output.
- Examine tool-by-genome or tool-by-ANI heterogeneity.
- Preserve the 10 simulated permutations as independent replicates only if their generation supports that assumption.

## 10. Treat LD50 as an uncertain, censored estimate

### Location

`notebooks/benchmark_stats.Rmd`, lines 243-328.

### Current behavior

- Per-genome fits have no saved convergence diagnostics or uncertainty.
- Never-reaching 50% is plotted at 2,000 reads.
- Curves already above 50% at minimum depth are assigned the minimum grid value.

### Expected correction

- Save model status and effective observations for every fit.
- Derive confidence intervals, preferably via profile/bootstrap appropriate to the model.
- Mark right-censored and left-censored LD50 values explicitly.
- Do not treat censoring display coordinates as measured LD50 values.

## 11. Separate conditional abundance calibration from overall recovery

### Location

`notebooks/benchmark_stats.Rmd`, lines 675-780.

### Current behavior

Calibration filters to detected taxa with positive assigned reads. This estimates calibration conditional on detection, not overall abundance performance.

The pooled stratum also includes genomes where the same clade-level assigned count is repeated across multiple truth genomes. The repository README correctly notes that only `n_truth_sharing_clade == 1` is a true per-genome comparison.

### Expected correction

- Label detected-only calibration as conditional.
- Make the resolvable, single-clade stratum primary for per-genome calibration.
- Add an overall endpoint that penalizes nondetection, such as fixed-denominator error or an explicitly defined zero-aware metric.
- Avoid treating repeated pooled clade totals as independent per-genome measurements.

## 12. Interpret negative-control inference narrowly

### Current committed result

There are seven heterogeneous negative samples. In the committed post-hoc tables:

- Metabuli differs in FP-read percentage from several tools after adjustment.
- No pairwise difference in false species count survives adjustment.

Relevant files:

- `stats/neg_true_negative_fp_pct_posthoc.csv`
- `stats/neg_true_negative_n_species_posthoc.csv`

Do not claim broad superiority in false-species count from the current tests. Recompute these after unit and failure-state corrections.

---

# Priority 3: timing, reproducibility, and provenance

## 13. Make the primary runtime comparison use matched settings

### Location

`notebooks/timing_and_memory.Rmd`, especially lines 66-155.

### Concern

Per-tool summaries mix all samples and all parameter settings. Tools have different numbers of tested settings, and some settings reuse upstream classification/sketch work. The resulting distribution is not a clean comparison of one standard run per tool.

### Expected correction

- Compare one declared default/headline setting per tool as the primary runtime result.
- Show parameter sweeps separately.
- State whether reusable sketch/index/classification costs are charged per run or amortized.
- Report Centrifuger's long-read failure as a failure rather than silently removing it from the comparison.

## 14. Make the repository reproducible

### Current blockers

- `data/` is ignored and acquisition/generation instructions are incomplete.
- `renv.lock` contains only `renv`, not the packages loaded by notebooks.
- There is no Python dependency lock/specification.
- Runner scripts contain HPC-specific absolute paths and named conda environments.
- Exact tool/database/taxonomy versions and checksums are not centrally recorded.
- The top-level `README.md` does not describe setup, data provenance, workflow order, or verification.

### Expected correction

At minimum add:

- a data manifest with paths, schemas, source accessions, generation commands, checksums where feasible, and privacy/licensing notes;
- complete R dependency capture;
- a Python environment specification;
- exact tool versions and database build identifiers/checksums;
- configurable paths instead of embedded `/data/tisza/...` paths;
- an executable workflow or ordered command list; and
- a small public/synthetic fixture dataset that exercises the pipeline without large databases.

Do not commit large or restricted raw data merely to satisfy reproducibility; provide generation/download instructions and manifests where appropriate.

## 15. Add tests before regenerating outputs

No test suite was found during the local audit. Prioritize tests for:

- all five native output parsers;
- empty, missing, and malformed output behavior;
- paired-read unit reconciliation;
- exact-at-node versus clade-cumulative rollup;
- accession version normalization;
- intermediate TaxID ancestry mapping;
- unclassified placeholders;
- ANI interval boundaries and non-enumerated ANI values;
- multiple genomes sharing a resolvable clade;
- capped truth/prediction transformations;
- detection precision/recall/F1 edge cases;
- Bray-Curtis against known vectors; and
- output completeness (expected tool x setting x sample combinations).

---

# Findings that require validation rather than immediate acceptance

These concerns are plausible from the code but need raw HPC evidence before changing results:

1. **Native count semantics:** confirm whether each installed tool reports pairs, mates, fragments, taxon-direct counts, or clade-cumulative counts in every consumed column.
2. **Centrifuger `numReads`:** verify exact-at-node versus cumulative behavior with a small known report.
3. **Taxonomy walk-up:** inspect the actual unmapped/dropped counts and ensure walking to the nearest EsViritu-mapped ancestor cannot create misleading partial lineages.
4. **Fixed versus random accession effects:** accession as a fixed block is not automatically wrong with seven deliberately selected genomes; the scientific generalization target determines the model.
5. **Random seeds:** the primary GLM/LM/emmeans calculations are deterministic. Seeds matter for jittered figures or any new bootstrap/permutation procedure, not as a blanket defect in existing deterministic statistics.
6. **Kraken/Metabuli root counts in negative controls:** root clade count may be exactly the desired total classified-read burden. Do not replace it with root taxon-direct count without checking semantics.

---

# Existing strengths to preserve

- All tools are intended to use the same viral reference universe, reducing database-content confounding.
- The Python mock pipeline is modular across taxonomy, parsing, and metrics.
- Fixed-denominator RPM and within-viral compositional metrics answer distinct useful questions.
- The repository explicitly recognizes that shared clades invalidate naive per-genome assigned-read comparisons.
- Negative-control tests use paired samples and BH adjustment.
- The analysis examines sensitivity, specificity, divergence, abundance, parameter robustness, runtime, and memory rather than relying on one score.

Preserve these design strengths while consolidating implementations.

---

# Suggested staged acceptance criteria

## Stage 1: parser and unit integrity

- Every expected run has an explicit status.
- No parse exception is silently converted to zero.
- Unit semantics are documented with raw-output examples.
- Synthetic parser/unit tests pass.

## Stage 2: mock scoring integrity

- One authoritative scoring implementation exists.
- ANI thresholds and subspecies handling are documented.
- Capped scoring transforms truth and predictions consistently.
- ANI-aware TP/FP/FN use a consistent unit.
- Metric tests pass on hand-calculated examples.

## Stage 3: statistical integrity

- Repeated/paired structure is represented in inference.
- Multiplicity families are declared and adjusted.
- Positive-control depth contrasts are labeled correctly.
- LD50 censoring and uncertainty are explicit.
- Conditional calibration is distinguished from overall recovery.

## Stage 4: reproducibility

- A clean environment can install dependencies.
- A documented command sequence regenerates tables and figures from declared inputs.
- Old and corrected outputs are compared in a review report.
- Large-data requirements and database versions are traceable.

---

# Suggested first prompt on the HPC

> Read `AUDIT_HANDOFF.md` and inspect the current repository and HPC data. Validate Priority 0 and Priority 1 findings against raw outputs and installed tool versions. Produce a plan that separates confirmed code defects, scientific design decisions, and reproducibility work. Do not edit analysis code or regenerate committed outputs until I approve the plan.
