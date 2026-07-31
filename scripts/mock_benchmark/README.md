# Mock-community tool-vs-truth benchmark

Compares taxonomic-profiling tool outputs (esviritu, sylph, kraken2,
centrifuger, metabuli) against the mock-community ground truth, in the
**EsViritu lineage backbone** (`k__ … t__`).

All five tools run against the same EsViritu viral DB, so bacterial/human reads
are unclassified everywhere and the comparison universe is purely viral.

## Pipeline

```bash
source /data/tisza/analyses/mjt_sandbox/conda_mjt.init && conda activate mut_virome

# 1. Parse every tool/param result file into one tidy per-taxon long table
python scripts/benchmark/parse_tool_outputs.py \
    --results-dir benchmarks_mock/results \
    --metadata    esviritu_DBs/v3.2.4b/virus_pathogen_database.all_metadata.tsv \
    --nodes-dmp   esviritu_DBs/v3.2.4/kraken2db/taxonomy/nodes.dmp \
    --out         benchmarks_mock/analysis/unified_long.tsv

# 2. Compute metrics vs ground truth
python scripts/benchmark/compute_benchmark_metrics.py \
    --unified    benchmarks_mock/analysis/unified_long.tsv \
    --truth-dir  mock_virome/mock_community_reads \
    --out-dir    benchmarks_mock/analysis
```

## Modules

| File | Role |
|------|------|
| `taxonomy.py` | Builds Accession/Assembly/NCBI-TaxID → EsViritu lineage maps from the metadata `LineageTaxIDs` column; optional `nodes.dmp` walk-up for intermediate TaxIDs. |
| `parsers.py` | One parser per tool format → unified per-taxon rows (`tool, sample, param_set, <8 ranks>, count, count_type`). Counts are per-taxon (exactly-at-node), not clade-cumulative. |
| `metrics.py` | Bray-Curtis, per-rank detection, ANI-aware detection, abundance accuracy. |
| `parse_tool_outputs.py` | CLI: auto-discovers result files, emits `unified_long.tsv`. |
| `compute_benchmark_metrics.py` | CLI: emits the metrics tables below. |

## Outputs (`benchmarks_mock/analysis/`)

- **`unified_long.tsv`** — every tool/param/sample taxon projected onto the lineage backbone.
- **`rank_metrics.tsv`** — per tool × param × sample × strict rank
  (family/genus/species/subspecies) plus one ANI-capped row (`rank=ani_capped`):
  `bray_curtis_rpm` (**primary**), `bray_curtis_rel`, `truth_total_rpm` /
  `pred_total_rpm` (recovery diagnostics), detection
  (`TP/FP/FN/precision/recall/F1`), abundance (`pearson/spearman/l1/l2`).
- **`ani_aware_detection.tsv`** — per tool × param × sample: each ground-truth
  genome scored at the deepest rank recoverable for its mutation ANI tier.
- **`per_genome_detection.tsv`** — per genome: `expected_rank`, `expected_value`,
  `detected`, `read_count` (truth, individual reads), `assigned_reads` (tool
  reads at that clade, unit-reconciled; `NaN` for sylph), and
  `n_truth_sharing_clade`. Useful for accuracy-vs-ANI plots.
  **Caveat:** taxon-level tools assign reads to a *clade*, not a specific genome,
  so when `n_truth_sharing_clade > 1`, `assigned_reads` is the pooled clade total
  repeated across those genomes (compare it to the *sum* of their `read_count`).
  Only `n_truth_sharing_clade == 1` rows are a true per-genome comparison.

## ANI-aware scoring

Strict rows evaluate every truth genome at each rank. The ANI-capped row
scores each ground-truth genome at the deepest rank it can plausibly be
recovered at, given its mutation ANI (edit `DEFAULT_ANI_TIERS` in `metrics.py`):

| ANI | expected rank |
|-----|---------------|
| 100, 98 | subspecies |
| 95, 90 | species |
| 85, 80 | genus |

ANI-capped scoring requires one truth genome per family. Prediction counts are
then mapped to the cap of the unique matching truth family; unmatched families
remain false-positive taxa. If this invariant is not met, metric calculation
fails rather than assigning shared clade mass ambiguously.

## Bray-Curtis: RPM (primary) vs relative

- **`bray_curtis_rpm` (primary)** uses RPM (reads per million; fixed per-sample
  denominator = ground-truth `total_reads`). Because the denominator is constant
  across truth and all tools, absolute magnitude is preserved, so the metric
  **penalizes both under-classification (low sensitivity) and over-classification
  (false-positive reads)**. There is no "unclassified" bin, so the ~99.9%
  unclassifiable background does not dominate the distance.
- **`bray_curtis_rel`** uses within-viral relative abundance (each vector sums to
  1) and measures compositional *shape* only — useful to decompose a high
  `bray_curtis_rpm` into "wrong composition" vs "wrong total".
- `truth_total_rpm` is constant per sample; compare `pred_total_rpm` to it to
  read off recovery (`pred < truth` = under-classified) or over-classification
  (`pred > truth`).

### Read-unit reconciliation

Ground-truth `read_count` counts **individual reads** (R1+R2). Tool counts are
converted to the same unit via `READ_UNIT_MULTIPLIER` in `metrics.py`:
esviritu ×1 (individual reads); kraken2 / centrifuger / metabuli ×2 (run
`--paired`, so a classified pair counts as 1). **Sylph is excluded from
Bray-Curtis** (`None`) because its `Sequence_abundance` is coverage-based, not a
read count; sylph still gets detection metrics (presence/absence is
unit-independent).

## Notes / conventions

- Per-taxon counts are rolled to coarser ranks by grouping on the lineage value,
  which is correct because counts are exactly-at-node (not clade-cumulative).
  Verified: centrifuger `numReads` is per-taxon (genus-row sum < species-row sum,
  impossible under cumulative counts).
- Taxa above family (root/Viruses/realm) and TaxIDs absent from the EsViritu DB
  are dropped; the parser prints a per-tool dropped-count summary.
- Lead interpretation with the **ANI-aware** detection table: at species rank,
  low-ANI (80/85%) genomes are unavoidable false-negatives by design.
- Underlying-run redundancies to be aware of (not analysis bugs): metabuli
  `ptax10239` ≡ `ptaxnone` (viral-only DB); centrifuger `--min-score` only
  changes results under `mhlauto`; kraken2 `--minimum-hit-groups` only matters at
  `conf0`.
