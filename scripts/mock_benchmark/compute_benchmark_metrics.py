#!/usr/bin/env python3
"""
compute_benchmark_metrics.py

Consume the unified long table (from parse_tool_outputs.py) plus the per-sample
ground-truth tables and emit benchmark metrics TSVs.

Outputs (written under --out-dir):
  rank_metrics.tsv        one row per tool x param_set x sample x rank:
                          bray_curtis_rpm (PRIMARY; RPM-normalized, captures
                          sensitivity & over-classification), bray_curtis_rel
                          (compositional), truth_total_rpm / pred_total_rpm
                          (recovery diagnostics), detection
                          (TP/FP/FN/precision/recall/F1), abundance
                          (pearson/spearman/l1/l2).
                          Read units are reconciled to individual reads
                          (paired tools x2); sylph is excluded from Bray-Curtis
                          (its output is coverage-based, not a read count).
  ani_aware_detection.tsv one row per tool x param_set x sample:
                          ANI-aware TP/FN/FP_aware/precision/recall/F1
  per_genome_detection.tsv one row per tool x param_set x sample x ground-truth
                          genome: expected_rank, expected_value, detected

Usage
-----
python scripts/benchmark/compute_benchmark_metrics.py \
    --unified    benchmarks_mock/analysis/unified_long.tsv \
    --truth-dir  mock_virome/mock_community_reads \
    --out-dir    benchmarks_mock/analysis
"""

import argparse
import glob
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from taxonomy import ESV_RANKS                          # noqa: E402
import metrics as M                                     # noqa: E402


def find_truth_file(truth_dir: str, sample: str):
    hits = sorted(glob.glob(os.path.join(truth_dir, f"{sample}.*ground_truth.tsv")))
    if not hits:
        hits = sorted(glob.glob(os.path.join(truth_dir, f"{sample}*ground_truth.tsv")))
    return hits[0] if hits else None


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--unified", required=True)
    ap.add_argument("--truth-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    return ap.parse_args()


def main():
    args = parse_args()
    long_df = pd.read_csv(args.unified, sep="\t", dtype=str).fillna("")
    long_df["count"] = pd.to_numeric(long_df["count"], errors="coerce").fillna(0.0)

    samples = sorted(long_df["sample"].unique())
    truth_cache = {}
    for s in samples:
        tf = find_truth_file(args.truth_dir, s)
        if tf is None:
            print(f"WARNING: no ground truth for sample {s}", file=sys.stderr)
            continue
        truth_cache[s] = pd.read_csv(tf, sep="\t", dtype=str).fillna("")
        print(f"Truth {s}: {tf}", file=sys.stderr)

    rank_rows, ani_rows, genome_rows = [], [], []

    keys = long_df[["tool", "param_set", "sample"]].drop_duplicates()
    for tool, param, sample in keys.itertuples(index=False):
        if sample not in truth_cache:
            continue
        pred = long_df[(long_df["tool"] == tool) &
                       (long_df["param_set"] == param) &
                       (long_df["sample"] == sample)]
        truth = truth_cache[sample]
        total_reads = pd.to_numeric(truth["total_reads"],
                                    errors="coerce").dropna()
        total_reads = float(total_reads.iloc[0]) if len(total_reads) else float("nan")

        rrows, arow, grows = M.compute_community_and_detection(
            truth, pred, tool, total_reads)
        meta = {"tool": tool, "param_set": param, "sample": sample}
        for r in rrows:
            rank_rows.append({**meta, **r})
        ani_rows.append({**meta, **arow})
        for g in grows:
            genome_rows.append({**meta, **g})

    os.makedirs(args.out_dir, exist_ok=True)
    rank_df = pd.DataFrame(rank_rows)
    ani_df = pd.DataFrame(ani_rows)
    genome_df = pd.DataFrame(genome_rows)

    rank_path = os.path.join(args.out_dir, "rank_metrics.tsv")
    ani_path = os.path.join(args.out_dir, "ani_aware_detection.tsv")
    genome_path = os.path.join(args.out_dir, "per_genome_detection.tsv")

    rank_df.to_csv(rank_path, sep="\t", index=False)
    ani_df.to_csv(ani_path, sep="\t", index=False)
    genome_df.to_csv(genome_path, sep="\t", index=False)

    print(f"Wrote {rank_path} ({len(rank_df)} rows)", file=sys.stderr)
    print(f"Wrote {ani_path} ({len(ani_df)} rows)", file=sys.stderr)
    print(f"Wrote {genome_path} ({len(genome_df)} rows)", file=sys.stderr)

    # console summary: best param per tool by mean species-level RPM Bray-Curtis
    # (sylph is excluded from Bray-Curtis -> NaN, dropped here)
    if not rank_df.empty:
        sp = rank_df[rank_df["rank"] == "species"]
        best = (sp.groupby(["tool", "param_set"])["bray_curtis_rpm"]
                .mean().reset_index().dropna(subset=["bray_curtis_rpm"])
                .sort_values(["tool", "bray_curtis_rpm"]))
        print("\nMean species-level RPM Bray-Curtis by tool/param "
              "(lower = closer to truth):", file=sys.stderr)
        print(best.to_string(index=False), file=sys.stderr)


if __name__ == "__main__":
    main()
