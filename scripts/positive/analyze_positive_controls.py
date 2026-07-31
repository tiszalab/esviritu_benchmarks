#!/usr/bin/env python3
"""
Analyze positive control benchmarks: detection sensitivity by tool, read count, and bestANI.
"""

import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
from plotnine import *

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from control_read_units import positive_count_to_read_pairs

# ─── Config ───────────────────────────────────────────────────────────────────
BASE_DIR = Path("/data/tisza/analyses/mjt_projects/esviritu_benchmarks/benchmarks_positive1")
RESULTS_DIR = BASE_DIR / "results"
GENOME_INFO = pl.read_csv(BASE_DIR / "genome_info.csv")
OUTPUT_DIR = BASE_DIR / "plots"
OUTPUT_DIR.mkdir(exist_ok=True)

# Map family names to how they appear in each tool's output
FAMILY_ALIASES = {
    "genomoviridae": ["Genomoviridae", "f__Genomoviridae"],
    "picornaviridae": ["Picornaviridae", "f__Picornaviridae"],
}


def parse_sample_name(dirname: str) -> dict | None:
    """Parse sample directory name like OM892290.1.1000.5 into components."""
    m = re.match(r"^([A-Z]{2}\d+\.\d+)\.(\d+)\.(\d+)$", dirname)
    if m:
        return {
            "accession": m.group(1),
            "read_count": int(m.group(2)),
            "permutation": int(m.group(3)),
        }
    return None


# ─── Parse EsViritu results ──────────────────────────────────────────────────
def parse_esviritu(mode: str = "sense"):
    """Parse EsViritu results for a specific mode (sense or sr)."""
    rows = []
    esv_dir = RESULTS_DIR / "esviritu"
    if not esv_dir.exists():
        return pl.DataFrame()

    for sample_dir in sorted(esv_dir.iterdir()):
        if not sample_dir.is_dir():
            continue
        info = parse_sample_name(sample_dir.name)
        if info is None:
            continue

        # Try new naming convention ({sample}_{mode}.detected_virus.info.tsv)
        info_file = sample_dir / f"{sample_dir.name}_{mode}.detected_virus.info.tsv"
        # Fall back to old naming if new doesn't exist
        if not info_file.exists() and mode == "sense":
            info_file = sample_dir / f"{sample_dir.name}.detected_virus.info.tsv"

        detected_reads = 0
        correct_family_reads = 0
        expected_family = GENOME_INFO.filter(
            pl.col("accession") == info["accession"]
        )["family"][0]

        if info_file.exists() and info_file.stat().st_size > 0:
            try:
                dt = pl.read_csv(info_file, separator="\t")
                if dt.height > 0 and "family" in dt.columns and "read_count" in dt.columns:
                    detected_reads = dt["read_count"].sum()
                    aliases = FAMILY_ALIASES.get(expected_family, [])
                    mask = dt["family"].is_in(aliases)
                    correct_family_reads = dt.filter(mask)["read_count"].sum()
            except Exception:
                pass

        rows.append({
            **info,
            "tool": f"esviritu_{mode}",
            "detected_reads": detected_reads,
            "correct_family_reads": correct_family_reads,
        })

    return pl.DataFrame(rows)


# ─── Parse Kraken2 results ───────────────────────────────────────────────────
def parse_kraken2(conf: str = "conf0", mhg: str = "mhg2"):
    """Parse kraken2 results for a specific parameter combo."""
    rows = []
    k2_dir = RESULTS_DIR / "kraken2"
    if not k2_dir.exists():
        return pl.DataFrame()

    for sample_dir in sorted(k2_dir.iterdir()):
        if not sample_dir.is_dir():
            continue
        info = parse_sample_name(sample_dir.name)
        if info is None:
            continue

        report_file = sample_dir / f"{sample_dir.name}_{conf}_{mhg}.report.txt"
        detected_reads = 0
        correct_family_reads = 0
        expected_family = GENOME_INFO.filter(
            pl.col("accession") == info["accession"]
        )["family"][0]

        if report_file.exists() and report_file.stat().st_size > 0:
            try:
                with open(report_file) as f:
                    for line in f:
                        parts = line.strip().split("\t")
                        if len(parts) >= 6:
                            rank = parts[3].strip()
                            reads_rooted = int(parts[1].strip())
                            name = parts[5].strip()
                            if rank == "R" and name == "root":
                                detected_reads = reads_rooted
                            if rank == "F":
                                family_name = name.strip()
                                aliases = FAMILY_ALIASES.get(expected_family, [])
                                if family_name in aliases:
                                    correct_family_reads = reads_rooted
            except Exception:
                pass

        rows.append({
            **info,
            "tool": f"kraken2_{conf}_{mhg}",
            "detected_reads": detected_reads,
            "correct_family_reads": correct_family_reads,
        })

    return pl.DataFrame(rows)


# ─── Parse Centrifuger results ───────────────────────────────────────────────
def parse_centrifuger(mhl: str = "mhl50", ms: str = "ms0"):
    """Parse centrifuger results for a specific parameter combo."""
    rows = []
    cf_dir = RESULTS_DIR / "centrifuger"
    if not cf_dir.exists():
        return pl.DataFrame()

    for sample_dir in sorted(cf_dir.iterdir()):
        if not sample_dir.is_dir():
            continue
        info = parse_sample_name(sample_dir.name)
        if info is None:
            continue

        report_file = sample_dir / f"{sample_dir.name}_{mhl}_{ms}.report.tsv"
        detected_reads = 0
        correct_family_reads = 0
        expected_family = GENOME_INFO.filter(
            pl.col("accession") == info["accession"]
        )["family"][0]

        if report_file.exists() and report_file.stat().st_size > 0:
            try:
                dt = pl.read_csv(report_file, separator="\t")
                if dt.height > 0:
                    # root row has total reads
                    root_rows = dt.filter(pl.col("name") == "root")
                    if root_rows.height > 0:
                        detected_reads = root_rows["numReads"][0]
                    # family row
                    aliases = FAMILY_ALIASES.get(expected_family, [])
                    family_rows = dt.filter(
                        (pl.col("taxRank") == "family") & pl.col("name").is_in(aliases)
                    )
                    if family_rows.height > 0:
                        correct_family_reads = family_rows["numReads"][0]
            except Exception:
                pass

        rows.append({
            **info,
            "tool": f"centrifuger_{mhl}_{ms}",
            "detected_reads": detected_reads,
            "correct_family_reads": correct_family_reads,
        })

    return pl.DataFrame(rows)


# ─── Parse Metabuli results ──────────────────────────────────────────────────
def parse_metabuli(prec: str = "prec0", ms: str = "ms0", ptax: str = "ptax10239"):
    """Parse metabuli results for a specific parameter combo."""
    rows = []
    mb_dir = RESULTS_DIR / "metabuli"
    if not mb_dir.exists():
        return pl.DataFrame()

    for sample_dir in sorted(mb_dir.iterdir()):
        if not sample_dir.is_dir():
            continue
        info = parse_sample_name(sample_dir.name)
        if info is None:
            continue

        setting_dir = sample_dir / f"{sample_dir.name}_{prec}_{ms}_{ptax}"
        report_file = setting_dir / f"{sample_dir.name}_{prec}_{ms}_{ptax}_report.tsv"
        detected_reads = 0
        correct_family_reads = 0
        expected_family = GENOME_INFO.filter(
            pl.col("accession") == info["accession"]
        )["family"][0]

        if report_file.exists() and report_file.stat().st_size > 0:
            try:
                with open(report_file) as f:
                    for line in f:
                        if line.startswith("#"):
                            continue
                        parts = line.strip().split("\t")
                        if len(parts) >= 6:
                            rank = parts[3].strip()
                            reads_rooted = int(parts[1].strip())
                            name = parts[5].strip()
                            if rank == "no rank" and name == "root":
                                detected_reads = reads_rooted
                            if rank == "family":
                                aliases = FAMILY_ALIASES.get(expected_family, [])
                                if name.strip() in aliases:
                                    correct_family_reads = reads_rooted
            except Exception:
                pass

        rows.append({
            **info,
            "tool": f"metabuli_{prec}_{ms}_{ptax}",
            "detected_reads": detected_reads,
            "correct_family_reads": correct_family_reads,
        })

    return pl.DataFrame(rows)


# ─── Parse Sylph results ─────────────────────────────────────────────────────
def parse_sylph(mnk: str = "mnk10"):
    """Parse sylph results for a specific min_n_kmers setting."""
    rows = []
    sylph_dir = RESULTS_DIR / "sylph"
    if not sylph_dir.exists():
        return pl.DataFrame()

    for sample_dir in sorted(sylph_dir.iterdir()):
        if not sample_dir.is_dir():
            continue
        info = parse_sample_name(sample_dir.name)
        if info is None:
            continue

        profile_file = sample_dir / f"{sample_dir.name}_{mnk}.profile.tsv"
        detected_reads = 0
        correct_family_reads = 0

        if profile_file.exists() and profile_file.stat().st_size > 0:
            try:
                dt = pl.read_csv(profile_file, separator="\t")
                if dt.height > 0:
                    # Sylph doesn't report read counts; use 1 as indicator of detection
                    detected_reads = 1
                    correct_family_reads = 1
            except Exception:
                pass

        rows.append({
            **info,
            "tool": f"sylph_{mnk}",
            "detected_reads": detected_reads,
            "correct_family_reads": correct_family_reads,
        })

    return pl.DataFrame(rows)


# ─── Main ────────────────────────────────────────────────────────────────────
def main():
    print("Parsing EsViritu results...")
    esv_dfs = [parse_esviritu(mode) for mode in ["sense", "sr"]]
    esv_df = pl.concat([df for df in esv_dfs if df.height > 0])

    print("Parsing Kraken2 results...")
    k2_dfs = []
    for conf in ["conf0", "conf0.2", "conf0.5"]:
        for mhg in ["mhg2", "mhg3", "mhg5"]:
            k2_dfs.append(parse_kraken2(conf, mhg))
    k2_df = pl.concat([df for df in k2_dfs if df.height > 0])

    print("Parsing Centrifuger results...")
    cf_dfs = []
    for mhl in ["mhl50", "mhl100", "mhlauto"]:
        for ms in ["ms0", "ms100", "ms500"]:
            cf_dfs.append(parse_centrifuger(mhl, ms))
    cf_df = pl.concat([df for df in cf_dfs if df.height > 0])

    print("Parsing Metabuli results...")
    mb_dfs = []
    for prec in ["prec0", "prec1"]:
        for ms in ["ms0", "ms0.15"]:
            for ptax in ["ptax10239", "ptaxnone"]:
                mb_dfs.append(parse_metabuli(prec, ms, ptax))
    mb_df = pl.concat([df for df in mb_dfs if df.height > 0])

    print("Parsing Sylph results...")
    sylph_dfs = []
    for mnk in ["mnk10", "mnk20", "mnk50"]:
        sylph_dfs.append(parse_sylph(mnk))
    sylph_df = pl.concat([df for df in sylph_dfs if df.height > 0])

    # Combine all results
    all_dfs = [df for df in [esv_df, k2_df, cf_df, mb_df, sylph_df] if df.height > 0]
    if not all_dfs:
        print("No results found!")
        return

    combined = pl.concat(all_dfs)

    # Join genome info
    combined = combined.join(GENOME_INFO, on="accession", how="left")

    # Compute detection metric: fraction of input read pairs classified to correct family
    combined = combined.with_columns(
        pl.col("read_count").alias("input_read_pairs"),
        pl.struct("tool", "detected_reads").map_elements(
            lambda x: positive_count_to_read_pairs(x["tool"], x["detected_reads"]),
            return_dtype=pl.Float64,
        ).alias("detected_read_pairs"),
        pl.struct("tool", "correct_family_reads").map_elements(
            lambda x: positive_count_to_read_pairs(x["tool"], x["correct_family_reads"]),
            return_dtype=pl.Float64,
        ).alias("correct_family_read_pairs"),
        (pl.col("correct_family_reads") > 0).cast(pl.Int8).alias("detected"),
    ).with_columns(
        (pl.col("correct_family_read_pairs") / pl.col("input_read_pairs")).alias("sensitivity"),
    )

    # Save combined table
    combined.write_csv(BASE_DIR / "positive_control_results.csv")
    print(f"Combined results: {combined.height} rows")
    print(combined.group_by("tool").agg(pl.col("detected").sum(), pl.len().alias("n")))

    # ─── Plots ────────────────────────────────────────────────────────────────

    # Use default parameter combos for a cleaner comparison plot
    plot_tools = ["esviritu_sense", "esviritu_sr", "kraken2_conf0_mhg2", "centrifuger_mhl50_ms0", "metabuli_prec0_ms0_ptax10239", "sylph_mnk10"]
    plot_pl = combined.filter(pl.col("tool").is_in(plot_tools))

    # Rename tools for display
    tool_labels = {
        "esviritu_sense": "EsViritu\n(sense)",
        "esviritu_sr": "EsViritu\n(sr)",
        "kraken2_conf0_mhg2": "Kraken2\n(conf0, mhg2)",
        "centrifuger_mhl50_ms0": "Centrifuger\n(mhl50, ms0)",
        "metabuli_prec0_ms0_ptax10239": "Metabuli\n(prec0, ms0)",
        "sylph_mnk10": "Sylph\n(mnk10)",
    }
    plot_pl = plot_pl.with_columns(
        pl.col("tool").replace(tool_labels).alias("tool_label")
    )
    plot_df = pd.DataFrame(plot_pl.to_dict())

    # Detection rate by read count and bestANI
    detection_summary = (
        plot_df.groupby(["tool_label", "read_count", "accession", "family", "bestANI"])
        .agg(detection_rate=("detected", "mean"), n=("detected", "count"))
        .reset_index()
    )

    # Plot 1: Detection rate vs read count, faceted by genome, colored by tool
    p1 = (
        ggplot(detection_summary, aes(x="factor(read_count)", y="detection_rate", fill="tool_label"))
        + geom_col(position="dodge")
        + facet_wrap("~accession", ncol=4)
        + labs(
            x="Read pairs (input)",
            y="Detection rate (fraction of permutations)",
            fill="Tool",
            title="Detection rate by read count and genome",
        )
        + theme_minimal()
        + theme(
            axis_text_x=element_text(rotation=45, hjust=1),
            figure_size=(12, 6),
        )
        + scale_y_continuous(limits=(0, 1))
    )
    p1.save(OUTPUT_DIR / "detection_rate_by_readcount.pdf", dpi=150)
    p1.save(OUTPUT_DIR / "detection_rate_by_readcount.png", dpi=150)
    print("Saved: detection_rate_by_readcount")

    # Plot 2: Detection rate vs bestANI, colored by tool, faceted by read count
    p2 = (
        ggplot(detection_summary, aes(x="bestANI", y="detection_rate", color="tool_label"))
        + geom_point(size=3, alpha=0.7)
        + geom_line(aes(group="tool_label"), alpha=0.5)
        + facet_wrap("~read_count", ncol=4, labeller="label_both")
        + labs(
            x="Best ANI to reference (%)",
            y="Detection rate",
            color="Tool",
            title="Detection rate vs. ANI to closest reference genome",
        )
        + theme_minimal()
        + theme(figure_size=(14, 5))
        + scale_y_continuous(limits=(0, 1))
    )
    p2.save(OUTPUT_DIR / "detection_rate_vs_bestANI.pdf", dpi=150)
    p2.save(OUTPUT_DIR / "detection_rate_vs_bestANI.png", dpi=150)
    print("Saved: detection_rate_vs_bestANI")

    # Plot 3: Sensitivity (fraction of reads recovered) - box/jitter
    p3 = (
        ggplot(plot_df, aes(x="factor(read_count)", y="sensitivity", color="tool_label"))
        + geom_jitter(width=0.2, alpha=0.4, size=1.5)
        + stat_summary(fun_y=np.mean, geom="point", size=4, shape="_")
        + facet_wrap("~accession", ncol=4)
        + labs(
            x="Read pairs (input)",
            y="Sensitivity (correct family reads / input reads)",
            color="Tool",
            title="Read recovery sensitivity by genome and read count",
        )
        + theme_minimal()
        + theme(
            axis_text_x=element_text(rotation=45, hjust=1),
            figure_size=(14, 7),
        )
    )
    p3.save(OUTPUT_DIR / "sensitivity_by_genome.pdf", dpi=150)
    p3.save(OUTPUT_DIR / "sensitivity_by_genome.png", dpi=150)
    print("Saved: sensitivity_by_genome")

    # Plot 4: Summary heatmap-like: detection rate by bestANI group and read count
    summary_for_heatmap = (
        plot_df.groupby(["tool_label", "read_count", "family", "bestANI"])
        .agg(detection_rate=("detected", "mean"))
        .reset_index()
    )
    summary_for_heatmap["label"] = (
        summary_for_heatmap["family"] + "\n(ANI=" + summary_for_heatmap["bestANI"].astype(str) + ")"
    )

    p4 = (
        ggplot(summary_for_heatmap, aes(x="factor(read_count)", y="label", fill="detection_rate"))
        + geom_tile(color="white", size=0.5)
        + geom_text(aes(label="detection_rate"), format_string="{:.1f}", size=8)
        + facet_wrap("~tool_label", ncol=3)
        + scale_fill_gradient(low="white", high="#2166AC", limits=(0, 1))
        + labs(
            x="Read pairs (input)",
            y="Genome (family, bestANI)",
            fill="Detection\nrate",
            title="Detection rate heatmap by tool, read depth, and ANI",
        )
        + theme_minimal()
        + theme(figure_size=(12, 6))
    )
    p4.save(OUTPUT_DIR / "detection_heatmap.pdf", dpi=150)
    p4.save(OUTPUT_DIR / "detection_heatmap.png", dpi=150)
    print("Saved: detection_heatmap")

    print("\nDone! All outputs in:", OUTPUT_DIR)


if __name__ == "__main__":
    main()
