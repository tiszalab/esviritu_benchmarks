#!/usr/bin/env python3
"""
Analyze negative control results across all tools and parameter settings.
Metrics: false-positive reads classified, number of species called, families called.
"""

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
from plotnine import *

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from control_read_units import negative_count_to_individual_reads

# ─── Config ───────────────────────────────────────────────────────────────────
BASE_DIR = Path("/data/tisza/analyses/mjt_projects/esviritu_benchmarks/benchmarks_negative1")
RESULTS_DIR = BASE_DIR / "results"
SRR_INFO = pl.read_csv(BASE_DIR / "srr_info.csv")
# Strip whitespace from columns
SRR_INFO = SRR_INFO.with_columns(
    pl.col("source_type").str.strip_chars(),
    pl.col("read_type").str.strip_chars(),
)
OUTPUT_DIR = BASE_DIR / "plots"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

SAMPLES = [s for s in SRR_INFO["accession"].to_list() if s.strip()]

# ─── Helpers ──────────────────────────────────────────────────────────────────
def get_total_reads(sample: str, mode: str = "sense") -> int:
    """Get total read count from EsViritu readstats."""
    # Try new naming convention first
    yaml_file = RESULTS_DIR / "esviritu" / sample / f"{sample}_{mode}_esviritu.readstats.yaml"
    if not yaml_file.exists():
        # Fall back to old naming
        yaml_file = RESULTS_DIR / "esviritu" / sample / f"{sample}_esviritu.readstats.yaml"
    if yaml_file.exists():
        with open(yaml_file) as f:
            text = f.read()
            m = re.search(r"(\d+)", text)
            if m:
                return int(m.group(1))
    return 0


# ─── Parse EsViritu ──────────────────────────────────────────────────────────
def parse_esviritu():
    rows = []
    esv_dir = RESULTS_DIR / "esviritu"
    for sample in SAMPLES:
        sample_dir = esv_dir / sample
        if not sample_dir.exists():
            continue

        for mode in ["sense", "sr"]:
            total_reads = get_total_reads(sample, mode)

            # Try new naming convention first, fall back to old for sense
            tax_file = sample_dir / f"{sample}_{mode}.tax_profile.tsv"
            if not tax_file.exists() and mode == "sense":
                tax_file = sample_dir / f"{sample}.tax_profile.tsv"
            info_file = sample_dir / f"{sample}_{mode}.detected_virus.info.tsv"
            if not info_file.exists() and mode == "sense":
                info_file = sample_dir / f"{sample}.detected_virus.info.tsv"

            fp_reads = 0
            n_species = 0
            n_families = 0
            species_list = []

            if tax_file.exists() and tax_file.stat().st_size > 0:
                try:
                    dt = pl.read_csv(tax_file, separator="\t")
                    if dt.height > 0:
                        fp_reads = dt["read_count"].sum()
                        n_species = dt.height
                        n_families = dt["family"].n_unique()
                        species_list = dt["species"].to_list()
                except Exception:
                    pass
            elif info_file.exists() and info_file.stat().st_size > 0:
                try:
                    dt = pl.read_csv(info_file, separator="\t")
                    if dt.height > 0:
                        fp_reads = dt["read_count"].sum()
                        n_species = dt["species"].n_unique() if "species" in dt.columns else dt.height
                        n_families = dt["family"].n_unique() if "family" in dt.columns else 0
                except Exception:
                    pass

            rows.append({
                "sample": sample,
                "tool": f"EsViritu ({mode})",
                "settings": mode,
                "total_reads": total_reads,
                "fp_reads": fp_reads,
                "n_species": n_species,
                "n_families": n_families,
                "species_list": "; ".join(str(s) for s in species_list),
            })
    return rows


# ─── Parse Kraken2 ────────────────────────────────────────────────────────────
def parse_kraken2():
    rows = []
    k2_dir = RESULTS_DIR / "kraken2"
    confs = ["conf0", "conf0.2", "conf0.5"]
    mhgs = ["mhg2", "mhg3", "mhg5"]

    for sample in SAMPLES:
        sample_dir = k2_dir / sample
        if not sample_dir.exists():
            continue
        total_reads = get_total_reads(sample)

        for conf in confs:
            for mhg in mhgs:
                report_file = sample_dir / f"{sample}_{conf}_{mhg}.report.txt"
                if not report_file.exists():
                    continue

                fp_reads = 0
                species_list = []
                families = set()

                try:
                    with open(report_file) as f:
                        for line in f:
                            parts = line.strip().split("\t")
                            if len(parts) >= 6:
                                rank = parts[3].strip()
                                clade_reads = int(parts[1].strip())
                                name = parts[5].strip()
                                if rank == "R" and name == "root":
                                    fp_reads = clade_reads
                                if rank == "S" and clade_reads > 0:
                                    species_list.append(name)
                                if rank == "F" and clade_reads > 0:
                                    families.add(name)
                except Exception:
                    pass

                rows.append({
                    "sample": sample,
                    "tool": "Kraken2",
                    "settings": f"{conf}_{mhg}",
                    "total_reads": total_reads,
                    "fp_reads": fp_reads,
                    "n_species": len(species_list),
                    "n_families": len(families),
                    "species_list": "; ".join(species_list[:20]),
                })
    return rows


# ─── Parse Centrifuger ────────────────────────────────────────────────────────
def parse_centrifuger():
    rows = []
    cf_dir = RESULTS_DIR / "centrifuger"
    mhls = ["mhl50", "mhl100", "mhlauto"]
    mss = ["ms0", "ms100", "ms500"]

    for sample in SAMPLES:
        sample_dir = cf_dir / sample
        if not sample_dir.exists():
            continue
        total_reads = get_total_reads(sample)

        for mhl in mhls:
            for ms in mss:
                report_file = sample_dir / f"{sample}_{mhl}_{ms}.report.tsv"
                if not report_file.exists():
                    continue

                fp_reads = 0
                species_list = []
                families = set()

                try:
                    dt = pl.read_csv(report_file, separator="\t")
                    root_rows = dt.filter(pl.col("name") == "root")
                    if root_rows.height > 0:
                        fp_reads = root_rows["numReads"][0]
                    sp_rows = dt.filter(
                        (pl.col("taxRank") == "species") & (pl.col("numReads") > 0)
                    )
                    species_list = sp_rows["name"].to_list()
                    fam_rows = dt.filter(
                        (pl.col("taxRank") == "family") & (pl.col("numReads") > 0)
                    )
                    families = set(fam_rows["name"].to_list())
                except Exception:
                    pass

                rows.append({
                    "sample": sample,
                    "tool": "Centrifuger",
                    "settings": f"{mhl}_{ms}",
                    "total_reads": total_reads,
                    "fp_reads": fp_reads,
                    "n_species": len(species_list),
                    "n_families": len(families),
                    "species_list": "; ".join(str(s).replace("_", " ") for s in species_list[:20]),
                })
    return rows


# ─── Parse Metabuli ───────────────────────────────────────────────────────────
def parse_metabuli():
    rows = []
    mb_dir = RESULTS_DIR / "metabuli"
    precs = ["prec0", "prec1"]
    mss = ["ms0", "ms0.15"]
    ptaxs = ["ptax10239", "ptaxnone"]

    for sample in SAMPLES:
        sample_dir = mb_dir / sample
        if not sample_dir.exists():
            continue
        total_reads = get_total_reads(sample)

        for prec in precs:
            for ms in mss:
                for ptax in ptaxs:
                    setting_dir = sample_dir / f"{sample}_{prec}_{ms}_{ptax}"
                    report_file = setting_dir / f"{sample}_{prec}_{ms}_{ptax}_report.tsv"
                    if not report_file.exists():
                        continue

                    fp_reads = 0
                    species_list = []
                    families = set()

                    try:
                        with open(report_file) as f:
                            for line in f:
                                if line.startswith("#"):
                                    continue
                                parts = line.strip().split("\t")
                                if len(parts) >= 6:
                                    rank = parts[3].strip()
                                    clade_reads = int(parts[1].strip())
                                    name = parts[5].strip()
                                    if rank == "no rank" and name == "root":
                                        fp_reads = clade_reads
                                    if rank == "species" and clade_reads > 0:
                                        species_list.append(name)
                                    if rank == "family" and clade_reads > 0:
                                        families.add(name)
                    except Exception:
                        pass

                    rows.append({
                        "sample": sample,
                        "tool": "Metabuli",
                        "settings": f"{prec}_{ms}_{ptax}",
                        "total_reads": total_reads,
                        "fp_reads": fp_reads,
                        "n_species": len(species_list),
                        "n_families": len(families),
                        "species_list": "; ".join(species_list[:20]),
                    })
    return rows


# ─── Parse Sylph ──────────────────────────────────────────────────────────────
def parse_sylph():
    rows = []
    sylph_dir = RESULTS_DIR / "sylph"
    mnks = ["mnk10", "mnk20", "mnk50"]

    for sample in SAMPLES:
        sample_dir = sylph_dir / sample
        if not sample_dir.exists():
            continue
        total_reads = get_total_reads(sample)

        for mnk in mnks:
            profile_file = sample_dir / f"{sample}_{mnk}.profile.tsv"
            fp_reads = 0
            n_species = 0
            species_list = []

            if profile_file.exists() and profile_file.stat().st_size > 0:
                try:
                    dt = pl.read_csv(profile_file, separator="\t")
                    if dt.height > 0:
                        n_species = dt.height
                        species_list = dt["Contig_name"].to_list()
                        fp_reads = dt.height  # no read count; use hit count
                except Exception:
                    pass

            rows.append({
                "sample": sample,
                "tool": "Sylph",
                "settings": mnk,
                "total_reads": total_reads,
                "fp_reads": fp_reads,
                "n_species": n_species,
                "n_families": 0,  # sylph doesn't report family
                "species_list": "; ".join(str(s) for s in species_list[:20]),
            })
    return rows


# ─── Main ─────────────────────────────────────────────────────────────────────
def main():
    print("Parsing negative control results...")

    all_rows = []
    for name, parser in [
        ("EsViritu", parse_esviritu),
        ("Kraken2", parse_kraken2),
        ("Centrifuger", parse_centrifuger),
        ("Metabuli", parse_metabuli),
        ("Sylph", parse_sylph),
    ]:
        r = parser()
        all_rows.extend(r)
        print(f"  {name}: {len(r)} entries")

    combined = pl.DataFrame(all_rows).rename({"fp_reads": "fp_reads_native"})
    combined = combined.join(
        SRR_INFO.rename({"accession": "sample"}), on="sample", how="left"
    ).with_columns(
        pl.struct("tool", "fp_reads_native", "read_type").map_elements(
            lambda x: negative_count_to_individual_reads(
                x["tool"], x["fp_reads_native"], x["read_type"]
            ),
            return_dtype=pl.Float64,
        ).alias("fp_reads"),
    ).with_columns(
        (pl.col("fp_reads") / pl.col("total_reads") * 100).alias("fp_pct"),
    )

    combined.write_csv(BASE_DIR / "negative_control_results.csv")
    print(f"\nTotal rows: {combined.height}")

    plot_df = pd.DataFrame(combined.to_dict())
    plot_df["sample_label"] = plot_df["sample"] + "\n(" + plot_df["source_type"] + ")"

    # ─── Plot 1: FP reads (%) by tool and settings ────────────────────────────
    p1 = (
        ggplot(plot_df, aes(x="sample_label", y="fp_pct", fill="settings"))
        + geom_col(position="dodge")
        + geom_text(
            aes(label="fp_pct"),
            position=position_dodge(width=0.9),
            angle=90, va="bottom", size=5,
            format_string="{:.4f}",
        )
        + facet_wrap("~tool", ncol=2)
        + labs(
            x="Sample",
            y="False-positive reads (%)",
            fill="Settings",
            title="Negative controls: Reads classified as viral (%)",
        )
        + theme_minimal()
        + theme(
            axis_text_x=element_text(rotation=45, hjust=1, size=7),
            figure_size=(14, 12),
            legend_text=element_text(size=7),
            legend_key_size=10,
        )
        + scale_y_continuous(labels=lambda l: [f"{v:.3f}" for v in l])
    )
    p1.save(OUTPUT_DIR / "neg_ctrl_fp_reads_pct.pdf", dpi=150)
    p1.save(OUTPUT_DIR / "neg_ctrl_fp_reads_pct.png", dpi=150)
    print("Saved: neg_ctrl_fp_reads_pct")

    # ─── Plot 2: Number of species called ─────────────────────────────────────
    p2 = (
        ggplot(plot_df, aes(x="sample_label", y="n_species", fill="settings"))
        + geom_col(position="dodge")
        + geom_text(
            aes(label="n_species"),
            position=position_dodge(width=0.9),
            angle=90, va="bottom", size=5,
            format_string="{:.0f}",
        )
        + facet_wrap("~tool", ncol=2)
        + labs(
            x="Sample",
            y="Number of species called",
            fill="Settings",
            title="Negative controls: Number of false-positive species calls",
        )
        + theme_minimal()
        + theme(
            axis_text_x=element_text(rotation=45, hjust=1, size=7),
            figure_size=(14, 12),
            legend_text=element_text(size=7),
            legend_key_size=10,
        )
    )
    p2.save(OUTPUT_DIR / "neg_ctrl_n_species.pdf", dpi=150)
    p2.save(OUTPUT_DIR / "neg_ctrl_n_species.png", dpi=150)
    print("Saved: neg_ctrl_n_species")

    # ─── Plot 3: Number of families called ────────────────────────────────────
    # Exclude sylph (no family info)
    fam_df = plot_df[plot_df["tool"] != "Sylph"]
    p3 = (
        ggplot(fam_df, aes(x="sample_label", y="n_families", fill="settings"))
        + geom_col(position="dodge")
        + geom_text(
            aes(label="n_families"),
            position=position_dodge(width=0.9),
            angle=90, va="bottom", size=5,
            format_string="{:.0f}",
        )
        + facet_wrap("~tool", ncol=2)
        + labs(
            x="Sample",
            y="Number of families called",
            fill="Settings",
            title="Negative controls: Number of false-positive family calls",
        )
        + theme_minimal()
        + theme(
            axis_text_x=element_text(rotation=45, hjust=1, size=7),
            figure_size=(14, 10),
            legend_text=element_text(size=7),
            legend_key_size=10,
        )
    )
    p3.save(OUTPUT_DIR / "neg_ctrl_n_families.pdf", dpi=150)
    p3.save(OUTPUT_DIR / "neg_ctrl_n_families.png", dpi=150)
    print("Saved: neg_ctrl_n_families")

    # ─── Plot 4: Cross-tool comparison (default settings) ─────────────────────
    default_settings = {
        "EsViritu (sense)": "sense",
        "EsViritu (sr)": "sr",
        "Kraken2": "conf0_mhg2",
        "Centrifuger": "mhl50_ms0",
        "Metabuli": "prec0_ms0_ptax10239",
        "Sylph": "mnk10",
    }
    default_df = plot_df[
        plot_df.apply(lambda r: r["settings"] == default_settings.get(r["tool"], ""), axis=1)
    ].copy()

    p4 = (
        ggplot(default_df, aes(x="sample_label", y="fp_pct", fill="tool"))
        + geom_col(position="dodge")
        + geom_text(
            aes(label="fp_pct"),
            position=position_dodge(width=0.9),
            angle=90, va="bottom", size=6,
            format_string="{:.4f}",
        )
        + labs(
            x="Sample",
            y="False-positive reads (%)",
            fill="Tool",
            title="Negative controls: FP rate comparison (default settings)",
        )
        + theme_minimal()
        + theme(
            axis_text_x=element_text(rotation=45, hjust=1),
            figure_size=(12, 6),
        )
        + scale_y_continuous(labels=lambda l: [f"{v:.3f}" for v in l])
    )
    p4.save(OUTPUT_DIR / "neg_ctrl_fp_comparison.pdf", dpi=150)
    p4.save(OUTPUT_DIR / "neg_ctrl_fp_comparison.png", dpi=150)
    print("Saved: neg_ctrl_fp_comparison")

    # ─── Plot 5: Cross-tool species count comparison (default settings) ───────
    p5 = (
        ggplot(default_df, aes(x="sample_label", y="n_species", fill="tool"))
        + geom_col(position="dodge")
        + geom_text(
            aes(label="n_species"),
            position=position_dodge(width=0.9),
            angle=90, va="bottom", size=6,
            format_string="{:.0f}",
        )
        + labs(
            x="Sample",
            y="Number of species called",
            fill="Tool",
            title="Negative controls: Species count comparison (default settings)",
        )
        + theme_minimal()
        + theme(
            axis_text_x=element_text(rotation=45, hjust=1),
            figure_size=(12, 6),
        )
    )
    p5.save(OUTPUT_DIR / "neg_ctrl_species_comparison.pdf", dpi=150)
    p5.save(OUTPUT_DIR / "neg_ctrl_species_comparison.png", dpi=150)
    print("Saved: neg_ctrl_species_comparison")

    # ─── Plot 6: Effect of parameters on FP rate per tool ─────────────────────
    for tool in ["Kraken2", "Centrifuger", "Metabuli"]:
        tool_df = plot_df[plot_df["tool"] == tool].copy()
        if len(tool_df) == 0:
            continue

        p6 = (
            ggplot(tool_df, aes(x="settings", y="fp_pct", fill="sample_label"))
            + geom_col(position="dodge")
            + geom_text(
                aes(label="fp_pct"),
                position=position_dodge(width=0.9),
                angle=90, va="bottom", size=5,
                format_string="{:.4f}",
            )
            + labs(
                x="Parameter settings",
                y="False-positive reads (%)",
                fill="Sample",
                title=f"{tool}: Effect of parameters on FP rate",
            )
            + theme_minimal()
            + theme(
                axis_text_x=element_text(rotation=45, hjust=1, size=7),
                figure_size=(12, 6),
            )
            + scale_y_continuous(labels=lambda l: [f"{v:.4f}" for v in l])
        )
        p6.save(OUTPUT_DIR / f"neg_ctrl_{tool.lower()}_params_fp.pdf", dpi=150)
        p6.save(OUTPUT_DIR / f"neg_ctrl_{tool.lower()}_params_fp.png", dpi=150)
        print(f"Saved: neg_ctrl_{tool.lower()}_params_fp")

        p7 = (
            ggplot(tool_df, aes(x="settings", y="n_species", fill="sample_label"))
            + geom_col(position="dodge")
            + geom_text(
                aes(label="n_species"),
                position=position_dodge(width=0.9),
                angle=90, va="bottom", size=5,
                format_string="{:.0f}",
            )
            + labs(
                x="Parameter settings",
                y="Number of species called",
                fill="Sample",
                title=f"{tool}: Effect of parameters on species calls",
            )
            + theme_minimal()
            + theme(
                axis_text_x=element_text(rotation=45, hjust=1, size=7),
                figure_size=(12, 6),
            )
        )
        p7.save(OUTPUT_DIR / f"neg_ctrl_{tool.lower()}_params_species.pdf", dpi=150)
        p7.save(OUTPUT_DIR / f"neg_ctrl_{tool.lower()}_params_species.png", dpi=150)
        print(f"Saved: neg_ctrl_{tool.lower()}_params_species")

    print(f"\nDone! All outputs in: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
