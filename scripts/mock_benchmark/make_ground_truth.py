#!/usr/bin/env python3
"""
make_ground_truth.py

Generate a ground-truth table for mock virome communities built with mut_virome
and sequenced with `iss generate`.

Read counts are scraped directly from FASTQ headers using the ISS naming
convention:  @{seqid}_{read_index}_0/{pair}

Inputs
------
--selected_genomes   mut_virome selected_genomes.tsv  (Assembly, ANI, taxonomy, targets)
--mutation_summary   mut_virome mutation_summary.tsv  (realized mutations, FASTA paths)
                     [default: <selected_genomes dir>/mutation_summary.tsv]
--esv_metadata       EsViritu virus_pathogen_database.all_metadata.tsv
                     (used as fallback seqid→Assembly map if FASTA files are absent)
--fastq_r1           ISS-generated R1 FASTQ  (required; counts here represent read PAIRS)
--fastq_r2           ISS-generated R2 FASTQ  (optional; when given read_count = R1+R2 reads)
--output             Output TSV path (default: stdout)

Output
------
One row per viral Assembly with columns from selected_genomes.tsv, realized
mutation counts from mutation_summary.tsv, and read-level statistics:
  read_count    – reads attributable to this Assembly (pairs if R1-only, else R1+R2)
  total_reads   – total reads in the provided FASTQ(s)
  read_fraction – read_count / total_reads
A stderr summary reports total / viral / background read counts.
"""

import argparse
import re
import sys
import os
from collections import defaultdict

import pandas as pd


# ── ISS header parsing ────────────────────────────────────────────────────────
# ISS read names: {seqid}_{read_index}_{cpu_id}/{pair_number}
_PAIR_RE = re.compile(r"/\d+$")
_ISS_SUFFIX_RE = re.compile(r"_\d+_\d+$")


def extract_seqid(read_name: str) -> str:
    """Strip ISS trailing fields to recover the source sequence ID."""
    name = _PAIR_RE.sub("", read_name)
    return _ISS_SUFFIX_RE.sub("", name)


# ── FASTA helpers ─────────────────────────────────────────────────────────────

def fasta_seqids(fasta_path: str) -> list:
    """Return all sequence IDs (text before first space on '>' lines)."""
    seqids = []
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                seqids.append(line[1:].split()[0])
    return seqids


# ── FASTQ read counter ────────────────────────────────────────────────────────

def count_reads_fastq(fastq_path: str) -> dict:
    """Stream a FASTQ and return {seqid: read_count}."""
    counts: dict = defaultdict(int)
    with open(fastq_path) as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                read_name = line[1:].rstrip()   # drop '@' and newline
                counts[extract_seqid(read_name)] += 1
    return counts


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    ap = argparse.ArgumentParser(
        description="Build ground-truth table for a mut_virome / iss generate mock community",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument("--selected_genomes", required=True,
                    help="mut_virome selected_genomes.tsv")
    ap.add_argument("--mutation_summary", default=None,
                    help="mut_virome mutation_summary.tsv "
                         "(default: <selected_genomes dir>/mutation_summary.tsv)")
    ap.add_argument("--esv_metadata", required=True,
                    help="EsViritu virus_pathogen_database.all_metadata.tsv")
    ap.add_argument("--fastq_r1", required=True,
                    help="ISS-generated R1 FASTQ")
    ap.add_argument("--fastq_r2", default=None,
                    help="ISS-generated R2 FASTQ (optional)")
    ap.add_argument("--output", default="-",
                    help="Output TSV path (default: stdout)")
    return ap.parse_args()


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    mutvir_dir = os.path.dirname(os.path.abspath(args.selected_genomes))
    mut_summary_path = (
        args.mutation_summary
        or os.path.join(mutvir_dir, "mutation_summary.tsv")
    )

    # ── Load input tables ─────────────────────────────────────────────────────
    sel = pd.read_csv(args.selected_genomes, sep="\t", dtype=str)
    mut = pd.read_csv(mut_summary_path, sep="\t", dtype=str)
    esv = pd.read_csv(
        args.esv_metadata, sep="\t",
        usecols=["Accession", "Assembly"],
        dtype=str,
    )

    # ── Build seqid → Assembly map ────────────────────────────────────────────
    # Primary source: FASTA files listed in mutation_summary
    seqid_to_asm: dict = {}
    missing_fastas: list = []

    for _, row in mut.iterrows():
        fasta_rel = row["genome_fasta"]
        fasta_abs = (
            fasta_rel
            if os.path.isabs(fasta_rel)
            else os.path.join(mutvir_dir, fasta_rel)
        )
        asm = row["Assembly"]
        if os.path.isfile(fasta_abs):
            for sid in fasta_seqids(fasta_abs):
                seqid_to_asm[sid] = asm
        else:
            missing_fastas.append((asm, fasta_abs))

    if missing_fastas:
        print(
            f"WARNING: {len(missing_fastas)} genome FASTA(s) not found; "
            "using EsViritu metadata as fallback for those assemblies.",
            file=sys.stderr,
        )
        for asm, path in missing_fastas:
            print(f"  missing: {path}  (Assembly: {asm})", file=sys.stderr)

    # Fallback: EsViritu metadata (Accession → Assembly)
    viral_asm_set = set(sel["Assembly"].values)
    for acc, asm in zip(esv["Accession"], esv["Assembly"]):
        if acc not in seqid_to_asm and asm in viral_asm_set:
            seqid_to_asm[acc] = asm

    n_asm_mapped = len(set(seqid_to_asm.values()))
    print(
        f"Mapped {len(seqid_to_asm):,} sequence IDs → {n_asm_mapped} assemblies",
        file=sys.stderr,
    )

    # ── Count reads ───────────────────────────────────────────────────────────
    print(f"Counting reads in {args.fastq_r1} ...", file=sys.stderr)
    seqid_counts = count_reads_fastq(args.fastq_r1)
    total_reads = sum(seqid_counts.values())

    if args.fastq_r2:
        print(f"Counting reads in {args.fastq_r2} ...", file=sys.stderr)
        r2_counts = count_reads_fastq(args.fastq_r2)
        total_reads += sum(r2_counts.values())
        for sid, cnt in r2_counts.items():
            seqid_counts[sid] += cnt

    # Aggregate to Assembly
    asm_reads: dict = defaultdict(int)
    background_reads = 0
    for sid, cnt in seqid_counts.items():
        asm = seqid_to_asm.get(sid)
        if asm:
            asm_reads[asm] += cnt
        else:
            background_reads += cnt

    viral_reads = sum(asm_reads.values())
    print(f"Total reads : {total_reads:,}", file=sys.stderr)
    print(f"Viral reads : {viral_reads:,}", file=sys.stderr)
    print(f"Background  : {background_reads:,}", file=sys.stderr)

    # ── Assemble output table ─────────────────────────────────────────────────
    mut_keep = [
        "Assembly",
        "realized_snp_count",
        "realized_indel_count",
        "genome_fasta",
        "notes",
    ]
    out = sel.merge(mut[mut_keep], on="Assembly", how="left")

    out["read_count"] = (
        out["Assembly"].map(asm_reads).fillna(0).astype(int)
    )
    out["total_reads"] = total_reads
    out["read_fraction"] = out["read_count"] / total_reads

    out = out.sort_values("read_count", ascending=False).reset_index(drop=True)

    # ── Write ─────────────────────────────────────────────────────────────────
    fh = open(args.output, "w") if args.output != "-" else sys.stdout
    out.to_csv(fh, sep="\t", index=False)
    if args.output != "-":
        fh.close()

    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
