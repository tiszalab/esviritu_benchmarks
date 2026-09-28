#!/usr/bin/env python3
"""Run BLASTN and summarize best non-overlapping, non-N alignment columns."""

import argparse
import csv
import shutil
import subprocess
import sys
from pathlib import Path


BLAST_FIELDS = (
    "qseqid",
    "sseqid",
    "qlen",
    "slen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qseq",
    "sseq",
)

SUMMARY_FIELDS = (
    "query_id",
    "subject_id",
    "query_length",
    "subject_length",
    "raw_hsp_count",
    "hsp_count",
    "raw_alignment_columns",
    "excluded_n_columns",
    "excluded_overlap_columns",
    "non_n_alignment_columns",
    "identities",
    "mismatches",
    "gap_columns",
    "percent_identity",
    "query_aligned_non_n_bases",
    "subject_aligned_non_n_bases",
    "query_coverage_percent",
    "subject_coverage_percent",
    "total_bit_score",
    "best_bit_score",
    "best_evalue",
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query", required=True, help="Ground-truth genome FASTA")
    parser.add_argument("--subject", required=True, help="Reconstructed-sequence FASTA")
    parser.add_argument("--blast-output", required=True, help="Raw BLASTN TSV output path")
    parser.add_argument("--summary-output", required=True, help="Query/subject summary TSV output path")
    return parser.parse_args()


def run_blastn(query, subject, output):
    blastn = shutil.which("blastn")
    if blastn is None:
        raise RuntimeError("blastn was not found on PATH")
    for label, path in (("query", query), ("subject", subject)):
        if not Path(path).is_file():
            raise FileNotFoundError(f"{label} FASTA does not exist: {path}")
    command = [
        blastn,
        "-query",
        query,
        "-subject",
        subject,
        "-out",
        output,
        "-outfmt",
        "6 " + " ".join(BLAST_FIELDS),
    ]
    subprocess.run(command, check=True)


def aligned_positions(start, end, aligned_sequence):
    step = 1 if end >= start else -1
    coordinate = start
    for base in aligned_sequence:
        if base == "-":
            yield None
        else:
            yield coordinate
            coordinate += step


def new_summary(row):
    return {
        "query_id": row["qseqid"],
        "subject_id": row["sseqid"],
        "query_length": int(row["qlen"]),
        "subject_length": int(row["slen"]),
        "raw_hsp_count": 0,
        "hsp_count": 0,
        "raw_alignment_columns": 0,
        "excluded_n_columns": 0,
        "excluded_overlap_columns": 0,
        "non_n_alignment_columns": 0,
        "identities": 0,
        "mismatches": 0,
        "gap_columns": 0,
        "query_positions": set(),
        "subject_positions": set(),
        "total_bit_score": 0.0,
        "best_bit_score": float("-inf"),
        "best_evalue": float("inf"),
    }


def add_hsp(summary, row, line_number):
    qseq = row["qseq"].upper()
    sseq = row["sseq"].upper()
    if len(qseq) != len(sseq):
        raise ValueError(f"line {line_number}: qseq and sseq have different aligned lengths")

    qpositions = aligned_positions(int(row["qstart"]), int(row["qend"]), qseq)
    spositions = aligned_positions(int(row["sstart"]), int(row["send"]), sseq)
    summary["raw_hsp_count"] += 1
    summary["raw_alignment_columns"] += len(qseq)
    contributed = False

    for qbase, sbase, qposition, sposition in zip(qseq, sseq, qpositions, spositions):
        if qbase == "N" or sbase == "N":
            summary["excluded_n_columns"] += 1
            continue
        query_overlaps = qposition is not None and qposition in summary["query_positions"]
        subject_overlaps = sposition is not None and sposition in summary["subject_positions"]
        if query_overlaps or subject_overlaps:
            summary["excluded_overlap_columns"] += 1
            continue
        contributed = True
        summary["non_n_alignment_columns"] += 1
        if qbase == "-" or sbase == "-":
            summary["gap_columns"] += 1
        elif qbase == sbase:
            summary["identities"] += 1
        else:
            summary["mismatches"] += 1
        if qposition is not None:
            summary["query_positions"].add(qposition)
        if sposition is not None:
            summary["subject_positions"].add(sposition)

    if contributed:
        bit_score = float(row["bitscore"])
        evalue = float(row["evalue"])
        summary["hsp_count"] += 1
        summary["total_bit_score"] += bit_score
        summary["best_bit_score"] = max(summary["best_bit_score"], bit_score)
        summary["best_evalue"] = min(summary["best_evalue"], evalue)


def parse_blast(blast_output):
    grouped_hsps = {}
    with open(blast_output, newline="") as handle:
        reader = csv.DictReader(handle, fieldnames=BLAST_FIELDS, delimiter="\t")
        for line_number, row in enumerate(reader, 1):
            if None in row.values():
                raise ValueError(f"line {line_number}: expected {len(BLAST_FIELDS)} tab-separated fields")
            key = (row["qseqid"], row["sseqid"])
            grouped_hsps.setdefault(key, []).append((line_number, row))

    summaries = {}
    for key, hsps in grouped_hsps.items():
        summary = new_summary(hsps[0][1])
        for line_number, row in hsps:
            if (
                summary["query_length"] != int(row["qlen"])
                or summary["subject_length"] != int(row["slen"])
            ):
                raise ValueError(f"line {line_number}: inconsistent sequence length for {key}")
        ranked_hsps = sorted(
            hsps,
            key=lambda item: (
                -float(item[1]["bitscore"]),
                float(item[1]["evalue"]),
                -len(item[1]["qseq"]),
            ),
        )
        for line_number, row in ranked_hsps:
            add_hsp(summary, row, line_number)
        summaries[key] = summary
    return summaries


def format_number(value):
    return f"{value:.10g}"


def write_summaries(summaries, summary_output):
    with open(summary_output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for key in sorted(summaries):
            summary = summaries[key]
            non_n_columns = summary["non_n_alignment_columns"]
            query_length = summary["query_length"]
            subject_length = summary["subject_length"]
            row = {field: summary[field] for field in SUMMARY_FIELDS if field in summary}
            row.update(
                {
                    "percent_identity": format_number(
                        100.0 * summary["identities"] / non_n_columns if non_n_columns else 0.0
                    ),
                    "query_aligned_non_n_bases": len(summary["query_positions"]),
                    "subject_aligned_non_n_bases": len(summary["subject_positions"]),
                    "query_coverage_percent": format_number(
                        100.0 * len(summary["query_positions"]) / query_length if query_length else 0.0
                    ),
                    "subject_coverage_percent": format_number(
                        100.0 * len(summary["subject_positions"]) / subject_length if subject_length else 0.0
                    ),
                    "total_bit_score": format_number(summary["total_bit_score"]),
                    "best_bit_score": format_number(summary["best_bit_score"]),
                    "best_evalue": format_number(summary["best_evalue"]),
                }
            )
            writer.writerow(row)


def main():
    args = parse_args()
    try:
        run_blastn(args.query, args.subject, args.blast_output)
        summaries = parse_blast(args.blast_output)
        write_summaries(summaries, args.summary_output)
    except (FileNotFoundError, OSError, RuntimeError, subprocess.CalledProcessError, ValueError) as error:
        sys.exit(f"ERROR: {error}")
    sys.stderr.write(
        f"Wrote {len(summaries)} query/subject summaries to {args.summary_output}\n"
        f"Raw BLASTN alignments: {args.blast_output}\n"
    )


if __name__ == "__main__":
    main()
