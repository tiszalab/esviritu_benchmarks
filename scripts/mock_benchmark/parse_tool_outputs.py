#!/usr/bin/env python3
"""
parse_tool_outputs.py

Discover every tool/parameter result file under a results directory and parse
them into a single unified, per-taxon long table projected onto the EsViritu
lineage backbone.

Usage
-----
python scripts/benchmark/parse_tool_outputs.py \
    --results-dir benchmarks_mock/results \
    --metadata    esviritu_DBs/v3.2.4b/virus_pathogen_database.all_metadata.tsv \
    --nodes-dmp   esviritu_DBs/v3.2.4/kraken2db/taxonomy/nodes.dmp \
    --out         benchmarks_mock/analysis/unified_long.tsv

Samples are auto-discovered from the per-tool subdirectories. The output is a
TSV with columns:
    tool, sample, param_set, <8 EsViritu ranks>, count, count_type
"""

import argparse
import glob
import os
import sys

import pandas as pd

# make sibling modules importable when run as a script
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from taxonomy import TaxonomyMaps, ESV_RANKS          # noqa: E402
import parsers as P                                    # noqa: E402


# (tool, subdir-glob for files, filename-suffix to strip, parser)
# metabuli is special-cased (results live in per-param sub-directories).
FILE_TOOLS = {
    "esviritu":    ("*.tax_profile.tsv", ".tax_profile.tsv", P.parse_esviritu),
    "sylph":       ("*.profile.tsv",     ".profile.tsv",     P.parse_sylph),
    "centrifuger": ("*.report.tsv",      ".report.tsv",      P.parse_centrifuger),
    "kraken2":     ("*.report.txt",      ".report.txt",      P.parse_kraken2),
}


def discover_samples(results_dir: str) -> list:
    samples = set()
    for tool in list(FILE_TOOLS) + ["metabuli"]:
        tdir = os.path.join(results_dir, tool)
        if not os.path.isdir(tdir):
            continue
        for entry in os.listdir(tdir):
            if os.path.isdir(os.path.join(tdir, entry)):
                samples.add(entry)
    return sorted(samples)


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--nodes-dmp", default=None,
                    help="NCBI nodes.dmp for rolling up intermediate TaxIDs")
    ap.add_argument("--out", required=True)
    return ap.parse_args()


def main():
    args = parse_args()
    print("Loading taxonomy maps ...", file=sys.stderr)
    maps = TaxonomyMaps(args.metadata, nodes_dmp=args.nodes_dmp)
    print(f"  accessions={len(maps.accession_to_lineage):,} "
          f"assemblies={len(maps.assembly_to_lineage):,} "
          f"taxids={len(maps.taxid_to_lineage):,} "
          f"nodes.dmp={len(maps.parent):,}", file=sys.stderr)

    samples = discover_samples(args.results_dir)
    print(f"Samples: {samples}", file=sys.stderr)

    all_rows = []
    unmapped = {}
    n_files = 0

    for sample in samples:
        # standard single-file tools
        for tool, (fglob, suffix, parser) in FILE_TOOLS.items():
            sdir = os.path.join(args.results_dir, tool, sample)
            for path in sorted(glob.glob(os.path.join(sdir, fglob))):
                stem = os.path.basename(path)[: -len(suffix)]
                param = P.strip_sample_param(stem, sample)
                rows = parser(path, sample, param, maps, unmapped)
                all_rows.extend(rows)
                n_files += 1

        # metabuli: results/metabuli/<sample>/<param_dir>/<param_dir>_report.tsv
        mdir = os.path.join(args.results_dir, "metabuli", sample)
        for path in sorted(glob.glob(os.path.join(mdir, "*", "*_report.tsv"))):
            param_dir = os.path.basename(os.path.dirname(path))
            param = P.strip_sample_param(param_dir, sample)
            rows = P.parse_metabuli(path, sample, param, maps, unmapped)
            all_rows.extend(rows)
            n_files += 1

    if not all_rows:
        print("ERROR: no rows parsed; check --results-dir layout.", file=sys.stderr)
        sys.exit(1)

    out_df = pd.DataFrame(all_rows,
                          columns=["tool", "sample", "param_set"] +
                                  ESV_RANKS + ["count", "count_type"])

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    out_df.to_csv(args.out, sep="\t", index=False)

    print(f"Parsed {n_files} files -> {len(out_df):,} rows", file=sys.stderr)
    print(f"Wrote {args.out}", file=sys.stderr)
    if unmapped:
        print(f"Unmapped taxa (dropped) per tool: {unmapped}", file=sys.stderr)
    # quick coverage summary
    summary = (out_df.groupby(["tool", "sample"])["param_set"]
               .nunique().reset_index(name="param_sets"))
    print(summary.to_string(index=False), file=sys.stderr)


if __name__ == "__main__":
    main()
