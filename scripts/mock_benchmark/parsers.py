"""
parsers.py

Parse each tool's native output into a unified, per-taxon long table projected
onto the EsViritu lineage backbone.

Unified row schema (one row per reported taxon):
    tool, sample, param_set,
    kingdom, phylum, tclass, order, family, genus, species, subspecies,
    count, count_type

`count` is a per-taxon (exactly-at-node) quantity, NOT clade-cumulative:
  - esviritu   : read_count            (leaf-level; lineage from its own columns)
  - sylph      : Sequence_abundance     (leaf-level; lineage via Accession)
  - kraken2    : taxon_reads  (col 3)   (lineage via TaxID)
  - centrifuger: numReads               (lineage via TaxID)
  - metabuli   : taxon_count  (col 3)   (lineage via TaxID)

Rolling up to a coarser rank is done downstream by grouping on that rank's
lineage value, which is correct for per-taxon counts.
"""

import os
import sys

import pandas as pd

from taxonomy import ESV_RANKS, TaxonomyMaps, empty_lineage


# ── filename -> (sample, param_set) helpers ───────────────────────────────────

def strip_sample_param(stem: str, sample: str) -> str:
    """Extract the parameter token from a file/dir stem like 'mock_comb1_mnk10'."""
    if stem.startswith(sample + "_"):
        return stem[len(sample) + 1:]
    if stem == sample:
        return "default"
    return stem


def _row(tool, sample, param, lineage, count, count_type):
    r = {"tool": tool, "sample": sample, "param_set": param}
    r.update({rank: lineage.get(rank, "") for rank in ESV_RANKS})
    r["count"] = float(count)
    r["count_type"] = count_type
    return r


# ── per-tool parsers ──────────────────────────────────────────────────────────

def parse_esviritu(path: str, sample: str, param: str, maps: TaxonomyMaps,
                   unmapped: dict) -> list:
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    rows = []
    for rec in df.itertuples(index=False):
        d = rec._asdict()
        lineage = {rank: d.get(rank, "") for rank in ESV_RANKS}
        try:
            cnt = float(d.get("read_count", 0) or 0)
        except ValueError:
            cnt = 0.0
        if cnt <= 0:
            continue
        rows.append(_row("esviritu", sample, param, lineage, cnt, "read_count"))
    return rows


def parse_sylph(path: str, sample: str, param: str, maps: TaxonomyMaps,
                unmapped: dict) -> list:
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    rows = []
    for rec in df.itertuples(index=False):
        d = rec._asdict()
        acc = d.get("Contig_name", "")
        lineage = maps.lineage_for_accession(acc)
        if lineage is None:
            unmapped["sylph"] = unmapped.get("sylph", 0) + 1
            continue
        try:
            cnt = float(d.get("Sequence_abundance", 0) or 0)
        except ValueError:
            cnt = 0.0
        if cnt <= 0:
            continue
        rows.append(_row("sylph", sample, param, lineage, cnt, "seq_abundance"))
    return rows


def _parse_taxid_report(path, sample, param, maps, tool, unmapped,
                        taxid_col, count_col, sep="\t", header=0,
                        names=None, comment=None):
    df = pd.read_csv(path, sep=sep, dtype=str, header=header,
                     names=names, comment=comment).fillna("")
    rows = []
    for rec in df.itertuples(index=False):
        d = rec._asdict()
        taxid = d.get(taxid_col, "")
        try:
            cnt = float(d.get(count_col, 0) or 0)
        except ValueError:
            cnt = 0.0
        if cnt <= 0:
            continue
        lineage = maps.lineage_for_taxid(taxid)
        if lineage is None:
            unmapped[tool] = unmapped.get(tool, 0) + 1
            continue
        rows.append(_row(tool, sample, param, lineage, cnt, "reads"))
    return rows


def parse_kraken2(path: str, sample: str, param: str, maps: TaxonomyMaps,
                  unmapped: dict) -> list:
    # No header. Columns: pct, clade_reads, taxon_reads, rank_code, taxid, name
    names = ["pct", "clade_reads", "taxon_reads", "rank_code", "taxid", "name"]
    return _parse_taxid_report(
        path, sample, param, maps, "kraken2", unmapped,
        taxid_col="taxid", count_col="taxon_reads",
        header=None, names=names,
    )


def parse_centrifuger(path: str, sample: str, param: str, maps: TaxonomyMaps,
                      unmapped: dict) -> list:
    # Header: name, taxID, taxRank, genomeSize, numReads, numUniqueReads, abundance
    return _parse_taxid_report(
        path, sample, param, maps, "centrifuger", unmapped,
        taxid_col="taxID", count_col="numReads", header=0,
    )


def parse_metabuli(path: str, sample: str, param: str, maps: TaxonomyMaps,
                   unmapped: dict) -> list:
    # Header line starts with '#':
    # clade_proportion, clade_count, taxon_count, rank, taxID, name
    names = ["clade_proportion", "clade_count", "taxon_count",
             "rank", "taxID", "name"]
    return _parse_taxid_report(
        path, sample, param, maps, "metabuli", unmapped,
        taxid_col="taxID", count_col="taxon_count",
        header=None, names=names, comment="#",
    )


PARSERS = {
    "esviritu": parse_esviritu,
    "sylph": parse_sylph,
    "kraken2": parse_kraken2,
    "centrifuger": parse_centrifuger,
    "metabuli": parse_metabuli,
}
