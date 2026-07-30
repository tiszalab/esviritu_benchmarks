"""
taxonomy.py

Build the EsViritu-style lineage backbone and the mappings needed to project
every tool's output into it.

The EsViritu metadata (virus_pathogen_database.all_metadata.tsv) is the Rosetta
Stone: it links Accession <-> Assembly <-> EsViritu lineage (k__..t__) <-> NCBI
TaxID.  The `LineageTaxIDs` column lists the NCBI TaxID for each EsViritu rank,
aligned position-by-position with the lineage name columns, which lets us map an
arbitrary NCBI TaxID (reported by kraken2/centrifuger/metabuli) to its EsViritu
lineage value at its own rank and all ancestor ranks.

An optional NCBI nodes.dmp parent map lets us walk intermediate TaxIDs (not
themselves EsViritu ranks) up to the nearest mapped ancestor so no classified
read is silently dropped.
"""

import os
import sys
from typing import Optional

import pandas as pd

# EsViritu rank columns (order matters: kingdom -> subspecies/strain)
ESV_RANKS = [
    "kingdom", "phylum", "tclass", "order",
    "family", "genus", "species", "subspecies",
]
# Single-letter prefixes used in both lineage values ("f__Kolmioviridae")
# and LineageTaxIDs tokens ("f__2842321").
RANK_LETTERS = ["k", "p", "c", "o", "f", "g", "s", "t"]


def _parse_lineage_taxids(field: str) -> list:
    """Return a list of 8 TaxID strings (or '') aligned with ESV_RANKS.

    Input looks like: 'k__2732396;p__2497569;c__;o__;f__1980418;...;t__'.
    """
    out = ["" for _ in ESV_RANKS]
    if not isinstance(field, str) or not field:
        return out
    for tok in field.split(";"):
        tok = tok.strip()
        if "__" not in tok:
            continue
        prefix, _, taxid = tok.partition("__")
        if prefix in RANK_LETTERS:
            out[RANK_LETTERS.index(prefix)] = taxid.strip()
    return out


class TaxonomyMaps:
    """Container for all lookups derived from the EsViritu metadata."""

    def __init__(self, metadata_path: str, nodes_dmp: Optional[str] = None):
        self.metadata_path = metadata_path
        self.accession_to_lineage: dict = {}
        self.assembly_to_lineage: dict = {}
        # NCBI TaxID -> {rank: value} for the taxid's own rank and all ancestors
        self.taxid_to_lineage: dict = {}
        # NCBI child TaxID -> parent TaxID (from nodes.dmp, optional)
        self.parent: dict = {}

        self._build_from_metadata()
        if nodes_dmp:
            self._load_nodes_dmp(nodes_dmp)

    # ── construction ──────────────────────────────────────────────────────────
    def _build_from_metadata(self) -> None:
        usecols = ESV_RANKS + ["Accession", "Assembly", "LineageTaxIDs"]
        df = pd.read_csv(
            self.metadata_path, sep="\t", usecols=usecols, dtype=str
        ).fillna("")

        for row in df.itertuples(index=False):
            d = row._asdict()
            vals = [d[r] for r in ESV_RANKS]
            lineage = {ESV_RANKS[i]: vals[i] for i in range(len(ESV_RANKS))}

            acc = d["Accession"]
            asm = d["Assembly"]
            if acc:
                self.accession_to_lineage.setdefault(acc, lineage)
            if asm:
                self.assembly_to_lineage.setdefault(asm, lineage)

            tids = _parse_lineage_taxids(d["LineageTaxIDs"])
            for i, tid in enumerate(tids):
                if not tid:
                    continue
                # ancestors + self: ranks 0..i
                anc = {ESV_RANKS[j]: vals[j] for j in range(i + 1)}
                # Keep the most specific (largest i) definition seen.
                prev = self.taxid_to_lineage.get(tid)
                if prev is None or len(anc) >= len(prev):
                    self.taxid_to_lineage[tid] = anc

    def _load_nodes_dmp(self, nodes_dmp: str) -> None:
        if not os.path.isfile(nodes_dmp):
            print(f"WARNING: nodes.dmp not found at {nodes_dmp}; "
                  "intermediate TaxIDs will not be rolled up.", file=sys.stderr)
            return
        with open(nodes_dmp) as fh:
            for line in fh:
                parts = [p.strip() for p in line.split("|")]
                if len(parts) >= 2:
                    self.parent[parts[0]] = parts[1]

    # ── lookups ─────────────────────────────────────────────────────────────--
    def lineage_for_taxid(self, taxid) -> Optional[dict]:
        """Resolve an NCBI TaxID to an EsViritu lineage dict.

        Direct hit if the TaxID is itself an EsViritu rank; otherwise walk up the
        NCBI tree (if nodes.dmp loaded) to the nearest mapped ancestor.
        """
        if taxid is None:
            return None
        tid = str(taxid).strip()
        if not tid or tid in ("0",):
            return None
        hit = self.taxid_to_lineage.get(tid)
        if hit is not None:
            return hit
        # walk ancestors
        seen = 0
        cur = self.parent.get(tid)
        while cur and seen < 100:
            hit = self.taxid_to_lineage.get(cur)
            if hit is not None:
                return hit
            nxt = self.parent.get(cur)
            if nxt == cur:
                break
            cur = nxt
            seen += 1
        return None

    def lineage_for_accession(self, accession: str) -> Optional[dict]:
        if accession is None:
            return None
        acc = str(accession).strip()
        if acc in self.accession_to_lineage:
            return self.accession_to_lineage[acc]
        # tolerate version-less accessions on either side
        base = acc.split(".")[0]
        for cand in (acc, base):
            for key in (cand, cand + ".1"):
                if key in self.accession_to_lineage:
                    return self.accession_to_lineage[key]
        return None


def empty_lineage() -> dict:
    return {r: "" for r in ESV_RANKS}
