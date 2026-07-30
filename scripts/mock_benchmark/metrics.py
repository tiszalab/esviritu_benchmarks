"""
metrics.py

Compute tool-vs-truth benchmark metrics in the EsViritu lineage backbone:

  * Bray-Curtis dissimilarity (primary) at family / genus / species, on
    within-viral relative abundances.
  * ANI-aware detection: each ground-truth genome is scored at the deepest rank
    it can plausibly be recovered at, given its mutation ANI tier.
  * Per-rank detection (precision / recall / F1) at family / genus / species.
  * Abundance accuracy: Pearson / Spearman / L1 / L2 on relative abundances.

Subspecies/strain (t__) is folded into species by default (EVAL_RANKS).
"""

import pandas as pd

from taxonomy import ESV_RANKS

# Ranks at which community vectors / detection are evaluated.
EVAL_RANKS = ["family", "genus", "species"]

# ANI tier -> deepest rank considered recoverable.
DEFAULT_ANI_TIERS = {
    100.0: "species", 98.0: "species",
    95.0: "genus",   90.0: "genus",
    85.0: "family",  80.0: "family",
}

# Read-unit reconciliation: multiplier to convert each tool's native per-taxon
# count into INDIVIDUAL reads (the unit of the ground-truth read_count = R1+R2).
#   - esviritu   counts individual reads               -> x1
#   - kraken2/centrifuger/metabuli run --paired, so a
#     classified read PAIR counts as 1                 -> x2
#   - sylph reports a coverage-based Sequence_abundance
#     (NOT a read count) and is EXCLUDED from Bray-Curtis (None)
READ_UNIT_MULTIPLIER = {
    "esviritu": 1.0,
    "kraken2": 2.0,
    "centrifuger": 2.0,
    "metabuli": 2.0,
    "sylph": None,
}

# rank -> how coarse it is (0 = kingdom ... 7 = subspecies)
_RANK_DEPTH = {r: i for i, r in enumerate(ESV_RANKS)}


# ── helpers ───────────────────────────────────────────────────────────────────

def is_unclassified(val) -> bool:
    """True for empty, prefix-only ('', 'g__'), or 'unclassified' lineage values."""
    if val is None:
        return True
    v = str(val).strip()
    if not v:
        return True
    if "unclassified" in v.lower():
        return True
    return len(v) == 3 and v.endswith("__")


def rank_vector(long_df: pd.DataFrame, rank: str) -> pd.Series:
    """Sum per-taxon `count` grouped by the lineage value at `rank`.

    Rows that are unclassified at `rank` are dropped (they were classified only
    to a coarser rank and cannot contribute to this rank's community vector).
    """
    if long_df.empty:
        return pd.Series(dtype=float)
    sub = long_df[[rank, "count"]].copy()
    mask = ~sub[rank].map(is_unclassified)
    sub = sub[mask]
    if sub.empty:
        return pd.Series(dtype=float)
    return sub.groupby(rank)["count"].sum()


def normalize(vec: pd.Series) -> pd.Series:
    total = vec.sum()
    if total <= 0:
        return vec * 0.0
    return vec / total


def bray_curtis(truth: pd.Series, pred: pd.Series) -> float:
    """Bray-Curtis dissimilarity over the union of taxa (missing = 0)."""
    idx = truth.index.union(pred.index)
    if len(idx) == 0:
        return float("nan")
    a = truth.reindex(idx).fillna(0.0)
    b = pred.reindex(idx).fillna(0.0)
    denom = (a + b).sum()
    if denom <= 0:
        return float("nan")
    return float((a - b).abs().sum() / denom)


def abundance_stats(truth: pd.Series, pred: pd.Series) -> dict:
    idx = truth.index.union(pred.index)
    a = truth.reindex(idx).fillna(0.0)
    b = pred.reindex(idx).fillna(0.0)
    out = {
        "l1": float((a - b).abs().sum()),
        "l2": float(((a - b) ** 2).sum() ** 0.5),
        "pearson": float("nan"),
        "spearman": float("nan"),
    }
    if len(idx) >= 2 and a.std() > 0 and b.std() > 0:
        out["pearson"] = float(a.corr(b, method="pearson"))
        out["spearman"] = float(a.corr(b, method="spearman"))
    return out


def detection_at_rank(truth_vals: set, pred_vals: set) -> dict:
    tp = len(truth_vals & pred_vals)
    fn = len(truth_vals - pred_vals)
    fp = len(pred_vals - truth_vals)
    precision = tp / (tp + fp) if (tp + fp) else float("nan")
    recall = tp / (tp + fn) if (tp + fn) else float("nan")
    f1 = (2 * precision * recall / (precision + recall)
          if precision and recall and (precision + recall) else float("nan"))
    return {"TP": tp, "FP": fp, "FN": fn,
            "precision": precision, "recall": recall, "F1": f1}


def _classified_value(lineage_row: dict, rank: str) -> str:
    """Value at `rank`, falling back to the deepest classified ancestor."""
    depth = _RANK_DEPTH[rank]
    for i in range(depth, -1, -1):
        r = ESV_RANKS[i]
        v = lineage_row.get(r, "")
        if not is_unclassified(v):
            return r, v
    return rank, ""


# ── per (tool, param, sample) computation ─────────────────────────────────────

def compute_community_and_detection(truth_df: pd.DataFrame,
                                    pred_df: pd.DataFrame,
                                    tool: str,
                                    total_reads: float,
                                    ani_tiers: dict = None) -> tuple:
    """Return (rank_metrics_rows, ani_rows, per_genome_rows).

    truth_df    : ground-truth table (Assembly, ESV_RANKS, read_count, ani_level)
    pred_df     : unified long rows for ONE (tool, param_set, sample)
    tool        : tool name (selects the read-unit multiplier)
    total_reads : per-sample total (individual) reads, the RPM denominator
    """
    ani_tiers = ani_tiers or DEFAULT_ANI_TIERS
    mult = READ_UNIT_MULTIPLIER.get(tool, 1.0)
    rpm_factor = (1e6 / total_reads) if (total_reads and total_reads > 0) \
        else float("nan")

    rank_rows = []
    pred_vals_by_rank = {}
    pred_counts_by_rank = {}
    for rank in EVAL_RANKS:
        t_counts = _truth_rank_vector(truth_df, rank)
        p_counts = rank_vector(pred_df, rank)
        pred_vals_by_rank[rank] = set(p_counts[p_counts > 0].index)
        pred_counts_by_rank[rank] = p_counts

        # detection is presence/absence -> unit-independent
        det = detection_at_rank(set(t_counts.index),
                                set(p_counts[p_counts > 0].index))

        truth_total_rpm = float(t_counts.sum() * rpm_factor)

        if mult is None:
            # sylph: coverage-based output is not a read count -> no Bray-Curtis
            bc_rpm = float("nan")
            bc_rel = float("nan")
            ab = {"l1": float("nan"), "l2": float("nan"),
                  "pearson": float("nan"), "spearman": float("nan")}
            pred_total_rpm = float("nan")
        else:
            # RPM-based Bray-Curtis (primary): fixed denominator preserves
            # magnitude, so sensitivity & over-classification are captured.
            t_rpm = t_counts * rpm_factor
            p_rpm = p_counts * (mult * rpm_factor)
            bc_rpm = bray_curtis(t_rpm, p_rpm)
            # relative/compositional Bray-Curtis (shape only)
            bc_rel = bray_curtis(normalize(t_counts), normalize(p_counts))
            ab = abundance_stats(t_rpm, p_rpm)
            pred_total_rpm = float(p_rpm.sum())

        rank_rows.append({
            "rank": rank,
            "bray_curtis_rpm": bc_rpm,
            "bray_curtis_rel": bc_rel,
            "truth_total_rpm": truth_total_rpm,
            "pred_total_rpm": pred_total_rpm,
            "n_truth_taxa": int((t_counts > 0).sum()),
            "n_pred_taxa": int((p_counts > 0).sum()),
            **det,
            **ab,
        })

    # ANI-aware detection (per ground-truth genome)
    per_genome = []
    tp = fn = 0
    expected_truth_by_rank = {r: set() for r in EVAL_RANKS}

    # how many truth genomes resolve to each (rank, lineage value): the tool's
    # assigned reads at a clade are POOLED, so they must be shared among all
    # truth genomes that map to that same value.
    genome_keys = []
    for rec in truth_df.itertuples(index=False):
        d = rec._asdict()
        try:
            ani = float(d.get("ani_level"))
        except (TypeError, ValueError):
            ani = 100.0
        used_rank, exp_val = _classified_value(d, ani_tiers.get(ani, "species"))
        genome_keys.append((d, ani, used_rank, exp_val))
    share_counts = {}
    for _, _, used_rank, exp_val in genome_keys:
        if not is_unclassified(exp_val):
            share_counts[(used_rank, exp_val)] = \
                share_counts.get((used_rank, exp_val), 0) + 1

    for d, ani, used_rank, exp_val in genome_keys:
        detected = (not is_unclassified(exp_val)) and \
                   (exp_val in pred_vals_by_rank.get(used_rank, set()))
        if not is_unclassified(exp_val):
            expected_truth_by_rank.setdefault(used_rank, set()).add(exp_val)
        if detected:
            tp += 1
        else:
            fn += 1

        # tool reads assigned to the clade this genome maps to, converted to
        # individual reads (the unit of ground-truth read_count). NaN for sylph
        # (coverage-based output, not a read count).
        n_sharing = share_counts.get((used_rank, exp_val), 0)
        if mult is None or is_unclassified(exp_val):
            assigned_reads = float("nan")
        else:
            raw = float(pred_counts_by_rank.get(used_rank,
                                                {}).get(exp_val, 0.0)) \
                if hasattr(pred_counts_by_rank.get(used_rank), "get") else 0.0
            assigned_reads = raw * mult

        per_genome.append({
            "Assembly": d.get("Assembly", ""),
            "ani_level": ani,
            "expected_rank": used_rank,
            "expected_value": exp_val,
            "read_count": d.get("read_count", 0),
            "detected": detected,
            "assigned_reads": assigned_reads,
            "n_truth_sharing_clade": n_sharing,
        })

    # ANI-aware false positives: predicted values not expected at each used rank
    fp_aware = 0
    for rank, exp_set in expected_truth_by_rank.items():
        if not exp_set:
            continue
        fp_aware += len(pred_vals_by_rank.get(rank, set()) - exp_set)

    precision_aware = tp / (tp + fp_aware) if (tp + fp_aware) else float("nan")
    recall_aware = tp / (tp + fn) if (tp + fn) else float("nan")
    f1_aware = (2 * precision_aware * recall_aware /
                (precision_aware + recall_aware)
                if precision_aware and recall_aware and
                (precision_aware + recall_aware) else float("nan"))

    ani_row = {
        "TP": tp, "FN": fn, "FP_aware": fp_aware,
        "precision_aware": precision_aware,
        "recall_aware": recall_aware,
        "F1_aware": f1_aware,
    }
    return rank_rows, ani_row, per_genome


def _truth_rank_vector(truth_df: pd.DataFrame, rank: str) -> pd.Series:
    sub = truth_df[[rank, "read_count"]].copy()
    sub["read_count"] = pd.to_numeric(sub["read_count"], errors="coerce").fillna(0.0)
    mask = ~sub[rank].map(is_unclassified)
    sub = sub[mask]
    if sub.empty:
        return pd.Series(dtype=float)
    return sub.groupby(rank)["read_count"].sum()
