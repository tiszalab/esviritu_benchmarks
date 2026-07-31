"""
Compute strict and ANI-capped tool-versus-truth metrics in the EsViritu lineage
backbone. Capped scoring requires at most one truth genome per family because
predicted counts are assigned to the unique matching truth family before being
rolled to that genome's permitted rank.
"""

import pandas as pd

from taxonomy import ESV_RANKS

EVAL_RANKS = ["family", "genus", "species", "subspecies"]
DEFAULT_ANI_TIERS = {
    100.0: "subspecies", 98.0: "subspecies",
    95.0: "species", 90.0: "species",
    85.0: "genus", 80.0: "genus",
}
READ_UNIT_MULTIPLIER = {
    "esviritu": 1.0,
    "kraken2": 2.0,
    "centrifuger": 2.0,
    "metabuli": 2.0,
    "sylph": None,
}
_RANK_DEPTH = {rank: index for index, rank in enumerate(ESV_RANKS)}


def is_unclassified(value) -> bool:
    if value is None:
        return True
    value = str(value).strip()
    return not value or "unclassified" in value.lower() or (len(value) == 3 and value.endswith("__"))


def normalize(vector: pd.Series) -> pd.Series:
    total = vector.sum()
    return vector * 0.0 if total <= 0 else vector / total


def bray_curtis(truth: pd.Series, prediction: pd.Series) -> float:
    index = truth.index.union(prediction.index)
    if len(index) == 0:
        return float("nan")
    truth = truth.reindex(index).fillna(0.0)
    prediction = prediction.reindex(index).fillna(0.0)
    denominator = (truth + prediction).sum()
    return float("nan") if denominator <= 0 else float((truth - prediction).abs().sum() / denominator)


def abundance_stats(truth: pd.Series, prediction: pd.Series) -> dict:
    index = truth.index.union(prediction.index)
    truth = truth.reindex(index).fillna(0.0)
    prediction = prediction.reindex(index).fillna(0.0)
    result = {"l1": float((truth - prediction).abs().sum()), "l2": float(((truth - prediction) ** 2).sum() ** 0.5), "pearson": float("nan"), "spearman": float("nan")}
    if len(index) >= 2 and truth.std() > 0 and prediction.std() > 0:
        result["pearson"] = float(truth.corr(prediction, method="pearson"))
        result["spearman"] = float(truth.corr(prediction, method="spearman"))
    return result


def detection_at_rank(truth_values: set, prediction_values: set) -> dict:
    tp = len(truth_values & prediction_values)
    fn = len(truth_values - prediction_values)
    fp = len(prediction_values - truth_values)
    precision = tp / (tp + fp) if tp + fp else float("nan")
    recall = tp / (tp + fn) if tp + fn else float("nan")
    f1 = 2 * precision * recall / (precision + recall) if precision and recall and precision + recall else float("nan")
    return {"TP": tp, "FP": fp, "FN": fn, "precision": precision, "recall": recall, "F1": f1}


def classified_value(lineage: dict, rank: str) -> tuple[str, str]:
    for index in range(_RANK_DEPTH[rank], -1, -1):
        candidate_rank = ESV_RANKS[index]
        value = lineage.get(candidate_rank, "")
        if not is_unclassified(value):
            return candidate_rank, value
    return rank, ""


def rank_vector(table: pd.DataFrame, rank: str, count_column: str) -> pd.Series:
    if table.empty:
        return pd.Series(dtype=float)
    subset = table[[rank, count_column]].copy()
    subset[count_column] = pd.to_numeric(subset[count_column], errors="coerce").fillna(0.0)
    subset = subset[~subset[rank].map(is_unclassified)]
    if subset.empty:
        return pd.Series(dtype=float)
    return subset.groupby(rank)[count_column].sum()


def metric_row(truth: pd.Series, prediction: pd.Series, rank: str, truth_mode: str, multiplier, rpm_factor: float) -> dict:
    detection = detection_at_rank(set(truth[truth > 0].index), set(prediction[prediction > 0].index))
    result = {
        "rank": rank,
        "truth_mode": truth_mode,
        "truth_total_rpm": float(truth.sum() * rpm_factor),
        "n_truth_taxa": int((truth > 0).sum()),
        "n_pred_taxa": int((prediction > 0).sum()),
        **detection,
    }
    if multiplier is None:
        return {**result, "bray_curtis_rpm": float("nan"), "bray_curtis_rel": float("nan"), "pred_total_rpm": float("nan"), "l1": float("nan"), "l2": float("nan"), "pearson": float("nan"), "spearman": float("nan")}
    truth_rpm = truth * rpm_factor
    prediction_rpm = prediction * multiplier * rpm_factor
    return {
        **result,
        "bray_curtis_rpm": bray_curtis(truth_rpm, prediction_rpm),
        "bray_curtis_rel": bray_curtis(normalize(truth), normalize(prediction)),
        "pred_total_rpm": float(prediction_rpm.sum()),
        **abundance_stats(truth_rpm, prediction_rpm),
    }


def capped_vectors(truth_df: pd.DataFrame, pred_df: pd.DataFrame, ani_tiers: dict) -> tuple[pd.Series, pd.Series, list[dict]]:
    truth_by_family = {}
    truth_records = []
    for record in truth_df.itertuples(index=False):
        lineage = record._asdict()
        try:
            ani = float(lineage.get("ani_level"))
        except (TypeError, ValueError):
            ani = 100.0
        cap_rank = ani_tiers.get(ani)
        if cap_rank is None:
            raise ValueError(f"ANI level {ani} is not defined in ANI tiers")
        family = lineage.get("family", "")
        if is_unclassified(family):
            raise ValueError("ANI-capped scoring requires a classified truth family")
        if family in truth_by_family:
            raise ValueError(f"ANI-capped scoring requires one truth genome per family; repeated family: {family}")
        used_rank, value = classified_value(lineage, cap_rank)
        truth_by_family[family] = (used_rank, value)
        truth_records.append({"lineage": lineage, "ani": ani, "expected_rank": used_rank, "expected_value": value})

    truth_counts = {}
    for item in truth_records:
        key = f"{item['expected_rank']}|{item['expected_value']}"
        truth_counts[key] = truth_counts.get(key, 0.0) + float(item["lineage"].get("read_count", 0))

    pred_counts = {}
    for record in pred_df.itertuples(index=False):
        lineage = record._asdict()
        family = lineage.get("family", "")
        if family in truth_by_family:
            cap_rank, value = classified_value(lineage, truth_by_family[family][0])
        else:
            cap_rank, value = classified_value(lineage, "subspecies")
        if is_unclassified(value):
            continue
        key = f"{cap_rank}|{value}"
        pred_counts[key] = pred_counts.get(key, 0.0) + float(lineage.get("count", 0.0))

    return pd.Series(truth_counts, dtype=float), pd.Series(pred_counts, dtype=float), truth_records


def compute_community_and_detection(truth_df: pd.DataFrame, pred_df: pd.DataFrame, tool: str, total_reads: float, ani_tiers: dict = None) -> tuple:
    ani_tiers = ani_tiers or DEFAULT_ANI_TIERS
    multiplier = READ_UNIT_MULTIPLIER.get(tool, 1.0)
    rpm_factor = 1e6 / total_reads if total_reads and total_reads > 0 else float("nan")
    rows = []
    for rank in EVAL_RANKS:
        rows.append(metric_row(rank_vector(truth_df, rank, "read_count"), rank_vector(pred_df, rank, "count"), rank, "strict", multiplier, rpm_factor))

    try:
        capped_truth, capped_prediction, truth_records = capped_vectors(truth_df, pred_df, ani_tiers)
    except ValueError as error:
        if "one truth genome per family" not in str(error):
            raise
        return rows, None, []
    rows.append(metric_row(capped_truth, capped_prediction, "ani_capped", "capped", multiplier, rpm_factor))

    expected = set(capped_truth[capped_truth > 0].index)
    observed = set(capped_prediction[capped_prediction > 0].index)
    capped_detection = detection_at_rank(expected, observed)
    family_counts = {}
    for item in truth_records:
        key = (item["expected_rank"], item["expected_value"])
        family_counts[key] = family_counts.get(key, 0) + 1
    per_genome = []
    for item in truth_records:
        key = f"{item['expected_rank']}|{item['expected_value']}"
        raw_count = capped_prediction.get(key, 0.0)
        per_genome.append({
            "Assembly": item["lineage"].get("Assembly", ""),
            "ani_level": item["ani"],
            "expected_rank": item["expected_rank"],
            "expected_value": item["expected_value"],
            "read_count": item["lineage"].get("read_count", 0),
            "detected": key in observed,
            "assigned_reads": float("nan") if multiplier is None else raw_count * multiplier,
            "n_truth_sharing_clade": family_counts[(item["expected_rank"], item["expected_value"])],
        })
    return rows, capped_detection, per_genome
