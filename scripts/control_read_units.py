from typing import Optional


PAIRED_TOOLS = {"kraken2", "centrifuger", "metabuli"}


def base_tool(tool: str) -> str:
    return tool.lower().split("_", 1)[0].split(" ", 1)[0]


def positive_count_to_read_pairs(tool: str, count: float) -> Optional[float]:
    name = base_tool(tool)
    if name == "sylph":
        return None
    return float(count) / 2 if name == "esviritu" else float(count)


def negative_count_to_individual_reads(
    tool: str, count: float, read_type: str
) -> Optional[float]:
    name = base_tool(tool)
    if name == "sylph":
        return None
    multiplier = 2 if name in PAIRED_TOOLS and read_type == "illumina_short" else 1
    return float(count) * multiplier
