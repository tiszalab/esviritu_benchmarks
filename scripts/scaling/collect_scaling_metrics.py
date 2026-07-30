#!/usr/bin/env python3
"""Parse .time.log files from the scaling benchmark into a single TSV.

Reads results_scaling/<tool>/<tool>_r<READCOUNT>_c<CPU>[.<step>].time.log
where each log holds a /usr/bin/time line of the form:

    [time]  0:03.36 real,  21.31 user,  0.14 sys,  0 amem, 144552 mmem

Emits one row per step, plus a synthesized "total" row per run
(sum of real/user/sys times, max of memory) for multi-step tools.
"""
import argparse
import re
import sys
from pathlib import Path

# results_scaling/<tool>/<tool>_r<RC>_c<CPU>[.<step>].time.log
NAME_RE = re.compile(
    r"^(?P<tool>[a-z0-9]+)_r(?P<rc>\d+)_c(?P<cpu>\d+)(?:\.(?P<step>[a-z]+))?\.time\.log$"
)
TIME_RE = re.compile(
    r"([\d:.]+)\s+real,\s+([\d.]+)\s+user,\s+([\d.]+)\s+sys,\s+(\d+)\s+amem,\s+(\d+)\s+mmem"
)


def parse_elapsed(text: str) -> float:
    """Convert %E ([h:]m:ss[.ss]) to seconds."""
    parts = text.split(":")
    parts = [float(p) for p in parts]
    if len(parts) == 3:
        h, m, s = parts
    elif len(parts) == 2:
        h, m, s = 0.0, parts[0], parts[1]
    else:
        h, m, s = 0.0, 0.0, parts[0]
    return h * 3600 + m * 60 + s


def parse_log(path: Path):
    m = TIME_RE.search(path.read_text())
    if not m:
        return None
    real, user, sys_, amem, mmem = m.groups()
    return {
        "real_s": round(parse_elapsed(real), 2),
        "real_raw": real,
        "user_s": float(user),
        "sys_s": float(sys_),
        "max_rss_kb": int(mmem),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", default="results_scaling",
                    help="directory containing per-tool subdirs (default: results_scaling)")
    ap.add_argument("-o", "--output", default="results_scaling/scaling_metrics.tsv",
                    help="output TSV path")
    args = ap.parse_args()

    root = Path(args.root)
    logs = sorted(root.glob("*/*.time.log"))
    if not logs:
        print(f"No .time.log files found under {root}", file=sys.stderr)
        return 1

    rows = []            # per-step rows
    runs: dict = {}      # (tool, rc, cpu) -> list of step metrics for totals

    for path in logs:
        nm = NAME_RE.match(path.name)
        if not nm:
            print(f"WARN: unrecognized name {path.name}", file=sys.stderr)
            continue
        tool = nm["tool"]
        rc = int(nm["rc"])
        cpu = int(nm["cpu"])
        step = nm["step"] or "run"

        metrics = parse_log(path)
        if metrics is None:
            print(f"WARN: no timing line in {path}", file=sys.stderr)
            continue

        rows.append({"tool": tool, "read_count": rc, "cpus": cpu,
                     "step": step, **metrics})
        runs.setdefault((tool, rc, cpu), []).append(metrics)

    # Synthesize total rows for runs with more than one step
    for (tool, rc, cpu), steps in runs.items():
        if len(steps) <= 1:
            continue
        rows.append({
            "tool": tool, "read_count": rc, "cpus": cpu, "step": "total",
            "real_s": round(sum(s["real_s"] for s in steps), 2),
            "real_raw": "",
            "user_s": round(sum(s["user_s"] for s in steps), 2),
            "sys_s": round(sum(s["sys_s"] for s in steps), 2),
            "max_rss_kb": max(s["max_rss_kb"] for s in steps),
        })

    cols = ["tool", "read_count", "cpus", "step",
            "real_s", "real_raw", "user_s", "sys_s", "max_rss_kb", "max_rss_mb"]

    def sort_key(r):
        step_order = {"run": 0, "classify": 0, "sketch": 0,
                      "quant": 1, "profile": 1, "total": 2}
        return (r["tool"], r["read_count"], r["cpus"], step_order.get(r["step"], 0))

    rows.sort(key=sort_key)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            r["max_rss_mb"] = round(r["max_rss_kb"] / 1024, 1)
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    print(f"Wrote {len(rows)} rows to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
