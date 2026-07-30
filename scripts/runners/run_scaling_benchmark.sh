#!/bin/bash
set -euo pipefail

# ============================================================
# Scaling benchmark launcher
#
# Submits one sbatch job per (tool x cpu-count x read-count).
# Each tool runs at its "best" parameters (fixed) while CPUs and
# read count vary, to measure runtime/memory scaling.
#
# Tools & best params (see run_scaling_worker.sh):
#   EsViritu     : -mmP sr
#   centrifuger  : --min-hitlen 100 , quant --min-score 0
#   kraken2      : --confidence 0.2 --minimum-hit-groups 2
#   metabuli     : precise=0 , --min-score 0.15
#   sylph        : --min-number-kmers 20 (sketch -c 100)
#
# Usage:
#   scripts/run_scaling_benchmark.sh [TOOL ...]
#     (no args = all tools)
#   DRY_RUN=1 scripts/run_scaling_benchmark.sh   # print sbatch cmds only
# ============================================================

PROJECT_DIR="/data/tisza/analyses/mjt_projects/esviritu_benchmarks"
cd "${PROJECT_DIR}"

WORKER="${PROJECT_DIR}/scripts/run_scaling_worker.sh"
DRY_RUN="${DRY_RUN:-0}"

# --- Grid ---
CPUS=(2 4 8 16 32)
READ_COUNTS=(1000000 10000000 50000000 100000000 200000000)
ALL_TOOLS=(esviritu centrifuger kraken2 metabuli sylph)

# Tools to run: CLI args override the full list
if [[ $# -gt 0 ]]; then
    TOOLS=("$@")
else
    TOOLS=("${ALL_TOOLS[@]}")
fi

MEM="128G"
TIME_LIMIT="24:00:00"

mkdir -p "${PROJECT_DIR}/logs" "${PROJECT_DIR}/results_scaling"

echo "Project dir : ${PROJECT_DIR}"
echo "Tools       : ${TOOLS[*]}"
echo "CPU counts  : ${CPUS[*]}"
echo "Read counts : ${READ_COUNTS[*]}"
echo "Dry run     : ${DRY_RUN}"
echo

N=0
for TOOL in "${TOOLS[@]}"; do
    for RC in "${READ_COUNTS[@]}"; do
        for CPU in "${CPUS[@]}"; do
            JOBNAME="scale_${TOOL}_r${RC}_c${CPU}"
            SBATCH_CMD=(sbatch
                --job-name="${JOBNAME}"
                --cpus-per-task="${CPU}"
                --mem="${MEM}"
                --time="${TIME_LIMIT}"
                --output="logs/${JOBNAME}_%j.out"
                --error="logs/${JOBNAME}_%j.err"
                "${WORKER}" "${TOOL}" "${CPU}" "${RC}")

            if [[ "${DRY_RUN}" == "1" ]]; then
                echo "${SBATCH_CMD[*]}"
            else
                "${SBATCH_CMD[@]}"
            fi
            N=$((N+1))
        done
    done
done

echo
echo "Submitted ${N} job(s)."
