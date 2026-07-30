#!/bin/bash
#SBATCH --job-name=kraken2_bench
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=logs/kraken2_%j.out
#SBATCH --error=logs/kraken2_%j.err

set -euo pipefail

# ============================================================
# Kraken2 (k2 CLI) parameter sweep
# Parameters: --confidence × --minimum-hit-groups
# ============================================================

SAMPLE_SHEET="${1:-scripts/sample_sheet.tsv}"
DB="/data/tisza/analyses/mjt_projects/esviritu_benchmarks/esviritu_DBs/v3.2.4/kraken2db"
OUTDIR="results/kraken2"
THREADS=16
LOGFILE="${OUTDIR}/kraken2_runs.log"
TIME_FMT='\n[time]\t%E real,\t%U user,\t%S sys,\t%K amem,\t%M mmem'

source /data/tisza/analyses/mjt_sandbox/conda_mjt.init && conda activate kraken

# Parameter sweep values
CONFIDENCE_VALUES=(0 0.2 0.5)
MIN_HIT_GROUPS_VALUES=(2 3 5)

mkdir -p "${OUTDIR}" logs

echo "========================================" >> "${LOGFILE}"
echo "Kraken2 benchmark run: $(date)" >> "${LOGFILE}"
echo "Sample sheet: ${SAMPLE_SHEET}" >> "${LOGFILE}"
echo "Database: ${DB}" >> "${LOGFILE}"
echo "Confidence values: ${CONFIDENCE_VALUES[*]}" >> "${LOGFILE}"
echo "Min hit groups values: ${MIN_HIT_GROUPS_VALUES[*]}" >> "${LOGFILE}"
echo "========================================" >> "${LOGFILE}"

while IFS=$'\t' read -r SAMPLE_ID READ1 READ2; do
    # Skip header/comment lines
    [[ "${SAMPLE_ID}" =~ ^#.*$ ]] && continue
    [[ -z "${SAMPLE_ID}" ]] && continue

    SAMPLE_OUTDIR="${OUTDIR}/${SAMPLE_ID}"
    mkdir -p "${SAMPLE_OUTDIR}"

    for CONF in "${CONFIDENCE_VALUES[@]}"; do
        for MHG in "${MIN_HIT_GROUPS_VALUES[@]}"; do
            OUTPREFIX="${SAMPLE_OUTDIR}/${SAMPLE_ID}_conf${CONF}_mhg${MHG}"
            REPORT="${OUTPREFIX}.report.txt"
            OUTPUT="${OUTPREFIX}.out"

            # Skip if output already exists
            if [[ -f "${REPORT}" ]]; then
                echo "SKIP: ${SAMPLE_ID} conf=${CONF} mhg=${MHG} — output exists" | tee -a "${LOGFILE}"
                continue
            fi

            echo "RUN: ${SAMPLE_ID} conf=${CONF} mhg=${MHG}" | tee -a "${LOGFILE}"

            CMD="k2 classify --db ${DB} --threads ${THREADS} --confidence ${CONF} --minimum-hit-groups ${MHG} --report ${REPORT} --output ${OUTPUT}"

            if [[ -n "${READ2}" ]]; then
                CMD="${CMD} --paired ${READ1} ${READ2}"
            else
                CMD="${CMD} ${READ1}"
            fi

            echo "  CMD: ${CMD}" >> "${LOGFILE}"

            TIME_LOG="${OUTPREFIX}.time.log"
            /usr/bin/time -f "${TIME_FMT}" -o "${TIME_LOG}" bash -c "${CMD}"
            echo "  TIME: $(cat ${TIME_LOG})" | tee -a "${LOGFILE}"
            echo "  DONE: ${SAMPLE_ID} conf=${CONF} mhg=${MHG}" | tee -a "${LOGFILE}"
        done
    done

done < "${SAMPLE_SHEET}"

echo "All Kraken2 runs complete: $(date)" | tee -a "${LOGFILE}"
