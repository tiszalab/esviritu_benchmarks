#!/bin/bash
#SBATCH --job-name=centrifuger_bench
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=logs/centrifuger_%j.out
#SBATCH --error=logs/centrifuger_%j.err

set -euo pipefail

# ============================================================
# Centrifuger parameter sweep
# Classification: --min-hitlen (auto, 50, 100)
# Quantification: --min-score (0, 100, 500)
# ============================================================

SAMPLE_SHEET="${1:-scripts/sample_sheet.tsv}"
DB="/data/tisza/analyses/mjt_projects/esviritu_benchmarks/centrifuger_tests/centrifuger_esv324"
OUTDIR="results/centrifuger"
THREADS=16
LOGFILE="${OUTDIR}/centrifuger_runs.log"
TIME_FMT='\n[time]\t%E real,\t%U user,\t%S sys,\t%K amem,\t%M mmem'

source /data/tisza/analyses/mjt_sandbox/conda_mjt.init && conda activate centrifuger

# Parameter sweep values
# "auto" means we don't pass --min-hitlen (use default)
MIN_HITLEN_VALUES=("auto" "50" "100")
MIN_SCORE_VALUES=(0 100 500)

mkdir -p "${OUTDIR}" logs

echo "========================================" >> "${LOGFILE}"
echo "Centrifuger benchmark run: $(date)" >> "${LOGFILE}"
echo "Sample sheet: ${SAMPLE_SHEET}" >> "${LOGFILE}"
echo "Database: ${DB}" >> "${LOGFILE}"
echo "Min hitlen values: ${MIN_HITLEN_VALUES[*]}" >> "${LOGFILE}"
echo "Min score values (quant): ${MIN_SCORE_VALUES[*]}" >> "${LOGFILE}"
echo "========================================" >> "${LOGFILE}"

while IFS=$'\t' read -r SAMPLE_ID READ1 READ2; do
    # Skip header/comment lines
    [[ "${SAMPLE_ID}" =~ ^#.*$ ]] && continue
    [[ -z "${SAMPLE_ID}" ]] && continue

    SAMPLE_OUTDIR="${OUTDIR}/${SAMPLE_ID}"
    mkdir -p "${SAMPLE_OUTDIR}"

    for MHL in "${MIN_HITLEN_VALUES[@]}"; do
        CLASSIFY_OUT="${SAMPLE_OUTDIR}/${SAMPLE_ID}_mhl${MHL}.classify.tsv"

        # --- Step 1: Classification ---
        if [[ ! -f "${CLASSIFY_OUT}" ]]; then
            echo "CLASSIFY: ${SAMPLE_ID} mhl=${MHL}" | tee -a "${LOGFILE}"

            HITLEN_FLAG=""
            if [[ "${MHL}" != "auto" ]]; then
                HITLEN_FLAG="--min-hitlen ${MHL}"
            fi

            if [[ -n "${READ2}" ]]; then
                CMD="centrifuger -x ${DB} -1 ${READ1} -2 ${READ2} -t ${THREADS} ${HITLEN_FLAG} > ${CLASSIFY_OUT}"
            else
                CMD="centrifuger -x ${DB} -1 ${READ1} -t ${THREADS} ${HITLEN_FLAG} > ${CLASSIFY_OUT}"
            fi

            echo "  CMD: ${CMD}" >> "${LOGFILE}"
            CLASSIFY_TIME_LOG="${SAMPLE_OUTDIR}/${SAMPLE_ID}_mhl${MHL}.classify.time.log"
            /usr/bin/time -f "${TIME_FMT}" -o "${CLASSIFY_TIME_LOG}" bash -c "${CMD}"
            echo "  CLASSIFY TIME: $(cat ${CLASSIFY_TIME_LOG})" | tee -a "${LOGFILE}"
            echo "  CLASSIFY DONE: ${SAMPLE_ID} mhl=${MHL}" | tee -a "${LOGFILE}"
        else
            echo "SKIP CLASSIFY: ${SAMPLE_ID} mhl=${MHL} — output exists" | tee -a "${LOGFILE}"
        fi

        # --- Step 2: Quantification with varying --min-score ---
        for MS in "${MIN_SCORE_VALUES[@]}"; do
            REPORT_OUT="${SAMPLE_OUTDIR}/${SAMPLE_ID}_mhl${MHL}_ms${MS}.report.tsv"

            if [[ -f "${REPORT_OUT}" ]]; then
                echo "SKIP QUANT: ${SAMPLE_ID} mhl=${MHL} ms=${MS} — output exists" | tee -a "${LOGFILE}"
                continue
            fi

            echo "QUANT: ${SAMPLE_ID} mhl=${MHL} ms=${MS}" | tee -a "${LOGFILE}"

            CMD="centrifuger-quant -x ${DB} -c ${CLASSIFY_OUT} --min-score ${MS} > ${REPORT_OUT}"
            echo "  CMD: ${CMD}" >> "${LOGFILE}"

            QUANT_TIME_LOG="${SAMPLE_OUTDIR}/${SAMPLE_ID}_mhl${MHL}_ms${MS}.quant.time.log"
            /usr/bin/time -f "${TIME_FMT}" -o "${QUANT_TIME_LOG}" bash -c "${CMD}"
            echo "  QUANT TIME: $(cat ${QUANT_TIME_LOG})" | tee -a "${LOGFILE}"
            echo "  QUANT DONE: ${SAMPLE_ID} mhl=${MHL} ms=${MS}" | tee -a "${LOGFILE}"
        done
    done

done < "${SAMPLE_SHEET}"

echo "All Centrifuger runs complete: $(date)" | tee -a "${LOGFILE}"
