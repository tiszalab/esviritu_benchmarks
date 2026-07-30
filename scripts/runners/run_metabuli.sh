#!/bin/bash
#SBATCH --job-name=metabuli_bench
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=logs/metabuli_%j.out
#SBATCH --error=logs/metabuli_%j.err

set -euo pipefail

# ============================================================
# Metabuli parameter sweep
# Parameters: --precise × --min-score × --priority-taxid
# ============================================================

SAMPLE_SHEET="${1:-scripts/sample_sheet.tsv}"
DB="/data/tisza/analyses/mjt_projects/esviritu_benchmarks/metabuli_tests/esv_metabuli_DB"
OUTDIR="results/metabuli"
THREADS=16
MAX_RAM=128
LOGFILE="${OUTDIR}/metabuli_runs.log"
TIME_FMT='\n[time]\t%E real,\t%U user,\t%S sys,\t%K amem,\t%M mmem'

source /data/tisza/analyses/mjt_sandbox/conda_mjt.init && conda activate metabuli

# Parameter sweep values
PRECISE_VALUES=(0 1)
MIN_SCORE_VALUES=(0 0.15)
PRIORITY_TAXID_VALUES=("none" "10239")

mkdir -p "${OUTDIR}" logs

echo "========================================" >> "${LOGFILE}"
echo "Metabuli benchmark run: $(date)" >> "${LOGFILE}"
echo "Sample sheet: ${SAMPLE_SHEET}" >> "${LOGFILE}"
echo "Database: ${DB}" >> "${LOGFILE}"
echo "Precise values: ${PRECISE_VALUES[*]}" >> "${LOGFILE}"
echo "Min score values: ${MIN_SCORE_VALUES[*]}" >> "${LOGFILE}"
echo "Priority taxid values: ${PRIORITY_TAXID_VALUES[*]}" >> "${LOGFILE}"
echo "========================================" >> "${LOGFILE}"

while IFS=$'\t' read -r SAMPLE_ID READ1 READ2; do
    # Skip header/comment lines
    [[ "${SAMPLE_ID}" =~ ^#.*$ ]] && continue
    [[ -z "${SAMPLE_ID}" ]] && continue

    SAMPLE_OUTDIR="${OUTDIR}/${SAMPLE_ID}"
    mkdir -p "${SAMPLE_OUTDIR}"

    for PREC in "${PRECISE_VALUES[@]}"; do
        for MS in "${MIN_SCORE_VALUES[@]}"; do
            for PTAX in "${PRIORITY_TAXID_VALUES[@]}"; do
                JOB_ID="${SAMPLE_ID}_prec${PREC}_ms${MS}_ptax${PTAX}"
                RUN_OUTDIR="${SAMPLE_OUTDIR}/${JOB_ID}"

                # Skip if output already exists
                if [[ -d "${RUN_OUTDIR}" ]] && [[ -f "${RUN_OUTDIR}/${JOB_ID}_report.tsv" ]]; then
                    echo "SKIP: ${JOB_ID} — output exists" | tee -a "${LOGFILE}"
                    continue
                fi

                mkdir -p "${RUN_OUTDIR}"

                echo "RUN: ${JOB_ID}" | tee -a "${LOGFILE}"

                # Build command
                CMD="metabuli classify"

                if [[ -n "${READ2}" ]]; then
                    CMD="${CMD} ${READ1} ${READ2}"
                else
                    CMD="${CMD} --seq-mode 1 ${READ1}"
                fi

                CMD="${CMD} ${DB} ${RUN_OUTDIR} ${JOB_ID}"
                CMD="${CMD} --threads ${THREADS} --max-ram ${MAX_RAM}"

                if [[ "${PREC}" -ne 0 ]]; then
                    CMD="${CMD} --precise ${PREC}"
                fi

                # Only pass --min-score if non-zero (avoid overriding precise preset with 0)
                if [[ "${MS}" != "0" ]]; then
                    CMD="${CMD} --min-score ${MS}"
                fi

                if [[ "${PTAX}" != "none" ]]; then
                    CMD="${CMD} --priority-taxid ${PTAX}"
                fi

                echo "  CMD: ${CMD}" >> "${LOGFILE}"

                TIME_LOG="${RUN_OUTDIR}/${JOB_ID}.time.log"
                /usr/bin/time -f "${TIME_FMT}" -o "${TIME_LOG}" bash -c "${CMD}"
                echo "  TIME: $(cat ${TIME_LOG})" | tee -a "${LOGFILE}"
                echo "  DONE: ${JOB_ID}" | tee -a "${LOGFILE}"
            done
        done
    done

done < "${SAMPLE_SHEET}"

echo "All Metabuli runs complete: $(date)" | tee -a "${LOGFILE}"
