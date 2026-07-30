#!/bin/bash
#SBATCH --job-name=sylph_bench
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=logs/sylph_%j.out
#SBATCH --error=logs/sylph_%j.err

set -euo pipefail

# ============================================================
# Sylph parameter sweep
# Sketch reads once at -c 100, then profile with varying --min-number-kmers
# ============================================================

SAMPLE_SHEET="${1:-scripts/sample_sheet.tsv}"
DB="/data/tisza/analyses/mjt_projects/esviritu_benchmarks/sylph_tests/esv_324_slyph_c100_DB.syldb"
OUTDIR="results/sylph"
THREADS=16
SKETCH_C=100
LOGFILE="${OUTDIR}/sylph_runs.log"
TIME_FMT='\n[time]\t%E real,\t%U user,\t%S sys,\t%K amem,\t%M mmem'

source /data/tisza/analyses/mjt_sandbox/conda_mjt.init && conda activate sylph

# Parameter sweep values
MIN_NUMBER_KMERS_VALUES=(10 20 50)

mkdir -p "${OUTDIR}" logs

echo "========================================" >> "${LOGFILE}"
echo "Sylph benchmark run: $(date)" >> "${LOGFILE}"
echo "Sample sheet: ${SAMPLE_SHEET}" >> "${LOGFILE}"
echo "Database: ${DB}" >> "${LOGFILE}"
echo "Sketch -c: ${SKETCH_C}" >> "${LOGFILE}"
echo "Min number kmers values: ${MIN_NUMBER_KMERS_VALUES[*]}" >> "${LOGFILE}"
echo "========================================" >> "${LOGFILE}"

while IFS=$'\t' read -r SAMPLE_ID READ1 READ2; do
    # Skip header/comment lines
    [[ "${SAMPLE_ID}" =~ ^#.*$ ]] && continue
    [[ -z "${SAMPLE_ID}" ]] && continue

    SAMPLE_OUTDIR="${OUTDIR}/${SAMPLE_ID}"
    mkdir -p "${SAMPLE_OUTDIR}"

    # --- Step 1: Sketch reads (once per sample) ---
    SKETCH_FILE="${SAMPLE_OUTDIR}/${SAMPLE_ID}.sylsp"

    if [[ -n "${READ2}" ]]; then
        EXPECTED_SKETCH="${SAMPLE_OUTDIR}/$(basename ${READ1}).paired.sylsp"
    else
        EXPECTED_SKETCH="${SAMPLE_OUTDIR}/$(basename ${READ1}).sylsp"
    fi

    if [[ ! -f "${EXPECTED_SKETCH}" ]]; then
        echo "SKETCH: ${SAMPLE_ID} -c ${SKETCH_C}" | tee -a "${LOGFILE}"
        SKETCH_TIME_LOG="${SAMPLE_OUTDIR}/${SAMPLE_ID}.sketch.time.log"

        if [[ -n "${READ2}" ]]; then
            /usr/bin/time -f "${TIME_FMT}" -o "${SKETCH_TIME_LOG}" \
                sylph sketch -1 "${READ1}" -2 "${READ2}" -c "${SKETCH_C}" -t "${THREADS}" -d "${SAMPLE_OUTDIR}"
        else
            /usr/bin/time -f "${TIME_FMT}" -o "${SKETCH_TIME_LOG}" \
                sylph sketch -r "${READ1}" -c "${SKETCH_C}" -t "${THREADS}" -d "${SAMPLE_OUTDIR}"
        fi

        echo "  SKETCH TIME: $(cat ${SKETCH_TIME_LOG})" | tee -a "${LOGFILE}"
        echo "  SKETCH DONE: ${SAMPLE_ID}" | tee -a "${LOGFILE}"
    else
        echo "SKIP SKETCH: ${SAMPLE_ID} — sketch exists" | tee -a "${LOGFILE}"
    fi

    # --- Step 2: Profile with varying --min-number-kmers ---
    for MNK in "${MIN_NUMBER_KMERS_VALUES[@]}"; do
        PROFILE_OUT="${SAMPLE_OUTDIR}/${SAMPLE_ID}_mnk${MNK}.profile.tsv"

        # Skip if output already exists
        if [[ -f "${PROFILE_OUT}" ]]; then
            echo "SKIP: ${SAMPLE_ID} mnk=${MNK} — output exists" | tee -a "${LOGFILE}"
            continue
        fi

        echo "RUN: ${SAMPLE_ID} mnk=${MNK}" | tee -a "${LOGFILE}"

        CMD="sylph profile ${DB} ${EXPECTED_SKETCH} --min-number-kmers ${MNK} -t ${THREADS} -o ${PROFILE_OUT}"
        echo "  CMD: ${CMD}" >> "${LOGFILE}"

        TIME_LOG="${SAMPLE_OUTDIR}/${SAMPLE_ID}_mnk${MNK}.time.log"
        /usr/bin/time -f "${TIME_FMT}" -o "${TIME_LOG}" bash -c "${CMD}"
        echo "  TIME: $(cat ${TIME_LOG})" | tee -a "${LOGFILE}"
        echo "  DONE: ${SAMPLE_ID} mnk=${MNK}" | tee -a "${LOGFILE}"
    done

done < "${SAMPLE_SHEET}"

echo "All Sylph runs complete: $(date)" | tee -a "${LOGFILE}"
