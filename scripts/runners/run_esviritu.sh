#!/bin/bash
#SBATCH --job-name=esviritu_bench
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=logs/esviritu_%j.out
#SBATCH --error=logs/esviritu_%j.err

set -euo pipefail
set +u
# ============================================================
# EsViritu classification — no parameter sweep
# ============================================================

SAMPLE_SHEET="${1:-scripts/sample_sheet.tsv}"
DB="/data/tisza/analyses/mjt_projects/esviritu_benchmarks/esviritu_DBs/v3.2.4b"
OUTDIR="results/esviritu"
THREADS=16
LOGFILE="${OUTDIR}/esviritu_runs.log"
TIME_FMT='\n[time]\t%E real,\t%U user,\t%S sys,\t%K amem,\t%M mmem'
TIME_CMD="/usr/bin/time -f '${TIME_FMT}'"

# Activate conda environment
source /data/tisza/analyses/mjt_sandbox/conda_mjt.init && conda activate esv_spec

#parameter values
MODES=("sr" "sense")

mkdir -p "${OUTDIR}" logs

echo "========================================" >> "${LOGFILE}"
echo "EsViritu benchmark run: $(date)" >> "${LOGFILE}"
echo "Sample sheet: ${SAMPLE_SHEET}" >> "${LOGFILE}"
echo "Database: ${DB}" >> "${LOGFILE}"
echo "Mode values: ${MODES[*]}" >> "${LOGFILE}"
echo "========================================" >> "${LOGFILE}"

while IFS=$'\t' read -r SAMPLE_ID READ1 READ2; do
    # Skip header/comment lines
    [[ "${SAMPLE_ID}" =~ ^#.*$ ]] && continue
    [[ -z "${SAMPLE_ID}" ]] && continue

    SAMPLE_OUTDIR="${OUTDIR}/${SAMPLE_ID}"

    for MODE in "${MODES[@]}" ; do

        SAMPMODE="${SAMPLE_ID}_${MODE}"

        # Skip if output already exists
        if [[ -d "${SAMPLE_OUTDIR}" ]] && [[ -f "${SAMPLE_OUTDIR}/${SAMPMODE}_esviritu.log" ]]; then
            echo "SKIP: ${SAMPMODE} — output exists" | tee -a "${LOGFILE}"
            continue
        fi

        mkdir -p "${SAMPLE_OUTDIR}"

        echo "RUN: ${SAMPMODE}" | tee -a "${LOGFILE}"
        TIME_LOG="${SAMPLE_OUTDIR}/${SAMPMODE}.time.log"

        if [[ -n "${READ2}" ]]; then
            CMD="EsViritu \
                --db ${DB} \
                --sample ${SAMPMODE} \
                -mmP ${MODE} \
                --output_dir ${SAMPLE_OUTDIR} \
                --cpu ${THREADS} \
                --reads ${READ1} ${READ2}"
        else
            if [[ "${MODE}" == "sr" ]] ; then MODE="lr:hq" ; fi
            CMD="EsViritu \
                --db ${DB} \
                --sample ${SAMPMODE} \
                --output_dir ${SAMPLE_OUTDIR} \
                --cpu ${THREADS} \
                --read_format single \
                -mmP ${MODE} \
                --reads ${READ1}"
        fi

        echo "  CMD: ${CMD}" >> "${LOGFILE}"
        /usr/bin/time -f "${TIME_FMT}" -o "${TIME_LOG}" bash -c "${CMD}"
        echo "  TIME: $(cat ${TIME_LOG})" | tee -a "${LOGFILE}"
        echo "  DONE: ${SAMPMODE}" | tee -a "${LOGFILE}"
    done

done < "${SAMPLE_SHEET}"

echo "All EsViritu runs complete: $(date)" | tee -a "${LOGFILE}"
