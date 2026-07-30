#!/bin/bash
#SBATCH --job-name=scale_bench
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=logs/scale_%j.out
#SBATCH --error=logs/scale_%j.err

set -euo pipefail
set +u

# ============================================================
# Scaling benchmark worker (one tool, one CPU count, one read count)
# Invoked by run_scaling_benchmark.sh via sbatch.
#
# Usage: run_scaling_worker.sh <TOOL> <CPU> <READCOUNT>
#   TOOL      : esviritu | centrifuger | kraken2 | metabuli | sylph
#   CPU       : number of threads / cpus-per-task
#   READCOUNT : read-count subset label (matches read file name)
# ============================================================

TOOL="${1:?TOOL required}"
CPU="${2:?CPU required}"
RC="${3:?READCOUNT required}"

# Resolve project root from this script's location (scripts/..)
PROJECT_DIR="/data/tisza/analyses/mjt_projects/esviritu_benchmarks"
cd "${PROJECT_DIR}"

# --- Databases ---
ESV_DB="${PROJECT_DIR}/esviritu_DBs/v3.2.4b"
CENT_DB="${PROJECT_DIR}/centrifuger_tests/centrifuger_esv324"
K2_DB="${PROJECT_DIR}/esviritu_DBs/v3.2.4/kraken2db"
MBULI_DB="${PROJECT_DIR}/metabuli_tests/esv_metabuli_DB"
SYLPH_DB="${PROJECT_DIR}/sylph_tests/esv_324_slyph_c100_DB.syldb"
MAX_RAM=128

# --- Reads (gzipped read-count subsets) ---
READ_BASE="${PROJECT_DIR}/viral_metagenome_reads/SRR33038190.fastp.${RC}"
READ1="${READ_BASE}.R1.fastq.gz"
READ2="${READ_BASE}.R2.fastq.gz"

if [[ ! -f "${READ1}" ]] || [[ ! -f "${READ2}" ]]; then
    echo "ERROR: read files not found for readcount ${RC}:" >&2
    echo "  ${READ1}" >&2
    echo "  ${READ2}" >&2
    exit 1
fi

OUTDIR="${PROJECT_DIR}/results_scaling/${TOOL}"
mkdir -p "${OUTDIR}" "${PROJECT_DIR}/logs"

TAG="${TOOL}_r${RC}_c${CPU}"
TIME_FMT='\n[time]\t%E real,\t%U user,\t%S sys,\t%K amem,\t%M mmem'
LOGFILE="${OUTDIR}/scaling_runs.log"

echo "========================================" >> "${LOGFILE}"
echo "Scaling run: ${TAG}  $(date)" >> "${LOGFILE}"
echo "TOOL=${TOOL} CPU=${CPU} READCOUNT=${RC}" >> "${LOGFILE}"
echo "Reads: ${READ1} ${READ2}" >> "${LOGFILE}"
echo "========================================" >> "${LOGFILE}"

source /data/tisza/analyses/mjt_sandbox/conda_mjt.init

run_timed () {
    # $1 = time-log path, remaining args = command string run via bash -c
    local tlog="$1"; shift
    echo "  CMD: $*" >> "${LOGFILE}"
    /usr/bin/time -f "${TIME_FMT}" -o "${tlog}" bash -c "$*"
    echo "  TIME: $(cat ${tlog})" | tee -a "${LOGFILE}"
}

case "${TOOL}" in

    esviritu)
        conda activate esv_spec
        SAMPLE_OUTDIR="${OUTDIR}/${TAG}"
        mkdir -p "${SAMPLE_OUTDIR}"
        TIME_LOG="${OUTDIR}/${TAG}.time.log"
        # best params: -mmP sr
        run_timed "${TIME_LOG}" \
            "EsViritu --db ${ESV_DB} --sample ${TAG} -mmP sr \
                --output_dir ${SAMPLE_OUTDIR} --cpu ${CPU} \
                --reads ${READ1} ${READ2}"
        ;;

    centrifuger)
        conda activate centrifuger
        # best params: --min-hitlen 100 (classify), --min-score 0 (quant)
        CLASSIFY_OUT="${OUTDIR}/${TAG}.classify.tsv"
        REPORT_OUT="${OUTDIR}/${TAG}.report.tsv"
        run_timed "${OUTDIR}/${TAG}.classify.time.log" \
            "centrifuger -x ${CENT_DB} -1 ${READ1} -2 ${READ2} -t ${CPU} --min-hitlen 100 > ${CLASSIFY_OUT}"
        run_timed "${OUTDIR}/${TAG}.quant.time.log" \
            "centrifuger-quant -x ${CENT_DB} -c ${CLASSIFY_OUT} --min-score 0 > ${REPORT_OUT}"
        ;;

    kraken2)
        conda activate kraken
        # best params: --confidence 0.2 --minimum-hit-groups 2
        REPORT="${OUTDIR}/${TAG}.report.txt"
        OUTPUT="${OUTDIR}/${TAG}.out"
        run_timed "${OUTDIR}/${TAG}.time.log" \
            "k2 classify --db ${K2_DB} --threads ${CPU} \
                --confidence 0.2 --minimum-hit-groups 2 \
                --report ${REPORT} --output ${OUTPUT} \
                --paired ${READ1} ${READ2}"
        ;;

    metabuli)
        conda activate metabuli
        # best params: precise=0 (no --precise flag), --min-score 0.15
        RUN_OUTDIR="${OUTDIR}/${TAG}"
        mkdir -p "${RUN_OUTDIR}"
        run_timed "${OUTDIR}/${TAG}.time.log" \
            "metabuli classify ${READ1} ${READ2} ${MBULI_DB} ${RUN_OUTDIR} ${TAG} \
                --threads ${CPU} --max-ram ${MAX_RAM} --min-score 0.15"
        ;;

    sylph)
        conda activate sylph
        # best params: --min-number-kmers 20 (sketch reads at -c 100)
        SKETCH_DIR="${OUTDIR}/${TAG}"
        mkdir -p "${SKETCH_DIR}"
        run_timed "${OUTDIR}/${TAG}.sketch.time.log" \
            "sylph sketch -1 ${READ1} -2 ${READ2} -c 100 -t ${CPU} -d ${SKETCH_DIR}"
        SKETCH_FILE="${SKETCH_DIR}/$(basename ${READ1}).paired.sylsp"
        PROFILE_OUT="${OUTDIR}/${TAG}.profile.tsv"
        run_timed "${OUTDIR}/${TAG}.profile.time.log" \
            "sylph profile ${SYLPH_DB} ${SKETCH_FILE} --min-number-kmers 20 -t ${CPU} -o ${PROFILE_OUT}"
        ;;

    *)
        echo "ERROR: unknown tool '${TOOL}'" >&2
        exit 1
        ;;
esac

echo "DONE: ${TAG}  $(date)" | tee -a "${LOGFILE}"
