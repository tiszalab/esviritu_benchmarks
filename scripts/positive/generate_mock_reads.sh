#!/bin/bash
# Generate mock reads with wgsim for benchmark genomes
# 150x150 nt paired-end reads
# 10 permutations of 1, 10, 100, 1000 read pairs per genome

set -euo pipefail

OUTDIR="/data/tisza/analyses/mjt_projects/esviritu_benchmarks/virus_genomes/mock_reads"
mkdir -p "$OUTDIR"

GENOME_DIRS=(
    "/data/tisza/analyses/mjt_projects/esviritu_benchmarks/virus_genomes/genomoviridae_for_test"
    "/data/tisza/analyses/mjt_projects/esviritu_benchmarks/virus_genomes/picornaviridae_for_test"
)

READ_COUNTS=(1 10 100 1000)
PERMUTATIONS=10
READ_LEN=150

for DIR in "${GENOME_DIRS[@]}"; do
    for GENOME in "$DIR"/*.fna; do
        # Extract genome name (accession) without extension
        GNAME=$(basename "$GENOME" .fna)
        
        for NREADS in "${READ_COUNTS[@]}"; do
            for PERM in $(seq 1 $PERMUTATIONS); do
                # Use a unique seed per permutation
                SEED=$((NREADS * 1000 + PERM))
                
                R1="${OUTDIR}/${GNAME}.${NREADS}.${PERM}.R1.fastq"
                R2="${OUTDIR}/${GNAME}.${NREADS}.${PERM}.R2.fastq"
                
                wgsim -1 ${READ_LEN} -2 ${READ_LEN} -N ${NREADS} -S ${SEED} \
                    "$GENOME" "$R1" "$R2" > /dev/null 2>&1
                
                echo "Done: ${GNAME} ${NREADS} reads, perm ${PERM}"
            done
        done
    done
done

echo "All mock reads generated in ${OUTDIR}"
