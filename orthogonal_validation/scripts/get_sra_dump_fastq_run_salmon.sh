#!/bin/bash
# Adapted from BHARATH_get_sra_dump_fastq_run_salmon.sh
# Change: single-end mode (GSE171975 samples are MARS-Seq, layout = SINGLE)
#   - replaced parallel-fastq-dump with fasterq-dump (built into sra-tools)
#   - replaced -1/-2 with -r in salmon quant
#
# USAGE: bash scripts/get_sra_dump_fastq_run_salmon.sh <SRP_accession> [path/to/salmon/index]
# Run from: orthogonal_validation/
#
# Example: bash scripts/get_sra_dump_fastq_run_salmon.sh SRP314631

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${SCRIPT_DIR}/.."

GENTROME="${BASE_DIR}/data/gentrome.fa"
DECOY="${BASE_DIR}/data/decoys.txt"
SALMON_INDEX="${BASE_DIR}/data/transcriptome_ext_index"
OUTDIR="${BASE_DIR}/results/${1}_salmon_$(date +"%Y-%m-%d")"
LOG="${OUTDIR}/${1}.log"
CORES=7

# If no index path given, build it from gentrome
if [ -z $2 ]; then
    salmon index -t $GENTROME \
        -i $SALMON_INDEX \
        --decoys $DECOY \
        -k 31 \
        --keepDuplicates
else
    SALMON_INDEX=$2
fi

mkdir -p $OUTDIR
echo "Run on $(date)" > $LOG
salmon --version >> $LOG
fasterq-dump --version >> $LOG

cat "${BASE_DIR}/meta/${1}_Acc_List.txt" | while read line
do
    # Fetch the .sra file
    prefetch -O "${BASE_DIR}/raw_data/" $line \
        >> $LOG

    # Extract the fastq file (single-end: fasterq-dump, produces one file)
    if [ ! -f "${BASE_DIR}/temp/${line}.fastq.gz" ]; then
        fasterq-dump --threads 6 \
            --outdir "${BASE_DIR}/temp/" \
            "${BASE_DIR}/raw_data/${line}/${line}.sra" \
            >> $LOG 2>&1
        gzip "${BASE_DIR}/temp/${line}.fastq"
    fi

    # Run salmon quantification (single-end: -r instead of -1/-2)
    salmon quant -i $SALMON_INDEX \
        -l A \
        -p $CORES \
        --seqBias --gcBias -q \
        -r "${BASE_DIR}/temp/${line}.fastq.gz" \
        -o "${OUTDIR}/${line}/" \
        >> $LOG
done
