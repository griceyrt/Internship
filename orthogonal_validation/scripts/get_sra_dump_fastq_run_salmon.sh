#!/bin/bash
# Adapted from BHARATH_get_sra_dump_fastq_run_salmon.sh
# Change: single-end mode (GSE130613 samples are MARS-Seq, layout = SINGLE)
#   - replaced parallel-fastq-dump with fasterq-dump (built into sra-tools)
#   - replaced -1/-2 with -r in salmon quant
#   - -l A detects -l U (stranded-reverse, confirmed by strand_mapping_bias=0.0)
#   - 2026-06-18: added cutadapt adapter trimming step
#   - 2026-06-19: updated trimming per Bharath's Zoom — full trim:
#     adapter AGATCGGAAGAG + polyA (A{10}) + polyT (T{10}) + polyG (G{10})
#     polyG is a NextSeq 500 two-color chemistry artifact
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
    echo ""
    echo "======================================"
    echo "Processing sample: $line"
    echo "======================================"

    # Fetch the .sra file
    echo "[1/4] Prefetching $line..."
    prefetch -O "${BASE_DIR}/raw_data/" $line \
        >> $LOG

    # Extract the fastq file (single-end: fasterq-dump, produces one file)
    if [ ! -f "${BASE_DIR}/temp/${line}.fastq.gz" ]; then
        echo "[2/4] Converting .sra to .fastq.gz..."
        fasterq-dump --threads 6 \
            --outdir "${BASE_DIR}/temp/" \
            "${BASE_DIR}/raw_data/${line}/${line}.sra" \
            >> $LOG 2>&1
        gzip "${BASE_DIR}/temp/${line}.fastq"
    else
        echo "[2/4] FASTQ already exists, skipping fasterq-dump."
    fi

    # Trim adapters + polyA/T/G with cutadapt (updated 2026-06-19 per Bharath)
    TRIMMED="${BASE_DIR}/temp/${line}_trimmed_full.fastq.gz"
    if [ ! -f "$TRIMMED" ]; then
        echo "[3/4] Trimming adapters with cutadapt..."
        cutadapt \
            -a AGATCGGAAGAG \
            -a "A{12}" \
            -a "T{12}" \
            -a "G{12}" \
            --minimum-length 20 \
            -j $CORES \
            -o "$TRIMMED" \
            "${BASE_DIR}/temp/${line}.fastq.gz" \
            >> $LOG 2>&1
    else
        echo "[3/4] Trimmed file already exists, skipping cutadapt."
    fi

    # Run salmon quantification (single-end: -r instead of -1/-2)
    echo "[4/4] Running Salmon quantification..."
    salmon quant -i $SALMON_INDEX \
        -l A \
        -p $CORES \
        --seqBias --gcBias -q \
        -r "$TRIMMED" \
        -o "${OUTDIR}/${line}/" \
        >> $LOG 2>&1
    echo "Done with $line!"
done

echo ""
echo "======================================"
echo "All samples processed successfully!"
echo "Results in: $OUTDIR"
echo "======================================"
