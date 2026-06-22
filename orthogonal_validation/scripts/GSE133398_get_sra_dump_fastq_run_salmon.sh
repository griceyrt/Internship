#!/bin/bash
# =============================================================================
# GSE133398 — SRSF3 KO hepatocytes (Illumina HiSeq 2500, paired-end)
# Adapted from: get_sra_dump_fastq_run_salmon.sh (GSE130613 single-end version)
#
# Key differences from GSE130613 script:
#   - PAIRED-END: fasterq-dump --split-files produces _1.fastq + _2.fastq
#   - Salmon: -1/-2 flags instead of -r
#   - TruSeq Stranded mRNA: Salmon auto-detects ISR (stranded reverse)
#   - HiSeq 2500 (not NextSeq 500): NO polyG artifact → trimming is adapter only
#   - Full-length mRNA: NO polyA/T tail trimming needed (not 3' end biased)
#
# SRP accession : SRP212166
# Groups:
#   CRE (SRSF3 KO): SRR9606468, SRR9606469, SRR9606470, SRR9606471
#   GFP (control) : SRR9606472, SRR9606473, SRR9606474, SRR9606475
#
# USAGE: bash scripts/GSE133398_get_sra_dump_fastq_run_salmon.sh [path/to/salmon/index]
# Run from: orthogonal_validation/
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${SCRIPT_DIR}/.."

SRP="SRP212166"
SALMON_INDEX="${BASE_DIR}/data/transcriptome_ext_index"
OUTDIR="${BASE_DIR}/results/${SRP}_salmon_$(date +"%Y-%m-%d")"
LOG="${OUTDIR}/${SRP}.log"
CORES=7

# Use existing index, or build it if path provided as $1
if [ ! -z "$1" ]; then
    SALMON_INDEX="$1"
fi

mkdir -p "${OUTDIR}"
echo "Run on $(date)" > "$LOG"
salmon --version >> "$LOG"
fasterq-dump --version >> "$LOG"

cat "${BASE_DIR}/meta/${SRP}_Acc_List.txt" | while read line
do
    echo ""
    echo "======================================"
    echo "Processing sample: $line"
    echo "======================================"

    # Fetch the .sra file
    echo "[1/4] Prefetching $line..."
    prefetch -O "${BASE_DIR}/raw_data/" "$line" >> "$LOG"

    # Extract PAIRED fastq files (_1 and _2)
    R1="${BASE_DIR}/temp/${line}_1.fastq.gz"
    R2="${BASE_DIR}/temp/${line}_2.fastq.gz"

    if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
        echo "[2/4] Converting .sra to paired .fastq.gz..."
        fasterq-dump --threads 6 \
            --split-files \
            --outdir "${BASE_DIR}/temp/" \
            "${BASE_DIR}/raw_data/${line}/${line}.sra" \
            >> "$LOG" 2>&1
        gzip "${BASE_DIR}/temp/${line}_1.fastq"
        gzip "${BASE_DIR}/temp/${line}_2.fastq"
    else
        echo "[2/4] FASTQ files already exist, skipping fasterq-dump."
    fi

    # Trim Illumina TruSeq adapter only — no polyG/polyA/T needed (HiSeq, full-length mRNA)
    TRIMMED_R1="${BASE_DIR}/temp/${line}_trimmed_R1.fastq.gz"
    TRIMMED_R2="${BASE_DIR}/temp/${line}_trimmed_R2.fastq.gz"

    if [ ! -f "$TRIMMED_R1" ] || [ ! -f "$TRIMMED_R2" ]; then
        echo "[3/4] Trimming Illumina adapters with cutadapt (paired-end)..."
        cutadapt \
            -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
            -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
            --minimum-length 20 \
            -j "$CORES" \
            -o "$TRIMMED_R1" \
            -p "$TRIMMED_R2" \
            "$R1" "$R2" \
            >> "$LOG" 2>&1
    else
        echo "[3/4] Trimmed files already exist, skipping cutadapt."
    fi

    # Salmon quantification — paired-end (-1/-2), auto library detection
    echo "[4/4] Running Salmon quantification (paired-end)..."
    SALMON_START=$(date +%s)
    salmon quant -i "$SALMON_INDEX" \
        -l A \
        -p "$CORES" \
        --seqBias --gcBias \
        -1 "$TRIMMED_R1" \
        -2 "$TRIMMED_R2" \
        -o "${OUTDIR}/${line}/" \
        1>> "$LOG"
    SALMON_END=$(date +%s)
    SALMON_ELAPSED=$(( (SALMON_END - SALMON_START) / 60 ))
    echo "Done with $line! Salmon took ${SALMON_ELAPSED} minutes."
    echo "Check lib_format_counts.json in ${OUTDIR}/${line}/ to verify ISR detection."
done

echo ""
echo "======================================"
echo "All samples processed successfully!"
echo "Results in: $OUTDIR"
echo "======================================"
echo ""
echo "IMPORTANT: Check Salmon's detected library type in each sample's"
echo "  aux_info/lib_format_counts.json — expected: ISR (stranded reverse)"
echo "  for TruSeq Stranded mRNA libraries."
