#!/bin/bash
#SBATCH --job-name=gse188840_download
#SBATCH --partition=Lake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=%j_gse188840_download.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biotools

# =============================================================================
# Downloads the 6 GSE188840 ONT direct RNA-seq runs (WT vs Fmr1 KO mouse
# brain cortex) via prefetch + fasterq-dump + pigz. Single-end (no
# --split-files) since this is nanopore, not paired-end Illumina.
#
# The "inspect archive" step logs each .sra file's size/contents so it's
# easy to tell whether the raw fast5 signal is bundled in or fastq-only.
# =============================================================================

BASE_DIR="$HOME/fmr1_polya_splicing"
RAW_DIR="$BASE_DIR/raw_data"
TEMP_DIR="$BASE_DIR/temp"
LOG="$BASE_DIR/results/GSE188840_download_$(date +"%Y-%m-%d").log"
CORES="$SLURM_CPUS_PER_TASK"

mkdir -p "$RAW_DIR" "$TEMP_DIR" "$(dirname "$LOG")"

echo "Run on $(date)" > "$LOG"
prefetch --version >> "$LOG"
fasterq-dump --version >> "$LOG"
pigz --version >> "$LOG" 2>&1

cat "$BASE_DIR/meta/GSE188840_Acc_List.txt" | while read -r SRX SAMPLE
do
    echo ""
    echo "======================================"
    echo "Processing $SAMPLE ($SRX)"
    echo "======================================"

    # prefetch names the output folder after the resolved RUN accession
    # (SRRxxxxxxx), not the SRX we searched for -- capture it from stdout.
    PREFETCH_OUT=$(prefetch -O "$RAW_DIR" "$SRX" 2>&1)
    echo "$PREFETCH_OUT" >> "$LOG"
    SRR=$(echo "$PREFETCH_OUT" | grep -oE 'SRR[0-9]+' | head -n 1)

    if [ -z "$SRR" ]; then
        echo "ERROR: could not determine run accession for $SAMPLE ($SRX) -- skipping" >> "$LOG"
        continue
    fi

    echo "[inspect] $SAMPLE resolved to run $SRR, archive contents:" >> "$LOG"
    du -sh "$RAW_DIR/$SRR" >> "$LOG" 2>&1
    find "$RAW_DIR/$SRR" -type f -exec ls -lh {} \; >> "$LOG" 2>&1

    # fasterq-dump has no built-in --gzip, so compress with pigz afterwards.
    if [ ! -f "$TEMP_DIR/${SAMPLE}.fastq.gz" ]; then
        fasterq-dump --threads "$CORES" \
            -O "$TEMP_DIR" \
            "$RAW_DIR/$SRR" \
            >> "$LOG" 2>&1

        # fasterq-dump names output after the run accession, not our sample
        # label -- find it, compress, and rename to match.
        DUMPED=$(find "$TEMP_DIR" -maxdepth 1 -name "*.fastq" | head -n 1)
        if [ -n "$DUMPED" ]; then
            pigz -p "$CORES" "$DUMPED"
            mv "${DUMPED}.gz" "$TEMP_DIR/${SAMPLE}.fastq.gz"
        else
            echo "WARNING: no fastq found for $SAMPLE after fasterq-dump -- check log" >> "$LOG"
        fi
    fi

    echo "Done with $SAMPLE."
done

echo ""
echo "======================================"
echo "All samples processed. Check $LOG for archive sizes."
echo "If any raw_data/<SRR> folder is only a few GB, that sample is almost"
echo "certainly fastq-only -- flag this before assuming Nanopolish can run."
echo "======================================"
