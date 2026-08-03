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
# Download script for GSE188840 (Uttam et al. 2022, PMID 35217597)
# Author: Gricey
# Description: Pulls the 6 ONT direct RNA-seq runs (WT vs Fmr1 KO mouse brain
#              cortex, PromethION, SQK-RNA002) for project 3. Uses SRA
#              accessions resolved from GEO GSM pages (see
#              meta/GSE188840_Acc_List.txt). Uses prefetch + fasterq-dump +
#              pigz (biotools env has these, not parallel-fastq-dump, which
#              BHARATH_get_sra_dump_fastq_run_salmon.sh in orthogonal_validation
#              used instead -- fasterq-dump is the modern official sra-tools
#              dumper anyway). Single-end (no --split-files) since this is
#              nanopore long-read data, not paired-end Illumina.
#
# IMPORTANT: this only confirms whether fastq is retrievable. Whether the
# raw fast5/pod5 signal (needed for Nanopolish poly(A) estimation in part 2)
# is bundled in the SRA object is NOT confirmed yet -- that's what the
# "inspect archive" step below is for. If the .sra files come back only a
# few GB each, they are almost certainly fastq-only and we'll need to ask
# Kiran/Lies where the raw signal lives.
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

    # Fetch the .sra object if not already present. NOTE: prefetch resolves
    # the SRX to its underlying run and names the output folder after the
    # RUN accession (SRRxxxxxxx), not the SRX we searched for -- capture that
    # from prefetch's own output rather than assuming $RAW_DIR/$SRX exists.
    PREFETCH_OUT=$(prefetch -O "$RAW_DIR" "$SRX" 2>&1)
    echo "$PREFETCH_OUT" >> "$LOG"
    SRR=$(echo "$PREFETCH_OUT" | grep -oE 'SRR[0-9]+' | head -n 1)

    if [ -z "$SRR" ]; then
        echo "ERROR: could not determine run accession for $SAMPLE ($SRX) -- skipping" >> "$LOG"
        continue
    fi

    # Inspect what we actually got (size is the quick tell for fastq-only vs.
    # fastq+signal) -- logged so we can eyeball it once the job finishes
    echo "[inspect] $SAMPLE resolved to run $SRR, archive contents:" >> "$LOG"
    du -sh "$RAW_DIR/$SRR" >> "$LOG" 2>&1
    find "$RAW_DIR/$SRR" -type f -exec ls -lh {} \; >> "$LOG" 2>&1

    # Extract fastq (single-end -- no --split-files, this is nanopore not
    # paired-end Illumina). fasterq-dump has no built-in --gzip, so compress
    # separately with pigz afterwards.
    if [ ! -f "$TEMP_DIR/${SAMPLE}.fastq.gz" ]; then
        fasterq-dump --threads "$CORES" \
            -O "$TEMP_DIR" \
            "$RAW_DIR/$SRR" \
            >> "$LOG" 2>&1

        # fasterq-dump names output after the run accession (e.g. SRRxxxxx.fastq),
        # not our SRX/sample label -- find it, compress, and rename
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
