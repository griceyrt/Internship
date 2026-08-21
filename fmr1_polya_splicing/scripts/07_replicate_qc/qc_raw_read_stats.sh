#!/bin/bash
#SBATCH --job-name=qc_raw_read_stats
#SBATCH --partition=Lake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=%j_qc_raw_read_stats.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biotools

# =============================================================================
# Independent, pipeline-free replicate check: raw fastq read counts +
# Minimap2 mapping rates, straight from the original fastq/BAM files
# already on disk. No new alignment or quantification -- upstream of both
# oarfish and DESeq2, so any outlier signal here is independent of
# quantifier/DESeq2-level/threshold choices made downstream.
# =============================================================================

BASE_DIR="$HOME/fmr1_polya_splicing"
SAMPLES=(WT_Rep1 WT_Rep2 WT_Rep3 KO_Rep1 KO_Rep2 KO_Rep3)

echo "=== Diagnostic: actual folder contents ==="
echo "-- temp/ (following symlink) --"
find -L "$BASE_DIR/temp" -maxdepth 2 2>&1
echo ""
echo "-- results/gse188840_minimap2_2026-08-03/ (following symlink) --"
find -L "$BASE_DIR/results/gse188840_minimap2_2026-08-03" -maxdepth 2 2>&1
echo ""
OUT="$BASE_DIR/results/qc_replicate_check/raw_read_stats.tsv"
mkdir -p "$(dirname "$OUT")"

echo -e "sample\traw_reads\ttotal_bam_reads\tmapped_reads\tmapping_pct" > "$OUT"

for SAMPLE in "${SAMPLES[@]}"; do
    echo "=== $SAMPLE ==="

    FASTQ=$(find -L "$BASE_DIR/temp" -iname "*${SAMPLE}*.fastq.gz" 2>/dev/null | head -1)
    if [ -z "$FASTQ" ]; then
        echo "  WARNING: no fastq found for $SAMPLE under temp/"
        RAW_READS="NA"
    else
        echo "  Counting raw reads in $FASTQ ..."
        RAW_READS=$(($(zcat "$FASTQ" | wc -l) / 4))
        echo "  Raw reads: $RAW_READS"
    fi

    BAM=$(find -L "$BASE_DIR/results/gse188840_minimap2_2026-08-03" -iname "*${SAMPLE}*sorted.bam" 2>/dev/null | head -1)
    if [ -z "$BAM" ]; then
        echo "  WARNING: no sorted BAM found for $SAMPLE"
        TOTAL="NA"; MAPPED="NA"; PCT="NA"
    else
        echo "  Running samtools flagstat on $BAM ..."
        FLAGSTAT=$(samtools flagstat "$BAM")
        TOTAL=$(echo "$FLAGSTAT" | grep "in total" | cut -d' ' -f1)
        MAPPED=$(echo "$FLAGSTAT" | grep " mapped (" | head -1 | cut -d' ' -f1)
        PCT=$(echo "$FLAGSTAT" | grep " mapped (" | head -1 | grep -oP '\(\K[0-9.]+(?=%)')
        echo "  Total: $TOTAL, Mapped: $MAPPED ($PCT%)"
    fi

    echo -e "${SAMPLE}\t${RAW_READS}\t${TOTAL}\t${MAPPED}\t${PCT}" >> "$OUT"
    echo ""
done

echo "=== DONE ==="
echo "Saved: $OUT"
cat "$OUT"
