#!/bin/bash
#SBATCH --job-name=star_align_count
#SBATCH --partition=Lake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=%j_star_align_count.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biotools

# Author: Gricey
# Aligns each of the 8 GSE130613 trimmed FASTQs with STAR (single-end), then
# counts with htseq-count (union mode, unstranded -- matching Asher's own
# usage). Reuses already-trimmed FASTQs from temp/, no re-download/trim needed.

BASE_DIR="$HOME/orthogonal_validation"
DATA_DIR="$BASE_DIR/data"
INDEX_DIR="$DATA_DIR/star_index_extended"
GTF="$DATA_DIR/transcriptome_productivity_extended.gtf"
OUTDIR="$BASE_DIR/results/GSE130613/SRP194523_star_htseq_$(date +"%Y-%m-%d")"
CORES="$SLURM_CPUS_PER_TASK"

mkdir -p "$OUTDIR"

cat "$BASE_DIR/meta/SRP194523_Acc_List.txt" | while read -r SRR
do
    echo ""
    echo "======================================"
    echo "Processing sample: $SRR"
    echo "======================================"

    TRIMMED="$BASE_DIR/temp/${SRR}_trimmed_full.fastq.gz"
    if [ ! -f "$TRIMMED" ]; then
        echo "ERROR: trimmed FASTQ not found for $SRR at $TRIMMED -- skipping."
        continue
    fi

    SAMPLE_OUT="$OUTDIR/$SRR"
    mkdir -p "$SAMPLE_OUT"

    echo "[1/2] STAR alignment..."
    STAR --runMode alignReads \
        --genomeDir "$INDEX_DIR" \
        --readFilesIn "$TRIMMED" \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix "${SAMPLE_OUT}/${SRR}_" \
        --runThreadN "$CORES"

    BAM="${SAMPLE_OUT}/${SRR}_Aligned.sortedByCoord.out.bam"

    echo "[2/2] htseq-count..."
    htseq-count \
        -f bam \
        -r pos \
        -s no \
        -m union \
        --additional-attr=gene_name \
        "$BAM" \
        "$GTF" \
        > "${SAMPLE_OUT}/${SRR}_counts.txt"

    echo "Done with $SRR."
done

echo ""
echo "======================================"
echo "All samples processed."
echo "Results in: $OUTDIR"
echo "======================================"
