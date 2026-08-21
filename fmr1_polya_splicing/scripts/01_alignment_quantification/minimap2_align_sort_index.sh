#!/bin/bash
#SBATCH --job-name=gse188840_align
#SBATCH --partition=Lake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --time=10:00:00
#SBATCH --output=%j_gse188840_align.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biotools

BASE_DIR="$HOME/fmr1_polya_splicing"
TEMP_DIR="$BASE_DIR/temp"
TRANSCRIPTOME="$BASE_DIR/data/reference/Mus_musculus.GRCm39.cdna.all.fa.gz"
OUTDIR="$BASE_DIR/results/gse188840_minimap2_$(date +"%Y-%m-%d")"
CORES="$SLURM_CPUS_PER_TASK"

mkdir -p "$OUTDIR"

awk '{print $2}' "$BASE_DIR/meta/GSE188840_Acc_List.txt" | while read -r SAMPLE
do
    echo ""
    echo "======================================"
    echo "Processing sample: $SAMPLE"
    echo "======================================"

    FASTQ="$TEMP_DIR/${SAMPLE}.fastq.gz"
    if [ ! -f "$FASTQ" ]; then
        echo "ERROR: fastq not found for $SAMPLE at $FASTQ -- skipping."
        continue
    fi

    SAMPLE_OUT="$OUTDIR/$SAMPLE"
    mkdir -p "$SAMPLE_OUT"

    # Skip if this sample already finished -- lets the job be safely
    # resubmitted after a cancel/timeout without redoing completed samples.
    if [ -f "${SAMPLE_OUT}/${SAMPLE}.sorted.bam.bai" ]; then
        echo "$SAMPLE already done (index found) -- skipping."
        continue
    fi

    echo "[1/3] Minimap2 alignment..."
    minimap2 -ax map-ont -N 100 -t "$CORES" "$TRANSCRIPTOME" "$FASTQ" \
        | samtools view -b -o "${SAMPLE_OUT}/${SAMPLE}.unsorted.bam" -

    echo "[2/3] Samtools sort..."
    samtools sort -@ "$CORES" -o "${SAMPLE_OUT}/${SAMPLE}.sorted.bam" "${SAMPLE_OUT}/${SAMPLE}.unsorted.bam"

    echo "[3/3] Samtools index..."
    samtools index "${SAMPLE_OUT}/${SAMPLE}.sorted.bam"

    echo "Done with $SAMPLE."
done

echo ""
echo "======================================"
echo "All samples processed."
echo "Results in: $OUTDIR"
echo "======================================"