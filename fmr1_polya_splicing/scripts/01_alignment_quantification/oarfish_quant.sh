#!/bin/bash
#SBATCH --job-name=gse188840_oarfish
#SBATCH --partition=Lake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=05:00:00
#SBATCH --output=%j_gse188840_oarfish.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biotools

# =============================================================================
# Full 6-sample oarfish quantification. Uses oarfish instead of Salmon,
# since Salmon's own --ont flag redirects ONT users to oarfish. Runs
# against the existing Minimap2 alignment, re-sorted by name (oarfish
# requires name/default order, not the coordinate-sorted BAM used for
# indexing). tximport has native support (type = "oarfish").
# =============================================================================

BASE_DIR="$HOME/fmr1_polya_splicing"
ALIGN_BASE="$BASE_DIR/results/gse188840_minimap2_2026-08-03"
OUTDIR="$BASE_DIR/results/gse188840_oarfish_$(date +"%Y-%m-%d")"
CORES="$SLURM_CPUS_PER_TASK"

mkdir -p "$OUTDIR"

awk '{print $2}' "$BASE_DIR/meta/GSE188840_Acc_List.txt" | while read -r SAMPLE
do
    echo ""
    echo "======================================"
    echo "Processing sample: $SAMPLE"
    echo "======================================"

    ALIGNED_BAM="$ALIGN_BASE/$SAMPLE/${SAMPLE}.sorted.bam"
    if [ ! -f "$ALIGNED_BAM" ]; then
        echo "ERROR: aligned BAM not found for $SAMPLE at $ALIGNED_BAM -- skipping."
        continue
    fi

    # Skip if this sample already finished
    if [ -f "${OUTDIR}/${SAMPLE}.quant" ]; then
        echo "$SAMPLE already done (quant file found) -- skipping."
        continue
    fi

    NAMESORTED_BAM="${OUTDIR}/${SAMPLE}.namesorted.bam"

    echo "[1/2] Re-sorting by name (oarfish requires name/default order)..."
    samtools sort -n -@ "$CORES" -o "$NAMESORTED_BAM" "$ALIGNED_BAM"

    echo "[2/2] Running oarfish..."
    oarfish -j "$CORES" \
        -a "$NAMESORTED_BAM" \
        -o "${OUTDIR}/${SAMPLE}" \
        --filter-group no-filters \
        --model-coverage

    # Clean up the intermediate name-sorted BAM immediately -- it's only
    # ever needed for the oarfish call above, no reason to keep 6 copies
    # of a multi-GB file around (this bit us with the unsorted Minimap2
    # BAMs earlier).
    echo "Cleaning up intermediate name-sorted BAM..."
    rm -f "$NAMESORTED_BAM"

    echo "Done with $SAMPLE."
done

echo ""
echo "======================================"
echo "All samples processed."
echo "Results in: $OUTDIR"
echo "======================================"
