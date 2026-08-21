#!/bin/bash
#SBATCH --job-name=oarfish_test
#SBATCH --partition=Lake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=%j_oarfish_test.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biotools

# =============================================================================
# One-sample validation test for oarfish, before committing to a full
# 6-sample run. Confirms oarfish works against the existing Minimap2 BAM
# once re-sorted by name (oarfish requires name/default order, not the
# coordinate-sorted BAM used for Salmon/indexing).
# =============================================================================

BASE_DIR="$HOME/fmr1_polya_splicing"
ALIGN_DIR="$BASE_DIR/results/gse188840_minimap2_2026-08-03/WT_Rep1"
TEST_OUT="$BASE_DIR/results/oarfish_test"
CORES="$SLURM_CPUS_PER_TASK"

mkdir -p "$TEST_OUT"

echo "[1/2] Re-sorting existing BAM by name (oarfish requires name order, not coordinate order)..."
samtools sort -n -@ "$CORES" \
    -o "$TEST_OUT/WT_Rep1.namesorted.bam" \
    "$ALIGN_DIR/WT_Rep1.sorted.bam"

echo "[2/2] Running oarfish..."
oarfish -j "$CORES" \
    -a "$TEST_OUT/WT_Rep1.namesorted.bam" \
    -o "$TEST_OUT/WT_Rep1" \
    --filter-group no-filters \
    --model-coverage

echo ""
echo "======================================"
echo "Done. Check for $TEST_OUT/WT_Rep1.quant"
echo "======================================"
