#!/bin/bash
#SBATCH --job-name=htseq_orig_gtf
#SBATCH --partition=Lake-short
#SBATCH --array=0-7
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:30:00
#SBATCH --output=%A_%a_htseq_original_gtf.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biotools

# Array version of htseq_count_original_gtf.sh -- each of the 8 samples runs
# as its own independent SLURM task instead of one long sequential loop.
# Faster to schedule (small, short tasks fit into queue gaps more easily
# than one big multi-hour job) and faster to run (all 8 in parallel instead
# of one after another).
#
# Still reruns htseq-count only, on the already-computed STAR BAM files,
# against the ORIGINAL non-extended GTF. No STAR re-alignment.

BASE_DIR="$HOME/orthogonal_validation"
DATA_DIR="$BASE_DIR/data"
STAR_RESULTS="$BASE_DIR/results/SRP194523_star_htseq_2026-07-17"
GTF_ORIGINAL="$DATA_DIR/transcriptome_productivity.gtf"
OUTDIR="$BASE_DIR/results/SRP194523_star_htseq_original_gtf_$(date +"%Y-%m-%d")"

mkdir -p "$OUTDIR"

# Pick this array task's sample (SLURM_ARRAY_TASK_ID is 0-indexed, sed is 1-indexed)
SRR=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$BASE_DIR/meta/SRP194523_Acc_List.txt")

echo "Array task $SLURM_ARRAY_TASK_ID -> sample $SRR"

BAM="$STAR_RESULTS/$SRR/${SRR}_Aligned.sortedByCoord.out.bam"
if [ ! -f "$BAM" ]; then
    echo "ERROR: BAM not found for $SRR at $BAM"
    exit 1
fi

SAMPLE_OUT="$OUTDIR/$SRR"
mkdir -p "$SAMPLE_OUT"

echo "htseq-count (original GTF)..."
htseq-count \
    -f bam \
    -r pos \
    -s no \
    -m union \
    --additional-attr=gene_name \
    "$BAM" \
    "$GTF_ORIGINAL" \
    > "${SAMPLE_OUT}/${SRR}_counts.txt"

echo "Done with $SRR."
