#!/bin/bash
#SBATCH --job-name=htseq_original_gtf
#SBATCH --partition=Lake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=%j_htseq_original_gtf.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biotools

# Author: Gricey
# Reruns htseq-count on the already-computed STAR BAM files (star_align_and_count.sh)
# against the ORIGINAL non-extended GTF, to isolate whether the high "ambiguous"
# rate was inflated by the 3'UTR/5' boundary extension specifically.

BASE_DIR="$HOME/orthogonal_validation"
DATA_DIR="$BASE_DIR/data"
STAR_RESULTS="$BASE_DIR/results/GSE130613/SRP194523_star_htseq_2026-07-17"
GTF_ORIGINAL="$DATA_DIR/transcriptome_productivity.gtf"
OUTDIR="$BASE_DIR/results/GSE130613/SRP194523_star_htseq_original_gtf_$(date +"%Y-%m-%d")"

mkdir -p "$OUTDIR"

cat "$BASE_DIR/meta/SRP194523_Acc_List.txt" | while read -r SRR
do
    echo ""
    echo "======================================"
    echo "Processing sample: $SRR"
    echo "======================================"

    BAM="$STAR_RESULTS/$SRR/${SRR}_Aligned.sortedByCoord.out.bam"
    if [ ! -f "$BAM" ]; then
        echo "ERROR: BAM not found for $SRR at $BAM -- skipping."
        continue
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
done

echo ""
echo "======================================"
echo "All samples processed."
echo "Results in: $OUTDIR"
echo "======================================"
