#!/bin/bash
#SBATCH --job-name=rmats_gse130613
#SBATCH --partition=Lake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=%j_rmats_gse130613.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biotools

# Author: Gricey
# rMATS-turbo splicing analysis for GSE130613, on the STAR BAMs from
# star_align_and_count.sh, ORIGINAL (non-extended) GTF, single-end,
# readLength=75. AF/AL event types not covered by rMATS (Pol II-dependent,
# not spliceosome-dependent, so not true AS events -- expected, not a gap).

BASE_DIR="$HOME/orthogonal_validation"
STAR_RESULTS="$BASE_DIR/results/GSE130613/SRP194523_star_htseq_2026-07-17"
GTF_ORIGINAL="$BASE_DIR/data/transcriptome_productivity.gtf"
# OUTDIR/TMPDIR are a fixed path (not date-based) so a resubmitted job after a
# timeout can reuse rMATS's --tmp cache instead of reprocessing every BAM.
OUTDIR="$BASE_DIR/results/GSE130613/SRP194523_rmats_2026-07-20"
TMPDIR="$OUTDIR/tmp"
CORES="$SLURM_CPUS_PER_TASK"

mkdir -p "$OUTDIR" "$TMPDIR"

WT_SAMPLES=(SRR9002567 SRR9002568 SRR9002569 SRR9002570)
KO_SAMPLES=(SRR9002583 SRR9002584 SRR9002585 SRR9002586)

# b1 = WT, b2 = KO/PerDKO -- rMATS's IncLevelDifference (b1-b2, WT-KO) is the
# OPPOSITE sign of DESeq2's log2FC (KO vs WT). Flip when comparing the two.
B1_FILE="$OUTDIR/b1_WT.txt"
B2_FILE="$OUTDIR/b2_KO.txt"

WT_BAMS=""
for SRR in "${WT_SAMPLES[@]}"; do
    BAM="$STAR_RESULTS/$SRR/${SRR}_Aligned.sortedByCoord.out.bam"
    if [ ! -f "$BAM" ]; then
        echo "ERROR: BAM not found for $SRR at $BAM"
        exit 1
    fi
    WT_BAMS="${WT_BAMS}${BAM},"
done
echo "${WT_BAMS%,}" > "$B1_FILE"

KO_BAMS=""
for SRR in "${KO_SAMPLES[@]}"; do
    BAM="$STAR_RESULTS/$SRR/${SRR}_Aligned.sortedByCoord.out.bam"
    if [ ! -f "$BAM" ]; then
        echo "ERROR: BAM not found for $SRR at $BAM"
        exit 1
    fi
    KO_BAMS="${KO_BAMS}${BAM},"
done
echo "${KO_BAMS%,}" > "$B2_FILE"

echo "b1 (WT): $(cat "$B1_FILE")"
echo "b2 (KO): $(cat "$B2_FILE")"

echo ""
echo "Running rMATS-turbo..."
rmats.py \
    --b1 "$B1_FILE" \
    --b2 "$B2_FILE" \
    --gtf "$GTF_ORIGINAL" \
    -t single \
    --readLength 75 \
    --nthread "$CORES" \
    --od "$OUTDIR" \
    --tmp "$TMPDIR"

echo ""
echo "======================================"
echo "rMATS done (or check log above for errors)."
echo "Results in: $OUTDIR"
echo "======================================"
