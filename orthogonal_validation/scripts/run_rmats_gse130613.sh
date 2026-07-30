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

# rMATS-turbo splicing analysis for GSE130613, run directly on the STAR BAMs
# already produced by star_align_and_count.sh (2026-07-17). Uses the ORIGINAL
# (non-extended) GTF, matching the adopted DE pipeline (deseq2_star_htseq.R,
# 390 sig genes). Single-end mode, readLength=75 (confirmed from SRA read
# headers, e.g. "@SRR9002584.1 1 length=75").
#
# Approved by Bharath/Kiran 2026-07-20 after a tool comparison (rMATS vs
# LeafCutter vs MAJIQ vs DEXSeq). Known caveats going in, do not be alarmed by:
#  - AF/AL event types are not covered by rMATS. Kiran: fine, those aren't
#    true alternative-splicing events anyway, they're Pol II-dependent
#    (transcription start/end site choice), not spliceosome-dependent.
#  - Bharath: comparison to the old SUPPA numbers won't be cleanly
#    interpretable (method changed AND data is 3'-biased, both at once) --
#    running this mainly so we can tell the reviewer we tried.
#  - Method may return very few events or effectively collapse if
#    junction-spanning reads are too sparse given the 3' bias -- expected
#    possible outcome per Bharath, not necessarily a bug.

BASE_DIR="$HOME/orthogonal_validation"
STAR_RESULTS="$BASE_DIR/results/SRP194523_star_htseq_2026-07-17"
GTF_ORIGINAL="$BASE_DIR/data/transcriptome_productivity.gtf"
# NOTE: OUTDIR/TMPDIR are a FIXED path (not $(date ...)) on purpose. rMATS
# caches each BAM's "prep" output into --tmp; if this job times out partway
# through and gets resubmitted, pointing at the SAME --tmp lets rMATS detect
# already-processed BAMs and skip redoing them, only finishing the rest. A
# date-based path would silently break that resumability on any rerun that
# crosses midnight. Bump the suffix manually only if starting a genuinely
# fresh/unrelated run.
OUTDIR="$BASE_DIR/results/SRP194523_rmats_2026-07-20"
TMPDIR="$OUTDIR/tmp"
CORES="$SLURM_CPUS_PER_TASK"

mkdir -p "$OUTDIR" "$TMPDIR"

WT_SAMPLES=(SRR9002567 SRR9002568 SRR9002569 SRR9002570)
KO_SAMPLES=(SRR9002583 SRR9002584 SRR9002585 SRR9002586)

# Build comma-separated BAM path lists.
# b1 = WT, b2 = KO/PerDKO -- matches deseq2_star_htseq.R's contrast direction
# (contrast = c("condition", "KO", "WT")), so rMATS's IncLevelDifference sign
# convention (b1 - b2, i.e. WT - KO) is the OPPOSITE sign of our DESeq2
# log2FC (KO vs WT). Flip the sign when comparing/plotting the two side by
# side, or re-derive dPSI as (b2 - b1) downstream to match.
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
