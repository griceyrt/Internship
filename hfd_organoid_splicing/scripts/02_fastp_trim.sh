#!/bin/bash
#SBATCH --job-name=hfd_organoid_fastp_trim
#SBATCH --partition=Lake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=%j_hfd_organoid_fastp_trim.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biotools

# Author: Gricey
#
# fastp trimming, all 15 raw fastq.gz files. Flags copied exactly from
# Desmond's brnaseq_pipeline.sh (poly-G removal, quality sliding window,
# adapter auto-detect, <35bp discard) so splicing results stay comparable
# to his DE/proteomics work on the same samples.
# Sequential loop, not GNU parallel (not installed in biotools) -- fastp
# gets multiple threads per sample instead.

set -euo pipefail

SCRATCH_DIR="/scratch/Lake/$USER/hfd_organoid_splicing"
FASTQ_DIR="$SCRATCH_DIR/raw_data/fastq"
TRIM_DIR="$SCRATCH_DIR/results/fastp/trimmed"
REPORT_DIR="$SCRATCH_DIR/results/fastp/reports_$(date +"%Y-%m-%d")"
mkdir -p "$TRIM_DIR" "$REPORT_DIR"

for fastq_file in "$FASTQ_DIR"/*.fastq.gz; do
    sample_name="$(basename "${fastq_file%_R1_001.fastq.gz}")"
    outfile="$TRIM_DIR/${sample_name}_trimmed.fastq.gz"

    echo "=== fastp processing ${sample_name} ==="
    fastp \
        -i "$fastq_file" \
        -o "$outfile" \
        --trim_poly_g \
        --qualified_quality_phred 28 \
        --length_required 35 \
        --cut_front --cut_tail \
        --cut_window_size 4 \
        --cut_mean_quality 28 \
        --thread "$SLURM_CPUS_PER_TASK" \
        -h "$REPORT_DIR/${sample_name}.html" \
        -j "$REPORT_DIR/${sample_name}.json"
done

echo ""
echo "=== Done. Trimmed reads: $TRIM_DIR ==="
echo "=== Per-sample fastp reports: $REPORT_DIR ==="
ls -lh "$TRIM_DIR"
