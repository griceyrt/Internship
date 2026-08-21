#!/bin/bash
#SBATCH --job-name=hfd_organoid_comprehensive_qc
#SBATCH --partition=Lake-short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=%j_hfd_organoid_comprehensive_qc.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biotools

# Author: Gricey
#
# Aggregates raw FastQC, fastp reports, trimmed FastQC, and salmon logs into
# one MultiQC report -- raw/trimmed rows sit side by side for direct before/
# after comparison. Mirrors Desmond's final "comprehensive multiqc" step in
# brnaseq_pipeline.sh. Pure re-aggregation, doesn't recompute anything.
#
# Update the *_DIR paths below if the date-stamped folder names differ from
# what's hardcoded here (check with `ls results/qc/` / `ls results/fastp/`).

set -euo pipefail

SCRATCH_DIR="/scratch/Lake/$USER/hfd_organoid_splicing"
RAW_QC_DIR="$SCRATCH_DIR/results/qc/raw_2026-08-03"
TRIMMED_QC_DIR="$SCRATCH_DIR/results/qc/trimmed_2026-08-11"
FASTP_REPORTS_DIR="$SCRATCH_DIR/results/fastp/reports_2026-08-04"
SALMON_DIR="$SCRATCH_DIR/results/salmon"
OUT_DIR="$SCRATCH_DIR/results/qc/comprehensive_$(date +"%Y-%m-%d")"
mkdir -p "$OUT_DIR"

multiqc "$RAW_QC_DIR" "$FASTP_REPORTS_DIR" "$TRIMMED_QC_DIR" "$SALMON_DIR" \
    -o "$OUT_DIR" -n complete_qc_report --flat

echo ""
echo "=== Done. Report: $OUT_DIR/complete_qc_report.html ==="
