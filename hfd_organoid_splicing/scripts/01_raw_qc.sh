#!/bin/bash
#SBATCH --job-name=hfd_organoid_raw_qc
#SBATCH --partition=Lake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=%j_hfd_organoid_raw_qc.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biotools

# =============================================================================
# Phase 1 — raw read QC (HFD organoid splicing project)
# Author: Gricey
# Description: FastQC on all 15 raw fastq.gz files + MultiQC aggregate
#              report. Mirrors Desmond's brnaseq_pipeline.sh QC step exactly
#              (same tools, same order) so this stays comparable to his run.
#
# STORAGE: raw fastq + QC output go on $SCRATCH, not $HOME (see note in
# 00_download_reference.sh -- $HOME quota is shared lab-wide and too small
# for this data; $SCRATCH has no quota and is where PSMN docs say large
# input/output files belong). Only scripts live in $HOME.
#
# PREREQUISITE: the 15 raw fastq.gz files must already be on the cluster
# under $SCRATCH_DIR/raw_data/fastq/ before submitting this job. They
# currently only live locally (Mac, under
# hfd_organoid_splicing/raw_data/fastq/) -- transfer them up first, e.g.:
#   ssh grangel@ssh.psmn.ens-lyon.fr "mkdir -p /scratch/Lake/grangel/hfd_organoid_splicing/raw_data/fastq"
#   rsync -avP --progress \
#     "/Users/gricey/Desktop/Internship/hfd_organoid_splicing/raw_data/fastq/" \
#     grangel@ssh.psmn.ens-lyon.fr:/scratch/Lake/grangel/hfd_organoid_splicing/raw_data/fastq/
# (run that from your Mac terminal, not on the cluster)
# =============================================================================

set -euo pipefail

SCRATCH_DIR="/scratch/Lake/$USER/hfd_organoid_splicing"
FASTQ_DIR="$SCRATCH_DIR/raw_data/fastq"
QC_DIR="$SCRATCH_DIR/results/qc/raw_$(date +"%Y-%m-%d")"
mkdir -p "$QC_DIR"

echo "=== Running FastQC on raw reads ==="
fastqc -t "$SLURM_CPUS_PER_TASK" -o "$QC_DIR" "$FASTQ_DIR"/*.fastq.gz

echo "=== Aggregating with MultiQC ==="
multiqc "$QC_DIR" -o "$QC_DIR" -n raw_mqc_report

echo ""
echo "=== Done. Report: $QC_DIR/raw_mqc_report.html ==="
