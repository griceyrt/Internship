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

# Author: Gricey
#
# FastQC + MultiQC on the 15 raw fastq.gz files. Mirrors Desmond's own QC
# step so results stay comparable to his run.
#
# Requires the 15 fastq.gz files already synced to $SCRATCH_DIR/raw_data/fastq/,
# e.g. from a local Mac:
#   rsync -avP raw_data/fastq/ cl5218comp1:/scratch/Lake/grangel/hfd_organoid_splicing/raw_data/fastq/
# (assumes an SSH config alias "cl5218comp1" set up for the gateway chain)

set -euo pipefail

SCRATCH_DIR="/scratch/Lake/$USER/hfd_organoid_splicing"
FASTQ_DIR="$SCRATCH_DIR/raw_data/fastq"
QC_DIR="$SCRATCH_DIR/results/qc/raw_$(date +"%Y-%m-%d")"
mkdir -p "$QC_DIR"

echo "=== Running FastQC on raw reads ==="
fastqc -t "$SLURM_CPUS_PER_TASK" -o "$QC_DIR" "$FASTQ_DIR"/*.fastq.gz

echo "=== Aggregating with MultiQC ==="
multiqc "$QC_DIR" -o "$QC_DIR" -n raw_mqc_report --flat

echo ""
echo "=== Done. Report: $QC_DIR/raw_mqc_report.html ==="
