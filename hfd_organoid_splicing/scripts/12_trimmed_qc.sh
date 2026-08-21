#!/bin/bash
#SBATCH --job-name=hfd_organoid_trimmed_qc
#SBATCH --partition=Lake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=%j_hfd_organoid_trimmed_qc.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biotools

# Author: Gricey
#
# FastQC + MultiQC on the 15 trimmed fastq files (mirrors 01_raw_qc.sh).
# --flat needed for MultiQC to render plots in this env -- if it comes back
# blank, `conda install -c conda-forge python-kaleido -y` first.

set -euo pipefail

SCRATCH_DIR="/scratch/Lake/$USER/hfd_organoid_splicing"
TRIM_DIR="$SCRATCH_DIR/results/fastp/trimmed"
QC_DIR="$SCRATCH_DIR/results/qc/trimmed_$(date +"%Y-%m-%d")"
mkdir -p "$QC_DIR"

echo "=== Running FastQC on trimmed reads ==="
fastqc -t "$SLURM_CPUS_PER_TASK" -o "$QC_DIR" "$TRIM_DIR"/*_trimmed.fastq.gz

echo "=== Aggregating with MultiQC ==="
multiqc "$QC_DIR" -o "$QC_DIR" -n trimmed_mqc_report --flat

echo ""
echo "=== Done. Report: $QC_DIR/trimmed_mqc_report.html ==="
