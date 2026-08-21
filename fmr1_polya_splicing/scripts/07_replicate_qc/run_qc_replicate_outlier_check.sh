#!/bin/bash
#SBATCH --job-name=qc_replicate_check
#SBATCH --partition=Lake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=%j_qc_replicate_check.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biotools

Rscript "$HOME/fmr1_polya_splicing/scripts/07_replicate_qc/qc_replicate_outlier_check.R"
