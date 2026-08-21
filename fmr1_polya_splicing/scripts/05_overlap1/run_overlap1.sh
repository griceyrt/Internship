#!/bin/bash
#SBATCH --job-name=gse188840_overlap1
#SBATCH --partition=Lake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=%j_gse188840_overlap1.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biotools

Rscript "$HOME/fmr1_polya_splicing/scripts/05_overlap1/overlap1_de_vs_dpsi.R"
