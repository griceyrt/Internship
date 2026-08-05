#!/bin/bash
#SBATCH --job-name=gse188840_deseq2
#SBATCH --partition=Lake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=%j_gse188840_deseq2.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biotools

Rscript "$HOME/fmr1_polya_splicing/scripts/normalisation_deseq2_oarfish.R"
