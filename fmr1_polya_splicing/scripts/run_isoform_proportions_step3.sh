#!/bin/bash
#SBATCH --job-name=isoform_prop_step3
#SBATCH --partition=Lake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=%j_isoform_prop_step3.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biotools

Rscript "$HOME/fmr1_polya_splicing/scripts/isoform_proportions_step3.R"
