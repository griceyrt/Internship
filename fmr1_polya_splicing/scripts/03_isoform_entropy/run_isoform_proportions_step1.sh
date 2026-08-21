#!/bin/bash
#SBATCH --job-name=isoform_prop_step1
#SBATCH --partition=Lake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=%j_isoform_prop_step1.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biotools

Rscript "$HOME/fmr1_polya_splicing/scripts/03_isoform_entropy/isoform_proportions_step1.R"
