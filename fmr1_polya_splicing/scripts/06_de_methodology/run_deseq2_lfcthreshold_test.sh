#!/bin/bash
#SBATCH --job-name=deseq2_lfcthreshold_test
#SBATCH --partition=Lake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=%j_deseq2_lfcthreshold_test.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biotools

Rscript "$HOME/fmr1_polya_splicing/scripts/06_de_methodology/deseq2_lfcthreshold_test.R"
