#!/bin/bash
#SBATCH --job-name=salmon_GSE130613_extended
#SBATCH --partition=Lake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=%j_salmon_GSE130613_extended.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biotools

cd ~/orthogonal_validation
bash scripts/get_sra_dump_fastq_run_salmon.sh SRP194523 ~/orthogonal_validation/data/transcriptome_ext_extended_index
