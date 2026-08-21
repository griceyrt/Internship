#!/bin/bash
#SBATCH --job-name=entropy_go_enrichment
#SBATCH --partition=Lake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=%j_entropy_go_enrichment.log

source ~/miniconda3/etc/profile.d/conda.sh

# Uses a separate env from biotools -- org.Mm.eg.db/GO.db conflict with
# biotools' existing R pins. Create once if missing:
#   conda create -n go_entropy -c bioconda -c conda-forge -y \
#       r-base r-ggplot2 bioconductor-org.mm.eg.db bioconductor-go.db
conda activate go_entropy

Rscript "$HOME/fmr1_polya_splicing/scripts/04_go_category_analysis/entropy_go_enrichment_dotplot.R"
