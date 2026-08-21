#!/bin/bash
#SBATCH --job-name=entropy_by_category
#SBATCH --partition=Lake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=%j_entropy_by_category.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biotools

# One-time prep: download the mouse housekeeping gene list if not present.
# Pulled from HRT Atlas's GitHub repo directly (the site's own download
# page is a dynamic app, not a stable static URL).
HK_FILE="$HOME/fmr1_polya_splicing/data/reference/Housekeeping_Genes_Mouse.RData"
if [ ! -f "$HK_FILE" ] || [ ! -s "$HK_FILE" ]; then
    echo "Downloading HRT Atlas mouse housekeeping gene RData..."
    rm -f "$HK_FILE"  # clear out any empty/partial file from an earlier failed attempt
    wget -O "$HK_FILE" https://raw.githubusercontent.com/Bidossessih/HRT_Atlas/master/www/Housekeeping_Genes_Mouse.RData
fi

Rscript "$HOME/fmr1_polya_splicing/scripts/04_go_category_analysis/entropy_by_gene_category.R"
