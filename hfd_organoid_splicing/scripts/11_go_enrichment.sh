#!/bin/bash
#SBATCH --job-name=hfd_organoid_go_enrichment
#SBATCH --partition=Lake-short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=%j_hfd_organoid_go_enrichment.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate go_entropy_hfd

# Author: Gricey
#
# Wrapper for 11_go_enrichment.R. Uses go_entropy_hfd, a clone of project
# 3's go_entropy env with clusterProfiler added -- cloned rather than
# installed directly into the shared env, to avoid risking project 3's
# reproducibility.

set -euo pipefail

SCRATCH_DIR="/scratch/Lake/$USER/hfd_organoid_splicing"

Rscript "$HOME/hfd_organoid_splicing/scripts/11_go_enrichment.R" "$SCRATCH_DIR"
