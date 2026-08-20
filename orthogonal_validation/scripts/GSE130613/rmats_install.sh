#!/bin/bash
# One-time rMATS-turbo install into the biotools conda env on the PSMN cluster.
# Author: Gricey
# Run INTERACTIVELY on cl5218comp1 (login/build node), NOT via sbatch.
#
# Usage:
#   ssh -A -J grangel@ssh.psmn.ens-lyon.fr grangel@allo-psmn.psmn.ens-lyon.fr
#   ssh -A cl5218comp1
#   bash ~/orthogonal_validation/scripts/rmats_install.sh

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biotools

echo "Installing rMATS-turbo into biotools env..."
conda install -y -c bioconda -c conda-forge rmats

echo ""
echo "Checking install..."
rmats.py --version

# NOTE: if conda's solver fails inside the existing biotools env (it already
# has salmon/STAR/htseq/etc, dependency conflicts are possible), fall back to
# a dedicated env instead:
#   conda create -y -n rmats_env -c bioconda -c conda-forge rmats
#   conda activate rmats_env
# and update the `conda activate` line in run_rmats_gse130613.sh to match.
