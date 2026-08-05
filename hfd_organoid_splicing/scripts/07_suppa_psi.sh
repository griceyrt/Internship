#!/bin/bash
#SBATCH --job-name=hfd_organoid_suppa_psi
#SBATCH --partition=Lake-short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=%j_hfd_organoid_suppa_psi.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biotools

# =============================================================================
# Phase 7 — SUPPA2 psiPerEvent (HFD organoid splicing project)
# Author: Gricey
# Description: Computes PSI (Percent Spliced In) per event per sample, for
#              each of the 7 event types, using the combined TPM matrix
#              (06_build_tpm_matrix.sh) against the .ioe event files
#              (05_suppa_generate_events.sh).
#
# Loops over all *.ioe files rather than hardcoding the 7 event type names,
# so it automatically adapts if the event set ever changes.
# =============================================================================

set -euo pipefail

CODE_DIR="$HOME/hfd_organoid_splicing"
SCRATCH_DIR="/scratch/Lake/$USER/hfd_organoid_splicing"
EVENTS_DIR="$SCRATCH_DIR/results/suppa/events"
TPM_MATRIX="$SCRATCH_DIR/results/suppa/tpm_matrix.tab"
PSI_DIR="$SCRATCH_DIR/results/suppa/psi"
mkdir -p "$PSI_DIR"

SUPPA="python3 $CODE_DIR/SUPPA/suppa.py"

for ioefile in "$EVENTS_DIR"/*.ioe; do
    event_name="$(basename "$ioefile" .ioe)"
    echo "=== psiPerEvent: ${event_name} ==="
    $SUPPA psiPerEvent \
        -i "$ioefile" \
        -e "$TPM_MATRIX" \
        -o "$PSI_DIR/${event_name}"
done

echo ""
echo "=== Done. PSI files: $PSI_DIR ==="
ls -lh "$PSI_DIR"
