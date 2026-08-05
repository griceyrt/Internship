#!/bin/bash
#SBATCH --job-name=hfd_organoid_salmon_quant
#SBATCH --partition=Lake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=%j_hfd_organoid_salmon_quant.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biotools

# =============================================================================
# Phase 5b — salmon quantification (HFD organoid splicing project)
# Author: Gricey
# Description: Quantifies each of the 15 trimmed samples against the
#              release-116 salmon index (03_salmon_index.sh). Flags match
#              Desmond's brnaseq_pipeline.sh exactly: library type SF
#              (single-end forward stranded, confirmed by Desmond via STAR +
#              RSeQC on Lib1 -- we did not re-run this check ourselves, see
#              notes/pipeline_plan.txt Phase 4), --seqBias --gcBias
#              --validateMappings.
#
# Sequential loop over samples (same reasoning as 02_fastp_trim.sh -- no
# proven GNU parallel in this env), salmon given all cpus per sample instead.
#
# STORAGE: everything on $SCRATCH, only this script lives in $HOME.
# =============================================================================

set -euo pipefail

SCRATCH_DIR="/scratch/Lake/$USER/hfd_organoid_splicing"
TRIM_DIR="$SCRATCH_DIR/results/fastp/trimmed"
INDEX_DIR="$SCRATCH_DIR/data/salmon_index"
QUANT_DIR="$SCRATCH_DIR/results/salmon"
mkdir -p "$QUANT_DIR"

for trimmed_fastq in "$TRIM_DIR"/*_trimmed.fastq.gz; do
    sample_name="$(basename "${trimmed_fastq%_trimmed.fastq.gz}")"
    outdir="$QUANT_DIR/${sample_name}"

    echo "=== salmon quant: ${sample_name} ==="
    salmon quant \
        -i "$INDEX_DIR" \
        -l SF \
        -r "$trimmed_fastq" \
        --seqBias --gcBias \
        --validateMappings \
        -p "$SLURM_CPUS_PER_TASK" \
        -o "$outdir"
done

echo ""
echo "=== Done. Per-sample quant.sf files under: $QUANT_DIR ==="
for d in "$QUANT_DIR"/*/; do
    echo -n "$(basename "$d"): "
    [ -f "${d}quant.sf" ] && echo "OK ($(wc -l < "${d}quant.sf") lines)" || echo "MISSING quant.sf"
done
