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

# Author: Gricey
#
# Quantifies each of the 15 trimmed samples against the salmon index.
# Flags match Desmond's brnaseq_pipeline.sh: -l SF (single-end forward
# stranded, confirmed via STAR+RSeQC on his end), --seqBias --gcBias
# --validateMappings.

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
