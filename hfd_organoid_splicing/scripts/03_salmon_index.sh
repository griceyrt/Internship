#!/bin/bash
#SBATCH --job-name=hfd_organoid_salmon_index
#SBATCH --partition=Lake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=%j_hfd_organoid_salmon_index.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biotools

# Author: Gricey
#
# Decoy-aware salmon index (genome as decoy, cDNA as target). Mirrors
# brnaseq_pipeline.sh exactly: plain `salmon index`, no -k/--keepDuplicates.
# Sized generously (32 cpus / 96G) since the GRCm39.116 cDNA is much larger
# than typical (963M vs ~240M decompressed for most Ensembl releases).

set -euo pipefail

SCRATCH_DIR="/scratch/Lake/$USER/hfd_organoid_splicing"
REF_DIR="$SCRATCH_DIR/data/reference"
INDEX_DIR="$SCRATCH_DIR/data/salmon_index"

GENOME_FA="$REF_DIR/Mus_musculus.GRCm39.dna.primary_assembly.fa"
CDNA_FA="$REF_DIR/Mus_musculus.GRCm39.cdna.all.fa"
DECOYS="$REF_DIR/decoys.txt"
GENTROME="$REF_DIR/gentrome.fa"

echo "=== Building decoys.txt from genome sequence names ==="
grep '^>' "$GENOME_FA" | cut -d ' ' -f1 | sed 's/>//' > "$DECOYS"
wc -l "$DECOYS"

echo "=== Building gentrome.fa (cDNA + genome, cDNA first) ==="
cat "$CDNA_FA" "$GENOME_FA" > "$GENTROME"
ls -lh "$GENTROME"

echo "=== Building salmon index ==="
salmon index \
    -t "$GENTROME" \
    -d "$DECOYS" \
    -i "$INDEX_DIR" \
    -p "$SLURM_CPUS_PER_TASK"

echo ""
echo "=== Done. Salmon index: $INDEX_DIR ==="
ls -lh "$INDEX_DIR"
