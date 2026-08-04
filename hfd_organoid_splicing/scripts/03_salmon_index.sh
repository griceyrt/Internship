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

# =============================================================================
# Phase 5 — salmon decoy-aware index (HFD organoid splicing project)
# Author: Gricey
# Description: Builds the decoy-aware salmon index against the corrected
#              release-116 reference (see notes/pipeline_plan.txt and
#              Desmond's 2026-08-04 emails for why 116, not 115). Steps and
#              flags mirror brnaseq_pipeline.sh exactly (decoy list from
#              genome headers, cDNA+genome concatenated into one "gentrome"
#              fasta, plain `salmon index` with no extra flags -- Desmond's
#              script didn't use -k/--keepDuplicates, so neither do we, for
#              consistency with his build).
#
# Mem/cpus/time sized generously (32 cpus / 96G / 12h -- a full Lake node
# has 32 cores and 192GB RAM, so this uses all the cores and half the RAM)
# since release 116's cDNA fasta turned out much larger than 115's (963M vs
# 238M decompressed) -- see memory notes on that size jump. Adjust down in
# future re-runs once we know the real indexing time from this first run.
#
# STORAGE: everything on $SCRATCH, only this script lives in $HOME.
# =============================================================================

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
