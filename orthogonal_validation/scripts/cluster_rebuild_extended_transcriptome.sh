#!/bin/bash
#SBATCH --job-name=rebuild_ext_transcriptome
#SBATCH --output=%x_%j.log
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#
# Rebuilds the Salmon transcriptome/index using the 3'UTR-extended GTF
# (transcriptome_productivity_extended.gtf -- extended 1000bp toward 5',
# 100bp toward 3', matching Asher's own MARS-seq processing approach).
#
# Steps:
#   1. Reconstruct the standalone genome FASTA from gentrome_ext.fa's decoy
#      sequences (no separate genome FASTA is staged on the cluster --
#      it's baked into gentrome_ext.fa, so we pull it back out instead of
#      re-transferring 2.7GB from the Mac).
#   2. gffread: extract transcript sequences from the extended GTF + genome.
#   3. Rebuild gentrome (extended transcriptome + same genome decoy).
#   4. Build new Salmon index.
#
# USAGE: sbatch cluster_rebuild_extended_transcriptome.sh
# Run from: ~/orthogonal_validation/
# Requires: biotools conda env active (gffread, samtools, salmon)

set -euo pipefail

DATA_DIR="$HOME/orthogonal_validation/data"
GENTROME_ORIG="$DATA_DIR/gentrome_ext.fa"
DECOYS="$DATA_DIR/decoys.txt"
GTF_EXT="$DATA_DIR/transcriptome_productivity_extended.gtf"

GENOME_ONLY="$DATA_DIR/genome_only.fa"
TRANSCRIPTOME_EXT_NEW="$DATA_DIR/transcriptome_ext_extended.fa"
GENTROME_NEW="$DATA_DIR/gentrome_ext_extended.fa"
SALMON_INDEX_NEW="$DATA_DIR/transcriptome_ext_extended_index"

echo "[1/4] Reconstructing standalone genome FASTA from gentrome_ext.fa decoys..."
if [ ! -f "$GENOME_ONLY" ]; then
    samtools faidx "$GENTROME_ORIG"   # builds gentrome_ext.fa.fai if missing
    samtools faidx "$GENTROME_ORIG" $(cat "$DECOYS") > "$GENOME_ONLY"
    samtools faidx "$GENOME_ONLY"     # index the extracted genome for gffread
else
    echo "  $GENOME_ONLY already exists, skipping."
fi

echo "[2/4] Running gffread to extract transcript sequences from extended GTF..."
gffread -w "$TRANSCRIPTOME_EXT_NEW" -g "$GENOME_ONLY" "$GTF_EXT"
echo "  Done: $TRANSCRIPTOME_EXT_NEW"
grep -c "^>" "$TRANSCRIPTOME_EXT_NEW"

echo "[3/4] Building new gentrome (extended transcriptome + genome decoy)..."
cat "$TRANSCRIPTOME_EXT_NEW" "$GENOME_ONLY" > "$GENTROME_NEW"
echo "  Done: $GENTROME_NEW"

echo "[4/4] Building new Salmon index..."
salmon index -t "$GENTROME_NEW" \
    -i "$SALMON_INDEX_NEW" \
    --decoys "$DECOYS" \
    -k 31 \
    --keepDuplicates \
    -p "$SLURM_CPUS_PER_TASK"

echo ""
echo "All done. New index at: $SALMON_INDEX_NEW"
echo "Next: re-quantify with the existing pipeline script, pointing at this new index:"
echo "  sbatch --wrap=\"bash scripts/get_sra_dump_fastq_run_salmon.sh SRP194523 $SALMON_INDEX_NEW\""
