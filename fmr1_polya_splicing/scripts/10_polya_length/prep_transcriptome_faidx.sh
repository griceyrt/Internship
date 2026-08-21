#!/bin/bash
#SBATCH --job-name=prep_transcriptome_faidx
#SBATCH --partition=Lake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=%j_prep_transcriptome_faidx.log

# One-time prep for Nanopolish: it needs a plain (not gzipped), samtools-
# faidx-indexed reference fasta -- the stock transcriptome we already have
# is gzipped and has never been indexed this way (Minimap2 reads .gz
# directly, so this was never needed until now).

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biotools

BASE_DIR="$HOME/fmr1_polya_splicing"
FASTA_GZ="$BASE_DIR/data/reference/Mus_musculus.GRCm39.cdna.all.fa.gz"
FASTA="$BASE_DIR/data/reference/Mus_musculus.GRCm39.cdna.all.fa"

if [ ! -f "$FASTA" ]; then
    echo "Decompressing transcriptome fasta..."
    gunzip -k "$FASTA_GZ"
else
    echo "$FASTA already exists, skipping decompression."
fi

echo "Running samtools faidx..."
samtools faidx "$FASTA"

echo "Done. Index file:"
ls -lh "${FASTA}.fai"
