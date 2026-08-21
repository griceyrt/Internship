#!/bin/bash
#SBATCH --job-name=grcm39_116_reference_download
#SBATCH --partition=Lake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=%j_grcm39_116_reference_download.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biotools

# Author: Gricey
#
# Downloads the genome FASTA, cDNA FASTA, and GTF used to build the salmon
# index. Release is hardcoded to 116 (confirmed with Desmond by email) --
# do not switch to "current_release", must match his DE analysis exactly.
# Genome assembly is primary_assembly (not toplevel), matching his index.
# Data goes on $SCRATCH, not $HOME (shared lab-group quota).

set -euo pipefail

ENSEMBL_RELEASE=116
SCRATCH_DIR="/scratch/Lake/$USER/hfd_organoid_splicing"
REF_DIR="$SCRATCH_DIR/data/reference"
mkdir -p "$REF_DIR"
cd "$REF_DIR"

# remove any stale release-115 files (unversioned filenames could otherwise survive)
rm -f "$REF_DIR"/Mus_musculus.GRCm39.dna.primary_assembly.fa* \
      "$REF_DIR"/Mus_musculus.GRCm39.cdna.all.fa* \
      "$REF_DIR"/Mus_musculus.GRCm39.115.gtf*

echo "=== Genome FASTA (primary_assembly), release ${ENSEMBL_RELEASE} ==="
wget -N "https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"

echo "=== cDNA (transcriptome) FASTA, release ${ENSEMBL_RELEASE} ==="
wget -N "https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz"

echo "=== GTF annotation, release ${ENSEMBL_RELEASE} ==="
wget -N "https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gtf/mus_musculus/Mus_musculus.GRCm39.${ENSEMBL_RELEASE}.gtf.gz"

echo ""
echo "=== Verifying downloads ==="
ls -lh "$REF_DIR"
for f in "$REF_DIR"/*.gz; do
    echo -n "$f : "
    gunzip -t "$f" && echo "OK" || echo "CORRUPT -- re-download"
done

echo ""
echo "=== Decompressing (salmon/STAR need unzipped FASTA/GTF) ==="
gunzip -k "$REF_DIR/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
gunzip -k "$REF_DIR/Mus_musculus.GRCm39.cdna.all.fa.gz"
gunzip -k "$REF_DIR/Mus_musculus.GRCm39.${ENSEMBL_RELEASE}.gtf.gz"

echo ""
echo "=== Done. Reference files in $REF_DIR ==="
