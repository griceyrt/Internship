#!/bin/bash
#SBATCH --job-name=mm39_reference_download
#SBATCH --partition=Lake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=%j_mm39_reference_download.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biotools

# =============================================================================
# Downloads the mm39 (GRCm39) genome fasta, GTF, and stock cDNA
# (transcriptome) fasta from Ensembl.
#
# Uses Ensembl's "current_release" URL pattern rather than a hardcoded
# release number, so this keeps working as Ensembl's release ticks up.
# Genome build is primary_assembly (not toplevel) -- excludes
# haplotype/patch scaffolds, the standard choice for alignment.
# =============================================================================

BASE_DIR="$HOME/fmr1_polya_splicing"
REF_DIR="$BASE_DIR/data/reference"
mkdir -p "$REF_DIR"
cd "$REF_DIR"

echo "=== Genome fasta (primary_assembly) ==="
wget -N "https://ftp.ensembl.org/pub/current_fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"

echo "=== GTF annotation ==="
# Ensembl doesn't expose a browsable "current_gtf" directory, so query the
# REST API for the current release number and construct the filename
# directly (always Mus_musculus.GRCm39.<release>.gtf.gz).
ENSEMBL_RELEASE=$(curl -s "https://rest.ensembl.org/info/software?content-type=application/json" \
    | grep -oE '"release":[0-9]+' | grep -oE '[0-9]+')

if [ -z "$ENSEMBL_RELEASE" ]; then
    echo "ERROR: could not determine Ensembl release number from the REST API."
    echo "Check https://rest.ensembl.org/info/software?content-type=application/json manually."
else
    echo "Ensembl release: $ENSEMBL_RELEASE"
    wget -N "https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gtf/mus_musculus/Mus_musculus.GRCm39.${ENSEMBL_RELEASE}.gtf.gz"
fi

echo "=== Stock cDNA (transcriptome) fasta ==="
wget -N "https://ftp.ensembl.org/pub/current_fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz"

echo ""
echo "=== Verifying downloads ==="
ls -lh "$REF_DIR"
for f in "$REF_DIR"/*.gz; do
    echo -n "$f : "
    gunzip -t "$f" && echo "OK" || echo "CORRUPT -- re-download"
done

echo ""
echo "=== Done. Files in $REF_DIR ==="
echo "NOTE: verify the GTF filename matched a single file (there should be"
echo "exactly one Mus_musculus.GRCm39.<release>.gtf.gz -- if wget pulled"
echo "anything else via the recursive listing, delete the extras)."
