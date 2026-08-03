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
# Download script for Phase 0 reference materials (mm39 Ensembl)
# Author: Gricey
# Description: Pulls the mm39 (GRCm39) genome fasta, GTF, and stock cDNA
#              (transcriptome) fasta directly from Ensembl -- decided
#              2026-07-29 to direct-download the transcriptome rather than
#              build it with gffread, for a cleaner/simpler path (equivalent
#              result, avoids any genome/GTF release-version mismatch risk).
#
# Uses Ensembl's "current_release" URL pattern (not a hardcoded release
# number) so this keeps working even as Ensembl's release number ticks up.
# Genome: primary_assembly (not toplevel) -- excludes haplotype/patch
# scaffolds, the standard choice for alignment to avoid duplicate regions.
# =============================================================================

BASE_DIR="$HOME/fmr1_polya_splicing"
REF_DIR="$BASE_DIR/data/reference"
mkdir -p "$REF_DIR"
cd "$REF_DIR"

echo "=== Genome fasta (primary_assembly) ==="
wget -N "https://ftp.ensembl.org/pub/current_fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"

echo "=== GTF annotation ==="
# NOTE (fixed 2026-08-03, take 2): both the recursive-wget approach and a
# plain curl-scrape of the "current_gtf" directory page failed -- that
# directory path returns a genuine 404, it doesn't exist as a listable
# path on this server (unlike current_fasta, which does). Rather than
# guess at directory structure a third time, query Ensembl's actual REST
# API for the current release number (a real, documented endpoint, not a
# scraped page), then construct the exact GTF filename directly --
# Ensembl's GTF files are always named Mus_musculus.GRCm39.<release>.gtf.gz.
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
