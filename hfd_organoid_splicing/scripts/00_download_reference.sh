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

# =============================================================================
# Phase 0 — reference download (HFD organoid splicing project)
# Author: Gricey
# Description: Pulls the genome FASTA, cDNA FASTA, and GTF that Desmond used
#              to build his salmon index (per brnaseq_pipeline.sh, shared by
#              Kiran 2026-08-03), so our splicing results stay comparable to
#              his DE results.
#
# RELEASE CORRECTED 2026-08-04: originally downloaded as release 115, going
# by the (as it turns out, unused/stale) "115" hardcoded into Desmond's
# script for a GTF he says he never actually used. Desmond confirmed by
# email (2026-08-04) that he does NOT use a GTF at all -- salmon quantifies
# directly against the cDNA transcriptome, no gene models needed for that.
# We need a GTF ourselves only because SUPPA2 requires one to define splicing
# events, which his pipeline never needed. He also confirmed his real cDNA
# download was from Ensembl release-116 (June/July), not 115 -- we can't
# verify this ourselves since his original reference files aren't in the
# shared drive, so we're going with his stated release. Downloading genome +
# cDNA + GTF all fresh from 116 (not just the GTF) for full consistency,
# since transcript models can shift slightly between releases even for a
# stable species like mouse.
#
# IMPORTANT: release number is HARDCODED, not "current_release" -- unlike
# project 3 (fmr1_polya_splicing), where any current mm39 build was fine,
# here we need the exact same annotation version Desmond indexed against.
# Do not change this without checking with Desmond first, since Ensembl
# releases can add/drop/rename transcripts between versions.
#
# Genome: primary_assembly (not toplevel) -- matches Desmond's file name
# (Mus_musculus.GRCm39.dna.primary_assembly.fa) and is the standard choice
# to avoid duplicate haplotype/patch scaffolds.
#
# STORAGE: reference files go on $SCRATCH, not $HOME. PSMN's $HOME quota is
# shared across the whole lab group (not per-user) and is not meant for
# large data -- $SCRATCH has no quota and is the documented place for large
# input/output files (see PSMN docs: "Where to store files?"). Only scripts
# live in $HOME.
#
# PREREQUISITE: run this directly on the cluster (not locally) -- it just
# downloads from Ensembl, no local data needed.
# =============================================================================

set -euo pipefail

ENSEMBL_RELEASE=116
SCRATCH_DIR="/scratch/Lake/$USER/hfd_organoid_splicing"
REF_DIR="$SCRATCH_DIR/data/reference"
mkdir -p "$REF_DIR"
cd "$REF_DIR"

echo "=== Removing any stale release-115 files (genome/cDNA filenames don't ==="
echo "=== embed a release number, so a stale copy could otherwise silently ==="
echo "=== survive and get reused instead of the real 116 download)         ==="
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
echo "Expected files:"
echo "  Mus_musculus.GRCm39.dna.primary_assembly.fa"
echo "  Mus_musculus.GRCm39.cdna.all.fa"
echo "  Mus_musculus.GRCm39.${ENSEMBL_RELEASE}.gtf"
