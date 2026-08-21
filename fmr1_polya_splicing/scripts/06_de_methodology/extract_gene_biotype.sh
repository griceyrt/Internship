#!/bin/bash
#SBATCH --job-name=extract_gene_biotype
#SBATCH --partition=Lake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=%j_extract_gene_biotype.log

# Extracts gene_id + gene_biotype from the same mm39 Ensembl GTF already
# used throughout this pipeline (SUPPA, DESeq2 gene models), for Panel A
# (differentially expressed, protein-coding vs non-coding split). Small,
# one-time, reproducible extraction -- not a manual lookup.

BASE_DIR="$HOME/fmr1_polya_splicing"
GTF="$BASE_DIR/data/reference/Mus_musculus.GRCm39.116.gtf"
OUT="$BASE_DIR/data/reference/gene_biotype_mm39.tsv"

if [ ! -f "$GTF" ]; then
    echo "ERROR: GTF not found at $GTF"
    echo "Check data/reference/ for the actual decompressed GTF filename and update this script."
    exit 1
fi

echo -e "gene_id\tgene_biotype" > "$OUT"
awk -F'\t' '$3=="gene"' "$GTF" | \
    sed -E 's/.*gene_id "([^"]+)".*gene_biotype "([^"]+)".*/\1\t\2/' \
    >> "$OUT"

echo "Saved: $OUT"
wc -l "$OUT"
echo "Biotype breakdown:"
tail -n +2 "$OUT" | cut -f2 | sort | uniq -c | sort -rn | head -10
