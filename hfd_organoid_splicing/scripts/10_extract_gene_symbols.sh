#!/bin/bash
#SBATCH --job-name=hfd_organoid_gene_symbols
#SBATCH --partition=Lake-short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=%j_hfd_organoid_gene_symbols.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biotools

# Author: Gricey
#
# Extracts gene_id -> gene_name from the GTF's own attributes (same source
# as org.Mm.eg.db) -- for labeling figures with symbols instead of Ensembl IDs.

set -euo pipefail

SCRATCH_DIR="/scratch/Lake/$USER/hfd_organoid_splicing"
GTF="$SCRATCH_DIR/data/reference/Mus_musculus.GRCm39.116.gtf"
OUT="$SCRATCH_DIR/data/reference/gene_id_to_symbol.tsv"

python3 - "$GTF" "$OUT" << 'PYEOF'
import re, sys

gtf_path, out_path = sys.argv[1], sys.argv[2]
seen = {}
gene_id_re = re.compile(r'gene_id "([^"]+)"')
gene_name_re = re.compile(r'gene_name "([^"]+)"')

with open(gtf_path) as f:
    for line in f:
        if line.startswith("#"):
            continue
        gid_m = gene_id_re.search(line)
        gname_m = gene_name_re.search(line)
        if gid_m and gname_m:
            seen[gid_m.group(1)] = gname_m.group(1)

with open(out_path, "w") as out:
    out.write("gene_id\tgene_name\n")
    for gid, gname in sorted(seen.items()):
        out.write(f"{gid}\t{gname}\n")

print(f"Wrote {len(seen)} gene_id -> gene_name mappings to {out_path}")
PYEOF

echo ""
echo "=== Done ==="
wc -l "$OUT"
head -5 "$OUT"
