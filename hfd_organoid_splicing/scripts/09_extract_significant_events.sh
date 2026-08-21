#!/bin/bash
#SBATCH --job-name=hfd_organoid_extract_sig
#SBATCH --partition=Lake-short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=%j_hfd_organoid_extract_sig.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biotools

# Author: Gricey
#
# Filters all 28 .dpsi files (4 timepoints x 7 event types) to significant
# events (|dPSI| > 0.1, p < 0.05 -- same cutoffs as the lab's PERIOD paper,
# Chikhaoui/Mamgain), and extracts gene lists for GO enrichment (10_/11_).
# event_id format: "<gene_id>;<event_type>:<chr>:<coords>:<strand>".
# Pure Python standard library, no pandas.

set -euo pipefail

SCRATCH_DIR="/scratch/Lake/$USER/hfd_organoid_splicing"
DIFF_DIR="$SCRATCH_DIR/results/suppa/diffsplice"
OUT_DIR="$SCRATCH_DIR/results/suppa/significant"
mkdir -p "$OUT_DIR"

PVAL_CUTOFF=0.05
DPSI_CUTOFF=0.1

python3 - "$DIFF_DIR" "$OUT_DIR" "$PVAL_CUTOFF" "$DPSI_CUTOFF" << 'PYEOF'
import sys, glob, os
from collections import defaultdict

diff_dir, out_dir, pval_cutoff, dpsi_cutoff = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4])

all_rows = []  # (comparison, event_type, gene_id, event_id, dpsi, pval)
sig_genes_by_comparison = defaultdict(set)

dpsi_files = sorted(glob.glob(f"{diff_dir}/*.dpsi"))
print(f"Found {len(dpsi_files)} .dpsi files")

for path in dpsi_files:
    fname = os.path.basename(path)[:-5]  # strip .dpsi
    # filename format: <comparison>_events_<TYPE>_strict, e.g.
    # W8_WHFD_vs_WT_SD_events_SE_strict
    comparison, event_type = fname.split("_events_")
    event_type = event_type.replace("_strict", "")

    n_tested = 0
    n_sig = 0
    with open(path) as f:
        f.readline()  # header
        for line in f:
            fields = line.rstrip("\n").split("\t")
            event_id, dpsi_str, pval_str = fields[0], fields[1], fields[2]
            if dpsi_str == "nan" or pval_str == "nan":
                continue
            n_tested += 1
            dpsi, pval = float(dpsi_str), float(pval_str)
            if abs(dpsi) > dpsi_cutoff and pval < pval_cutoff:
                n_sig += 1
                gene_id = event_id.split(";")[0]
                all_rows.append((comparison, event_type, gene_id, event_id, dpsi, pval))
                sig_genes_by_comparison[comparison].add(gene_id)

    print(f"  {comparison} / {event_type}: {n_sig} significant / {n_tested} tested")

# Master table: every significant event, all comparisons, all event types
master_path = f"{out_dir}/significant_events.tsv"
with open(master_path, "w") as out:
    out.write("comparison\tevent_type\tgene_id\tevent_id\tdPSI\tp_value\n")
    for row in all_rows:
        out.write("\t".join(str(x) for x in row) + "\n")
print(f"\nWrote master table: {master_path} ({len(all_rows)} significant event rows)")

# Per-comparison gene lists (deduplicated), for GO enrichment input
for comparison, genes in sig_genes_by_comparison.items():
    gene_list_path = f"{out_dir}/significant_genes_{comparison}.txt"
    with open(gene_list_path, "w") as out:
        for g in sorted(genes):
            out.write(g + "\n")
    print(f"Wrote gene list: {gene_list_path} ({len(genes)} unique genes)")

# Union gene list across all 4 comparisons (any significant splicing change vs baseline)
union_genes = set()
for genes in sig_genes_by_comparison.values():
    union_genes |= genes
union_path = f"{out_dir}/significant_genes_union_all_comparisons.txt"
with open(union_path, "w") as out:
    for g in sorted(union_genes):
        out.write(g + "\n")
print(f"Wrote union gene list: {union_path} ({len(union_genes)} unique genes)")
PYEOF

echo ""
echo "=== Done. Output: $OUT_DIR ==="
ls -lh "$OUT_DIR"
