"""
Combines the 8 per-sample htseq-count output files into a single gene x
sample count matrix, ready for DESeq2. htseq-count output (with
--additional-attr=gene_name) has gene_id/gene_name/count columns plus special
QC summary rows starting with "__" (no_feature, ambiguous, etc.), split out
separately here since they're useful multi-mapping diagnostics.
Author: Gricey

Usage:
    python3 combine_htseq_counts.py <results_dir> <output_prefix>

Example:
    python3 combine_htseq_counts.py \
        results/GSE130613/SRP194523_star_htseq_2026-07-16 \
        results/GSE130613/SRP194523_star_htseq_2026-07-16/combined

Run from: orthogonal_validation/
"""
import sys
import os
import pandas as pd

SAMPLES_WT = ["SRR9002567", "SRR9002568", "SRR9002569", "SRR9002570"]
SAMPLES_KO = ["SRR9002583", "SRR9002584", "SRR9002585", "SRR9002586"]
ALL_SAMPLES = SAMPLES_WT + SAMPLES_KO

SPECIAL_ROWS = {"__no_feature", "__ambiguous", "__too_low_aQual",
                "__not_aligned", "__alignment_not_unique"}


def load_one(results_dir, srr):
    path = os.path.join(results_dir, srr, f"{srr}_counts.txt")
    df = pd.read_csv(path, sep="\t", header=None,
                      names=["gene_id", "gene_name", srr])
    special = df[df["gene_id"].isin(SPECIAL_ROWS)].set_index("gene_id")[srr]
    genes = df[~df["gene_id"].isin(SPECIAL_ROWS)]
    return genes, special


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    results_dir, out_prefix = sys.argv[1], sys.argv[2]

    gene_names = None
    counts = {}
    specials = {}

    for srr in ALL_SAMPLES:
        genes, special = load_one(results_dir, srr)
        genes = genes.set_index("gene_id")
        if gene_names is None:
            gene_names = genes["gene_name"]
        counts[srr] = genes[srr]
        specials[srr] = special
        print(f"{srr}: {genes[srr].sum():,} reads assigned to genes")

    count_matrix = pd.DataFrame(counts)
    count_matrix.insert(0, "gene_name", gene_names)
    count_matrix.to_csv(f"{out_prefix}_gene_counts.tsv", sep="\t")
    print(f"\nSaved combined count matrix: {out_prefix}_gene_counts.tsv "
          f"({count_matrix.shape[0]} genes x {len(ALL_SAMPLES)} samples)")

    special_df = pd.DataFrame(specials)
    special_df.to_csv(f"{out_prefix}_htseq_qc_counters.tsv", sep="\t")
    print(f"Saved QC counters (ambiguous/multimapped/unassigned reads per sample): "
          f"{out_prefix}_htseq_qc_counters.tsv")
    print()
    print("=== QC counters summary ===")
    print(special_df)


if __name__ == "__main__":
    main()
