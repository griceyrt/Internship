"""
Combines the 7 per-sample htseq-count output files for GSE171975 CT20 into a
single gene x sample count matrix, ready for DESeq2. Same logic as
combine_htseq_counts.py (GSE130613), different sample set/naming.
Author: Gricey

Usage:
    python3 GSE171975_combine_htseq_counts.py <results_dir> <output_prefix>

Example:
    python3 GSE171975_combine_htseq_counts.py \
        results/GSE171975/GSE171975_CT20_2026-08-17 \
        results/GSE171975/GSE171975_CT20_2026-08-17/combined

Run from: orthogonal_validation/
"""
import sys
import os
import pandas as pd

SAMPLES_WT = ["WT_A", "WT_B", "WT_C", "WT_D"]
SAMPLES_KO = ["PerDKO_A", "PerDKO_B", "PerDKO_C"]
ALL_SAMPLES = SAMPLES_WT + SAMPLES_KO

SPECIAL_ROWS = {"__no_feature", "__ambiguous", "__too_low_aQual",
                "__not_aligned", "__alignment_not_unique"}


def load_one(results_dir, label):
    path = os.path.join(results_dir, label, f"{label}_counts.txt")
    df = pd.read_csv(path, sep="\t", header=None,
                      names=["gene_id", "gene_name", label])
    special = df[df["gene_id"].isin(SPECIAL_ROWS)].set_index("gene_id")[label]
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

    for label in ALL_SAMPLES:
        genes, special = load_one(results_dir, label)
        genes = genes.set_index("gene_id")
        if gene_names is None:
            gene_names = genes["gene_name"]
        counts[label] = genes[label]
        specials[label] = special
        print(f"{label}: {genes[label].sum():,} reads assigned to genes")

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
