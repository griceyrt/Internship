"""
Computes the percentage breakdown of htseq-count's QC counters
(assigned / ambiguous / not_unique / no_feature / too_low_aQual / not_aligned)
per sample, from the *_htseq_qc_counters.tsv file produced by
combine_htseq_counts.py.

These counters matter for this project specifically because a high
"ambiguous" or "not_unique" rate means a large fraction of reads are being
discarded rather than counted -- directly relevant to whether STAR+htseq-count
is giving a "clean" gene-level count, and whether the 3'UTR/5' boundary
extension (which makes neighboring genes' regions more likely to overlap)
is inflating the ambiguous rate specifically.

Usage:
    python3 analyze_htseq_qc.py <qc_counters.tsv> <gene_counts.tsv>

Example:
    python3 analyze_htseq_qc.py \
        results/SRP194523_star_htseq_2026-07-17/combined_htseq_qc_counters.tsv \
        results/SRP194523_star_htseq_2026-07-17/combined_gene_counts.tsv

Run from: orthogonal_validation/
"""
import sys
import pandas as pd


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)

    qc_path, counts_path = sys.argv[1], sys.argv[2]

    qc = pd.read_csv(qc_path, sep="\t", index_col=0)

    counts = pd.read_csv(counts_path, sep="\t", index_col=0)
    sample_cols = [c for c in counts.columns if c != "gene_name"]
    assigned = counts[sample_cols].sum()

    full = pd.concat([assigned.rename("assigned"), qc.T], axis=1)
    totals = full.sum(axis=1)

    print("=== htseq-count QC breakdown per sample ===\n")
    for sample in full.index:
        total = totals[sample]
        print(f"{sample}  (total processed: {total:,.0f})")
        for col in full.columns:
            val = full.loc[sample, col]
            print(f"    {col:20s}: {val:>12,.0f}  ({100*val/total:5.1f}%)")
        print()

    print("=== Average across all samples ===")
    avg_pct = (full.div(totals, axis=0) * 100).mean()
    for col, pct in avg_pct.items():
        print(f"    {col:20s}: {pct:5.1f}%")


if __name__ == "__main__":
    main()
