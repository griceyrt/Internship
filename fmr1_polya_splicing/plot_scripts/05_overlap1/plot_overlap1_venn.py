#!/usr/bin/env python3
# =============================================================================
# Overlap #1 Venn diagrams: DESeq2-significant genes vs SUPPA dPSI-significant
# genes, one panel per event type (SE, A3, A5, MX, RI, AF, AL).
# Author: Gricey
#
# 7 separate panels rather than one pooled Venn, matching this project's
# established convention of testing each SUPPA event type independently
# (never pooled into one "any event type" bucket) -- same design already
# used for the SUPPA significance counts and Overlap #1 itself.
#
# Inputs (both produced by overlap1_de_vs_dpsi.R -- nothing computed here
# that isn't already a saved pipeline output, for reproducibility):
#   results/overlap1/overlap1_summary_by_event_type.tsv
#     -> per event type: n_genes_significant (SUPPA circle size),
#        n_genes_overlap_with_DE (intersection size)
#   results/overlap1/overlap1_deseq2_total_sig_genes.tsv
#     -> single number, DESeq2 circle size (same across all 7 panels)
#
# Requires: pandas, matplotlib, matplotlib-venn
#   pip install pandas matplotlib matplotlib-venn
# =============================================================================

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"]

BASE = "/Users/gricey/Desktop/Internship/fmr1_polya_splicing/results/overlap1"
OUT_DIR = BASE

summary = pd.read_csv(f"{BASE}/overlap1_summary_by_event_type.tsv", sep="\t")
deseq_total_df = pd.read_csv(f"{BASE}/overlap1_deseq2_total_sig_genes.tsv", sep="\t")
n_deseq_sig = int(deseq_total_df.loc[deseq_total_df["metric"] == "n_deseq2_significant_genes_padj0.05", "value"].iloc[0])

print(f"DESeq2-significant genes (padj<0.05), same for every panel: {n_deseq_sig}")

EVENT_TYPES = ["SE", "A3", "A5", "MX", "RI", "AF", "AL"]

fig, axes = plt.subplots(2, 4, figsize=(16, 8))
axes = axes.flatten()

for ax, event in zip(axes, EVENT_TYPES):
    row = summary[summary["event_type"] == event]
    if row.empty:
        ax.axis("off")
        ax.set_title(f"{event}\n(no data)", fontsize=11)
        continue

    n_suppa_sig = int(row["n_genes_significant"].iloc[0])
    n_overlap = int(row["n_genes_overlap_with_DE"].iloc[0])

    # matplotlib_venn subset sizes are (only-set-A, only-set-B, both)
    deseq_only = n_deseq_sig - n_overlap
    suppa_only = n_suppa_sig - n_overlap

    venn2(subsets=(deseq_only, suppa_only, n_overlap),
          set_labels=("DESeq2", f"SUPPA {event}"),
          ax=ax)
    ax.set_title(f"{event}\n(DESeq2 n={n_deseq_sig}, SUPPA n={n_suppa_sig})", fontsize=10)

# 8th subplot unused (7 event types in a 2x4 grid) -- hide it
axes[7].axis("off")

fig.suptitle(
    "Overlap #1: DESeq2-significant genes ∩ SUPPA dPSI-significant genes, by event type",
    fontsize=14, fontweight="bold", y=1.02,
)
fig.tight_layout()

png_path = f"{OUT_DIR}/overlap1_venn_by_event_type.png"
svg_path = f"{OUT_DIR}/overlap1_venn_by_event_type.svg"
fig.savefig(png_path, dpi=200, bbox_inches="tight")
fig.savefig(svg_path, bbox_inches="tight")
plt.close(fig)
print("Saved:", png_path)
print("Saved:", svg_path)
