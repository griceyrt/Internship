#!/usr/bin/env python3
# =============================================================================
# Isoform usage stacked-bar grids, WT vs Fmr1 KO
#
# For each expression threshold (10/20/50 min TPM), plots the top 15 genes
# by entropy shift as a grid of stacked-bar panels, one panel per gene, WT
# bar next to KO bar. Each bar segment is one isoform, sized to its share
# of that gene's total expression (0-1 scale). Colors use a continuous
# Set3-based gradient built per gene from its own isoform count, so colors
# never repeat within a gene's bars regardless of isoform count. Same
# isoform keeps the same color across its own WT/KO pair; color identity
# is not meant to be compared across different genes.
#
# Inputs (produced by isoform_proportions_step1-4.R):
#   results/isoform_proportions/isoform_proportions_long.tsv
#   results/isoform_proportions/entropy_shift_ranked_min{10,20,50}TPM.tsv
#
# Run locally (only needs the small TSV files, not the cluster).
# Requires: pandas, matplotlib
# =============================================================================

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# Arial throughout; falls back to a common substitute if Arial isn't
# installed on this machine's matplotlib font cache (e.g. Helvetica/Liberation
# Sans render nearly identically)
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"]

BASE = "/Users/gricey/Desktop/Internship/fmr1_polya_splicing/results/isoform_proportions"
OUT_DIR = BASE  # PNG + SVG land alongside the source tables; change if you want them elsewhere

long_df = pd.read_csv(f"{BASE}/isoform_proportions_long.tsv", sep="\t")

# Set3 (ColorBrewer) is natively a fixed 12-color palette, but some genes
# have far more isoforms than that (e.g. up to 37) -- using it literally
# would force color repeats within a single gene's bar. Instead, Set3's 12
# official hex values are used as anchors for a continuous colormap,
# sampled at N evenly-spaced points for whatever N isoforms a gene has, so
# colors never repeat within a gene regardless of isoform count.
SET3_PALETTE_ANCHORS = ["#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
                         "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
                         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F"]
ISOFORM_CMAP = LinearSegmentedColormap.from_list("gricey_set3", SET3_PALETTE_ANCHORS)


def isoform_colors(n_isoforms):
    """n evenly-spaced hues across the full rainbow, scaled to THIS gene's
    own isoform count -- every gene uses the full spectrum regardless of
    how many isoforms it has, so colors never repeat within one gene."""
    if n_isoforms == 1:
        return [ISOFORM_CMAP(0.0)]
    return [ISOFORM_CMAP(t) for t in np.linspace(0.0, 1.0, n_isoforms)]


def plot_gene(ax, gene_id, gene_name):
    sub = long_df[long_df["gene_id"] == gene_id]

    # Stable isoform order (by combined WT+KO proportion, descending) so
    # each isoform gets the same color/stack-position in both bars.
    order = (sub.groupby("transcript_id")["proportion"].sum()
             .sort_values(ascending=False).index.tolist())
    colors = isoform_colors(len(order))
    color_map = dict(zip(order, colors))

    # Isoforms absent from a condition (proportion=0) are skipped when
    # drawing that bar -- a color present in one bar and missing from the
    # other reflects a real drop to zero usage, not a plotting artifact.
    # The "n=" label on each bar makes the isoform count explicit.
    n_present = {}
    for x, cond in enumerate(["WT", "KO"]):
        cond_data = sub[sub["condition"] == cond].set_index("transcript_id")
        n_present[cond] = len(cond_data)
        bottom = 0
        for tx in order:
            if tx in cond_data.index:
                val = cond_data.loc[tx, "proportion"]
                ax.bar(x, val, bottom=bottom, width=0.6, color=color_map[tx],
                       edgecolor="white", linewidth=0.5)
                bottom += val

    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"WT\n(n={n_present['WT']})", f"KO\n(n={n_present['KO']})"], fontsize=7)
    ax.set_ylim(0, 1)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_yticklabels(["0", "0.25", "0.5", "0.75", "1"], fontsize=7)
    ax.set_title(gene_name, fontsize=10, fontweight="bold")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


def make_grid(threshold, n_genes=15, ncols=5):
    df = pd.read_csv(f"{BASE}/entropy_shift_ranked_min{threshold}TPM.tsv", sep="\t")
    df = df.sort_values("abs_entropy_diff", ascending=False).head(n_genes)

    nrows = -(-n_genes // ncols)  # ceil division
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 2.6, nrows * 2.8))
    axes = axes.flatten()

    for ax, (_, row) in zip(axes, df.iterrows()):
        plot_gene(ax, row["gene_id"], row["gene_name"])
    for ax in axes[len(df):]:
        ax.axis("off")

    fig.suptitle(
        f"Isoform usage proportions, WT vs Fmr1 KO (min TPM > {threshold}, top {n_genes} by entropy shift)",
        fontsize=13, fontweight="bold", y=1.02,
    )
    fig.tight_layout()

    png_path = f"{OUT_DIR}/stacked_bars_min{threshold}TPM.png"
    svg_path = f"{OUT_DIR}/stacked_bars_min{threshold}TPM.svg"
    fig.savefig(png_path, dpi=200, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")  # vector version -- no dpi needed, scales cleanly for print/editing
    plt.close(fig)
    print("Saved:", png_path)
    print("Saved:", svg_path)


if __name__ == "__main__":
    for t in (10, 20, 50):
        make_grid(t)
