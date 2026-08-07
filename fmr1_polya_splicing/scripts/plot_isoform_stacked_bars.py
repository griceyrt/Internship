#!/usr/bin/env python3
# =============================================================================
# Isoform usage stacked-bar grids, WT vs Fmr1 KO
# Author: Gricey
#
# For each expression threshold (10/20/50 min TPM), plots the top 15 genes
# by entropy shift as a grid of small stacked-bar panels, one panel per
# gene, WT bar next to KO bar. Each bar segment is one isoform, sized to
# its share of that gene's total expression (0-1 scale). Colors are a
# continuous blue gradient built fresh per gene from that gene's own
# isoform count (dark = most-dominant isoform, light = least-dominant),
# so every segment is a genuine shade of blue -- no separate "other"
# bucket, even for genes with many isoforms. Same isoform keeps the same
# color across its own WT/KO pair, but color identity isn't meant to be
# compared across different genes.
#
# Inputs (already produced by isoform_proportions_step1-4.R):
#   results/isoform_proportions/isoform_proportions_long.tsv
#   results/isoform_proportions/entropy_shift_ranked_min{10,20,50}TPM.tsv
#
# This only needs the small TSV files, not the cluster -- run it locally.
# Requires: pandas, matplotlib (pip install pandas matplotlib if needed)
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

# Gricey's full 5-color palette pick, ordered darkest to lightest by actual
# perceptual luminance (Prussian Blue was darkest, Pale Sky was lightest --
# the first version of this script only used 4 of the 5 and left Pale Sky
# unused). Used as anchors for a CONTINUOUS gradient rather than a fixed
# list, so a gene with 3 isoforms and a gene with 37 isoforms each get their
# own dark-to-light blue spread, no separate "other"/grey bucket needed.
BLUE_PALETTE_ANCHORS = ["#061A40", "#003559", "#0353A4", "#006DAA", "#B9D6F2"]
BLUE_CMAP = LinearSegmentedColormap.from_list("gricey_blues", BLUE_PALETTE_ANCHORS)


def isoform_colors(n_isoforms):
    """n evenly-spaced shades from dark (most-dominant isoform) to light
    (least-dominant), scaled to THIS gene's own isoform count -- every
    gene uses the full dark-to-light range regardless of how many isoforms
    it has."""
    if n_isoforms == 1:
        return [BLUE_CMAP(0.0)]
    return [BLUE_CMAP(t) for t in np.linspace(0.0, 1.0, n_isoforms)]


def plot_gene(ax, gene_id, gene_name):
    sub = long_df[long_df["gene_id"] == gene_id]

    # Stable isoform order (by combined WT+KO proportion, descending) so
    # each isoform gets the same color/stack-position in both bars.
    order = (sub.groupby("transcript_id")["proportion"].sum()
             .sort_values(ascending=False).index.tolist())
    colors = isoform_colors(len(order))
    color_map = dict(zip(order, colors))

    # Isoforms absent from a condition (proportion=0, no row in the long
    # table) are simply skipped when drawing that bar -- correct, since a
    # 0-height segment can't be shown anyway. But that means a color used
    # in one bar can be entirely missing from the other bar right next to
    # it (e.g. a dominant isoform in WT that drops to fully undetected in
    # KO) -- this is a real biological result, not a plotting artifact,
    # but a reader scanning left-to-right could easily misread "color
    # disappeared" as "we ran out of colors" rather than "this isoform's
    # usage genuinely went to zero". Labeling each bar with how many
    # isoforms are actually present in it (n=) makes that explicit instead
    # of relying on the reader to infer it from color continuity alone.
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
    fig.text(
        0.5, -0.02,
        "Y-axis: proportion of gene's total expression per isoform (0-1). "
        "Each color = one isoform, consistent within a gene's own WT/KO pair. "
        "\"n=\" below each bar is how many isoforms were actually detected in that "
        "condition -- a color present in one bar and missing from the other means "
        "that isoform's usage genuinely dropped to zero, not a plotting artifact.",
        ha="center", fontsize=8, style="italic", wrap=True,
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
