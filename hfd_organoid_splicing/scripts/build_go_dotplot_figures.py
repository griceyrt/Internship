#!/usr/bin/env python3
"""
GO enrichment dot plots (HFD organoid splicing project)
Author: Gricey

For each of the 5 GO result sets (4 timepoints + union), plots the top 20
significant GO Biological Process terms: term on the y-axis,
-log10(adjusted p-value) on the x-axis, dot size = gene count, dot color =
adjusted p-value. Inspired by Fig 3G of Chikhaoui_Mamgain_MS.pdf (same lab,
PERIOD paper); size/color use Count and p.adjust instead of their "delta
entropy" metric, following the standard clusterProfiler dotplot() convention.

Run locally, not on the cluster -- lightweight plotting.
USAGE: python3 scripts/build_go_dotplot_figures.py
Run from: hfd_organoid_splicing/
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")

BASE = "."
GO_DIR = f"{BASE}/results/suppa/go_enrichment"
FIGURES_DIR = f"{BASE}/figures"

COMPARISONS = [
    ("W8_WHFD_vs_WT_SD", "W8_WHFD vs WT_SD"),
    ("W18_WHFD_vs_WT_SD", "W18_WHFD vs WT_SD"),
    ("W24_WHFD_vs_WT_SD", "W24_WHFD vs WT_SD"),
    ("W42_WHFD_vs_WT_SD", "W42_WHFD vs WT_SD"),
    ("union_all_comparisons", "Union across all timepoints"),
]

TOP_N_TERMS = 20

for file_key, title in COMPARISONS:
    csv_path = f"{GO_DIR}/go_enrichment_{file_key}.csv"
    df = pd.read_csv(csv_path)
    df = df[df["p.adjust"] < 0.05].sort_values("p.adjust").head(TOP_N_TERMS)

    if df.empty:
        print(f"Skipping {file_key} — no significant terms")
        continue

    # Sort so the most significant term is at the TOP of the plot
    df = df.sort_values("p.adjust", ascending=False)
    df["neg_log10_padj"] = -np.log10(df["p.adjust"])

    fig, ax = plt.subplots(figsize=(7, max(5, 0.35 * len(df))))
    scatter = ax.scatter(
        df["neg_log10_padj"], df["Description"],
        s=df["Count"] * 8, c=df["p.adjust"],
        cmap="viridis_r", edgecolors="black", linewidths=0.6, alpha=0.9,
    )

    ax.set_xlabel("-log10(adjusted p-value)", fontsize=11)
    ax.set_title(f"GO Biological Process — {title}", fontsize=13, fontweight="bold")
    ax.tick_params(axis="y", labelsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.margins(x=0.18)  # extra room so the largest bubbles don't clip at the edges

    # Colorbar (p.adjust)
    cbar = fig.colorbar(scatter, ax=ax, pad=0.02, fraction=0.05)
    cbar.set_label("adjusted p-value", fontsize=9)

    # Size legend (Count) -- proxy handles at a few representative sizes
    count_min, count_max = df["Count"].min(), df["Count"].max()
    legend_counts = sorted(set([count_min, int((count_min + count_max) / 2), count_max]))
    size_handles = [
        plt.scatter([], [], s=c * 8, color="grey", edgecolors="black", linewidths=0.6, label=str(c))
        for c in legend_counts
    ]
    ax.legend(handles=size_handles, title="no. of genes", loc="lower right",
              bbox_to_anchor=(1.32, 0), fontsize=8, title_fontsize=9,
              frameon=False, labelspacing=1.2)

    plt.tight_layout()
    png_path = f"{FIGURES_DIR}/go_dotplot_{file_key}.png"
    svg_path = f"{FIGURES_DIR}/go_dotplot_{file_key}.svg"
    plt.savefig(png_path, dpi=150, bbox_inches="tight")
    plt.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {png_path}")
    print(f"Saved: {svg_path}")

print("\nDone.")
