#!/usr/bin/env python3
"""
Summary figures: event-type bar chart + gene count over time
(HFD organoid splicing project)
Author: Gricey

Two figures: (1) mirrored bar chart of significant event counts per event
type, split up/down, one panel per timepoint (Fig 3B equivalent from
Chikhaoui_Mamgain_MS.pdf); (2) significant unique-gene count per timepoint
as a dot plot, highlighting the W18 peak.

Run locally, not on the cluster -- lightweight plotting.
USAGE: python3 scripts/build_summary_figures.py
Run from: hfd_organoid_splicing/
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")

BASE = "."
SIG_TABLE = f"{BASE}/results/suppa/significant/significant_events.tsv"
FIGURES_DIR = f"{BASE}/figures"

EVENT_TYPE_ORDER = ["SE", "MX", "RI", "A5", "A3", "AL", "AF"]
COMPARISON_ORDER = ["W8_WHFD", "W18_WHFD", "W24_WHFD", "W42_WHFD"]

sig_df = pd.read_csv(SIG_TABLE, sep="\t")
sig_df["direction"] = sig_df["dPSI"].apply(lambda x: "up" if x > 0 else "down")

# --- Figure 1: mirrored bar chart per comparison ---
fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharey=True)
fig.suptitle("Significant splicing events per type — up vs down, per timepoint",
             fontsize=15, fontweight="bold")

for ax, comparison in zip(axes.flat, COMPARISON_ORDER):
    comp_df = sig_df[sig_df["comparison"] == f"{comparison}_vs_WT_SD"]
    counts = comp_df.groupby(["event_type", "direction"]).size().unstack(fill_value=0)
    counts = counts.reindex(EVENT_TYPE_ORDER)

    x = range(len(EVENT_TYPE_ORDER))
    up_counts = counts.get("up", pd.Series(0, index=EVENT_TYPE_ORDER))
    down_counts = counts.get("down", pd.Series(0, index=EVENT_TYPE_ORDER))

    ax.bar(x, up_counts, color="mediumturquoise", label="Up", width=0.6)
    ax.bar(x, -down_counts, color="mediumpurple", label="Down", width=0.6)
    ax.axhline(0, color="black", linewidth=0.8)

    ax.set_xticks(list(x))
    ax.set_xticklabels(EVENT_TYPE_ORDER)
    ax.set_title(f"{comparison} vs WT_SD", fontsize=12, fontweight="bold")
    ax.set_ylabel("Number of events")

    for i, (u, d) in enumerate(zip(up_counts, down_counts)):
        if u > 0:
            ax.text(i, u, str(int(u)), ha="center", va="bottom", fontsize=8)
        if d > 0:
            ax.text(i, -d, str(int(d)), ha="center", va="top", fontsize=8)

axes[0, 0].legend(loc="upper right")
plt.tight_layout(rect=[0, 0, 1, 0.95])
png_path = f"{FIGURES_DIR}/significant_events_by_type.png"
svg_path = f"{FIGURES_DIR}/significant_events_by_type.svg"
plt.savefig(png_path, dpi=150)
plt.savefig(svg_path)
plt.close(fig)
print(f"Saved: {png_path}")
print(f"Saved: {svg_path}")

# --- Figure 2: significant unique-gene count per timepoint, dot plot ---
gene_counts = {}
for comparison in COMPARISON_ORDER:
    comp_df = sig_df[sig_df["comparison"] == f"{comparison}_vs_WT_SD"]
    gene_counts[comparison] = comp_df["gene_id"].nunique()

labels = [c.replace("_WHFD", "") for c in COMPARISON_ORDER]  # W8, W18, W24, W42
values = [gene_counts[c] for c in COMPARISON_ORDER]

fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(labels, values, linestyle="-", linewidth=1.2, color="grey", zorder=1)
ax.scatter(labels, values, s=140, color="orchid", zorder=2, edgecolors="black", linewidths=0.8)

for x, y in zip(labels, values):
    ax.annotate(str(y), (x, y), fontsize=11, fontweight="bold",
                ha="center", va="bottom", xytext=(0, 10), textcoords="offset points")

ax.set_xlabel("HFD timepoint (vs WT_SD baseline)", fontsize=11)
ax.set_ylabel("Number of unique genes with significant splicing change", fontsize=10)
ax.set_title("Significant splicing genes over time", fontsize=13, fontweight="bold")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_ylim(0, max(values) * 1.2)

plt.tight_layout()
png_path = f"{FIGURES_DIR}/significant_gene_count_over_time.png"
svg_path = f"{FIGURES_DIR}/significant_gene_count_over_time.svg"
plt.savefig(png_path, dpi=150)
plt.savefig(svg_path)
plt.close(fig)
print(f"Saved: {png_path}")
print(f"Saved: {svg_path}")

print("\nDone.")
