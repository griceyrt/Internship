#!/usr/bin/env python3
# =============================================================================
# Overlap #1 scatter: log2FoldChange (DESeq2) vs dPSI (SUPPA) for the genes
# that overlap between the two, testing whether splicing-change size
# predicts expression-change size.
# Author: Gricey
#
# Input (produced by overlap1_de_vs_dpsi.R, nothing computed here that
# isn't already a pipeline output): results/overlap1/overlap1_gene_list.tsv
# One row per (gene, event_type) overlap pair -- a gene significant in 2
# different event types gets 2 rows, which is expected given this
# project's per-event-type-separate convention, not a duplication bug.
#
# Requires: pandas, matplotlib
#   pip install pandas matplotlib
# =============================================================================

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"]

BASE = "/Users/gricey/Desktop/Internship/fmr1_polya_splicing/results/overlap1"
OUT_DIR = BASE

df = pd.read_csv(f"{BASE}/overlap1_gene_list.tsv", sep="\t")
print(f"Loaded {len(df)} (gene, event_type) overlap rows, "
      f"{df['gene_id'].nunique()} unique genes.")

EVENT_TYPES = ["SE", "A3", "A5", "MX", "RI", "AF", "AL"]
# Okabe-Ito colorblind-safe palette, consistent with earlier figures in
# this project that needed a discrete categorical palette
COLORS = ["#E69F00", "#56B4E9", "#009E73", "#F0E442",
          "#0072B2", "#D55E00", "#CC79A7"]
color_map = dict(zip(EVENT_TYPES, COLORS))

fig, ax = plt.subplots(figsize=(7, 6))

for event in EVENT_TYPES:
    sub = df[df["event_type"] == event]
    if sub.empty:
        continue
    ax.scatter(sub["dpsi_value"], sub["log2FoldChange"],
               label=f"{event} (n={len(sub)})", color=color_map[event], s=50,
               edgecolor="white", linewidth=0.5, alpha=0.85)

# Pearson correlation across ALL overlap rows pooled (matches how this
# project reports correlations elsewhere, e.g. isoforms-vs-entropy) --
# reported here as a simple summary stat, not a claim about any single
# event type's own relationship
corr = df["dpsi_value"].corr(df["log2FoldChange"])
n = len(df)
ax.text(0.02, 0.98, f"Pearson R = {corr:.2f} (n={n})",
        transform=ax.transAxes, ha="left", va="top", fontsize=10,
        style="italic")

ax.axhline(0, color="grey", linewidth=0.5, linestyle="--")
ax.axvline(0, color="grey", linewidth=0.5, linestyle="--")
ax.set_xlabel("dPSI (SUPPA)")
ax.set_ylabel("log2FoldChange (DESeq2)")
ax.set_title("Overlap #1: splicing change vs expression change", fontsize=13, fontweight="bold")
ax.legend(title="Event type", fontsize=8, title_fontsize=9, loc="best")
for spine in ("top", "right"):
    ax.spines[spine].set_visible(False)

fig.tight_layout()

png_path = f"{OUT_DIR}/overlap1_scatter_dpsi_vs_log2fc.png"
svg_path = f"{OUT_DIR}/overlap1_scatter_dpsi_vs_log2fc.svg"
fig.savefig(png_path, dpi=200, bbox_inches="tight")
fig.savefig(svg_path, bbox_inches="tight")
plt.close(fig)
print("Saved:", png_path)
print("Saved:", svg_path)
