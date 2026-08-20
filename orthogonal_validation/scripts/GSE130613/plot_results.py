#!/usr/bin/env python3
"""
Orthogonal Validation — Steps 11 & 13
Author: Gricey

Step 11: dPSI volcano plot  (Illumina SUPPA2 diffSplice, |dPSI|>0.1, p<0.3)
Step 13: DE volcano plot    (DESeq2, |log2FC|>0.5, padj<0.05)

Run from: orthogonal_validation/
"""

import os
import glob
import re
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from adjustText import adjust_text

# STYLE — DESeq2-inspired clean aesthetic
# White background, open axes (no top/right spines), subtle grid behind data.
matplotlib.rcParams.update({
    "font.family":          "sans-serif",
    "font.sans-serif":      ["Arial", "DejaVu Sans", "Helvetica", "Verdana"],
    "font.size":            11,
    "axes.spines.top":      False,
    "axes.spines.right":    False,
    "axes.spines.left":     True,
    "axes.spines.bottom":   True,
    "axes.edgecolor":       "#AAAAAA",
    "axes.linewidth":       0.8,
    "axes.grid":            True,
    "grid.color":           "#E5E5E5",
    "grid.linewidth":       0.6,
    "grid.linestyle":       "--",
    "axes.axisbelow":       True,
    "figure.facecolor":     "white",
    "axes.facecolor":       "white",
    "xtick.color":          "#444444",
    "ytick.color":          "#444444",
    "xtick.labelsize":      10,
    "ytick.labelsize":      10,
    "axes.labelcolor":      "#222222",
    "axes.titlepad":        10,
    "legend.frameon":       True,
    "legend.framealpha":    0.9,
    "legend.edgecolor":     "#DDDDDD",
})

# PATHS
BASE_DIR    = "/Users/gricey/Desktop/Internship/orthogonal_validation"
GTF_PATH    = "/Users/gricey/Desktop/Internship/boundary_analysis/data/transcriptome_productivity.gtf"
SUPPA_DATE  = "2026-06-20"
SUPPA_DIR   = os.path.join(BASE_DIR, "results", "GSE130613", f"suppa_{SUPPA_DATE}", "diff")
SUPPA_PSI   = os.path.join(BASE_DIR, "results", "GSE130613", f"suppa_{SUPPA_DATE}")
DESEQ_FILE  = os.path.join(BASE_DIR, "results", "GSE130613", "normalisation", "deseq2_KO_vs_WT.tsv")
OUT_DIR     = os.path.join(BASE_DIR, "results", "GSE130613", "figures")
os.makedirs(OUT_DIR, exist_ok=True)

# STEP 11 — dPSI Volcano Plot
# Thresholds: |dPSI| > 0.1, p < 0.3
print("=== Step 11: dPSI Volcano Plot ===")

DPSI_THRESH  = 0.1
PVAL_SPLICE  = 0.3

EVENT_COLORS = {
    "A3": "#E69F00",
    "A5": "#56B4E9",
    "AF": "#009E73",
    "AL": "#F0E442",
    "MX": "#0072B2",
    "RI": "#D55E00",
    "SE": "#CC79A7",
}

dfs = []
for dpsi_file in sorted(glob.glob(os.path.join(SUPPA_DIR, "diff_*_strict.dpsi.temp.0"))):
    event_type = os.path.basename(dpsi_file).replace("diff_", "").replace("_strict.dpsi.temp.0", "")
    df = pd.read_csv(dpsi_file, sep="\t", index_col=0)
    df.columns = ["dPSI", "pval"]
    df["event_type"] = event_type
    dfs.append(df)

dpsi_df = pd.concat(dfs)
dpsi_df["dPSI"] = -dpsi_df["dPSI"]   # negate: file is WT-KO, we want KO-WT
dpsi_df = dpsi_df.dropna(subset=["dPSI", "pval"])
dpsi_df = dpsi_df[dpsi_df["pval"] > 0]
dpsi_df["neg_log10_p"] = -np.log10(dpsi_df["pval"])
dpsi_df["significant"] = (abs(dpsi_df["dPSI"]) > DPSI_THRESH) & (dpsi_df["pval"] < PVAL_SPLICE)

print(f"Total events with data : {len(dpsi_df)}")
print(f"Significant events     : {dpsi_df['significant'].sum()}")

fig, ax = plt.subplots(figsize=(8, 6))

ns = dpsi_df[~dpsi_df["significant"]]
ax.scatter(ns["dPSI"], ns["neg_log10_p"],
           color="lightgrey", alpha=0.4, s=15, linewidths=0, zorder=1)

sig = dpsi_df[dpsi_df["significant"]]
for etype, color in EVENT_COLORS.items():
    sub = sig[sig["event_type"] == etype]
    if len(sub) > 0:
        ax.scatter(sub["dPSI"], sub["neg_log10_p"],
                   color=color, alpha=0.9, s=40, linewidths=0,
                   label=f"{etype} (n={len(sub)})", zorder=2)

ax.axvline(x= DPSI_THRESH,  color="black", linestyle="--", linewidth=0.8, alpha=0.6)
ax.axvline(x=-DPSI_THRESH,  color="black", linestyle="--", linewidth=0.8, alpha=0.6)
ax.axhline(y=-np.log10(PVAL_SPLICE), color="black", linestyle="--", linewidth=0.8, alpha=0.6)

ax.set_xlabel("dPSI (PerDKO − WT)", fontsize=12)
ax.set_ylabel("−log₁₀(p-value)", fontsize=12)
ax.set_title("Differentially spliced events — Illumina only (PerDKO vs WT, CT16-20, Liver)\nGSE130613 — Illumina short-read", fontsize=11)
ax.legend(title="Event type", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=9)
plt.tight_layout()

out_path = os.path.join(OUT_DIR, "step11_dPSI_volcano.png")
plt.savefig(out_path, dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved: {out_path}")

# STEP 13 — DE Volcano Plot
# Thresholds: |log2FoldChange| > 0.5, padj < 0.05
print("\n=== Step 13: DE Volcano Plot ===")

LFC_THRESH  = 0.5
PADJ_THRESH = 0.05

# Build gene_id -> gene_name map
gene_map = {}
with open(GTF_PATH) as fh:
    for line in fh:
        if "\tgene\t" not in line:
            continue
        gid   = re.search(r'gene_id "([^"]+)"', line)
        gname = re.search(r'gene_name "([^"]+)"', line)
        if gid and gname:
            gene_map[gid.group(1)] = gname.group(1)
name_to_id = {v: k for k, v in gene_map.items()}

# Genes to label (selected by Kiran)
CIRCADIAN_GENES = [
    "Ren1", "Orm2", "Prtn3", "Cyp2b9", "Dbp", "Cxcl1", "Saa1", "Ighg2b",
    "Celf4", "Egr1", "Ciart", "Nr1d1", "Nr1d2", "Btg2", "Per3", "Ephx1",
    "Gstm2", "Clock", "Cyp2a5", "Arntl", "Cry1", "Kmt2d", "Gstm3",
    "Prkca", "Cyp7a1", "Fancm", "Cspg5"
]

# New stageR output: transcript-level — collapse to gene level for volcano
# Use padj_gene_stageR (stage 1 gene screening) for significance
# Take the transcript with largest |log2FC| per gene as representative
de_tx = pd.read_csv(DESEQ_FILE, sep="\t")
de_tx = de_tx.dropna(subset=["log2FoldChange", "padj_gene_stageR"])
de_tx = de_tx[de_tx["padj_gene_stageR"] > 0]

# Collapse: one row per gene (representative = max |log2FC| transcript)
de = (de_tx.sort_values("log2FoldChange", key=abs, ascending=False)
           .drop_duplicates(subset="gene_id")
           .set_index("gene_id"))

de["neg_log10_padj"] = -np.log10(de["padj_gene_stageR"])

de["category"] = "not significant"
de.loc[(de["log2FoldChange"] >  LFC_THRESH) & (de["padj_gene_stageR"] < PADJ_THRESH), "category"] = "up in PerDKO"
de.loc[(de["log2FoldChange"] < -LFC_THRESH) & (de["padj_gene_stageR"] < PADJ_THRESH), "category"] = "down in PerDKO"

cat_colors = {
    "not significant": "lightgrey",
    "up in PerDKO":    "#D55E00",
    "down in PerDKO":  "#0072B2",
}

print(f"Total genes: {len(de)}")
for cat, count in de["category"].value_counts().items():
    print(f"  {cat}: {count}")

fig, ax = plt.subplots(figsize=(13, 8))
ax.grid(False)   # no gridlines — consistent with figures 12 and 15

# Cap y-axis just above the highest significant point (~12.5)
Y_CAP = 13.5
de["neg_log10_padj_clipped"] = de["neg_log10_padj"].clip(upper=Y_CAP)
de["clipped"] = de["neg_log10_padj"] > Y_CAP

for cat, color in cat_colors.items():
    sub = de[(de["category"] == cat) & (~de["clipped"])]
    # not significant: tiny, very light — same as figs 12/15
    alpha = 0.85  if cat != "not significant" else 0.3
    size  = 30    if cat != "not significant" else 6
    c     = color if cat != "not significant" else "#CCCCCC"
    ax.scatter(sub["log2FoldChange"], sub["neg_log10_padj_clipped"],
               color=c, alpha=alpha, s=size, linewidths=0,
               label=f"{cat} (n={len(de[de['category']==cat])})",
               zorder=2 if cat != "not significant" else 1)
    sub_clip = de[(de["category"] == cat) & (de["clipped"])]
    if len(sub_clip) > 0:
        ax.scatter(sub_clip["log2FoldChange"], [Y_CAP - 0.3] * len(sub_clip),
                   color=color, alpha=alpha, s=50, marker="^", linewidths=0, zorder=3)

ax.axvline(x= LFC_THRESH,  color="#BBBBBB", linestyle="--", linewidth=0.8)
ax.axvline(x=-LFC_THRESH,  color="#BBBBBB", linestyle="--", linewidth=0.8)
ax.axhline(y=-np.log10(PADJ_THRESH), color="#BBBBBB", linestyle="--", linewidth=0.8)
ax.text(0.01, 0.98, f"▲ clipped (y > {Y_CAP})", transform=ax.transAxes,
        fontsize=7, color="#AAAAAA", va="top")

# Label selected genes — only if significant
# gene_name is now a column in de; build reverse map from gene_name -> gene_id
name_to_id_stageR = {row["gene_name"]: gid for gid, row in de.iterrows()
                     if pd.notna(row.get("gene_name", None))}
# Separate genes into blue (down) and orange (up) for different placement strategies
blue_genes = []   # (gene, x, y) — will be stacked in a left-side column
orange_texts = []

for gene in CIRCADIAN_GENES:
    gid = name_to_id_stageR.get(gene) or name_to_id.get(gene)
    if gid and gid in de.index:
        if de.loc[gid, "category"] == "not significant":
            continue
        x = de.loc[gid, "log2FoldChange"]
        y = de.loc[gid, "neg_log10_padj_clipped"]
        dot_color = cat_colors.get(de.loc[gid, "category"], "lightgrey")
        ax.scatter(x, y, color=dot_color, s=70, zorder=5,
                   edgecolors="black", linewidths=0.8)
        if y < 1.35:
            continue
        if de.loc[gid, "category"] == "down in PerDKO":
            blue_genes.append((gene, x, y))
        else:
            orange_texts.append(ax.text(
                x + 1.5, y, gene,
                fontsize=8.5, fontstyle="italic", zorder=6, clip_on=False,
            ))

# Orange: use adjust_text to find good positions, then redraw as ax.annotate
# so lines reliably touch both the dot and the label box.
adjust_text(
    orange_texts, ax=ax,
    expand=(2.0, 2.0),
    force_text=(1.2, 1.2),
    force_points=(1.0, 1.0),
    lim=500,
)

# Collect final positions, remove the text objects, redraw as annotate
orange_gene_list = []
for gene in CIRCADIAN_GENES:
    gid = name_to_id_stageR.get(gene) or name_to_id.get(gene)
    if gid and gid in de.index and de.loc[gid, "category"] == "up in PerDKO":
        x = de.loc[gid, "log2FoldChange"]
        y = de.loc[gid, "neg_log10_padj_clipped"]
        if y >= 1.35:
            orange_gene_list.append((gene, x, y))

for txt, (gene, dot_x, dot_y) in zip(orange_texts, orange_gene_list):
    lx, ly = txt.get_position()
    txt.remove()
    ax.annotate(
        gene,
        xy=(dot_x, dot_y),
        xytext=(lx, ly),
        fontsize=8.5, fontstyle="italic", zorder=6,
        ha=txt.get_ha(), va="center",
        arrowprops=dict(arrowstyle="-", color="#888888", lw=0.7,
                        shrinkA=4, shrinkB=4),
    )

# Blue: place each label just to the left of its own dot
blue_genes_sorted = sorted(blue_genes, key=lambda t: t[2], reverse=True)

MIN_YGAP = 0.7
label_positions = []
for gene, dot_x, dot_y in blue_genes_sorted:
    lx = dot_x - 2.0
    ly = dot_y
    for _, prev_lx, prev_ly in label_positions:
        if abs(ly - prev_ly) < MIN_YGAP:
            ly = prev_ly - MIN_YGAP
    label_positions.append((gene, lx, ly))

for (gene, dot_x, dot_y), (_, lx, ly) in zip(blue_genes_sorted, label_positions):
    ax.annotate(
        gene,
        xy=(dot_x, dot_y),
        xytext=(lx, ly),
        fontsize=8.5, fontstyle="italic", zorder=6,
        ha="right", va="center",
        arrowprops=dict(arrowstyle="-", color="#888888", lw=0.7,
                        shrinkA=4, shrinkB=4),
    )

ax.set_xlabel("log₂ fold change (PerDKO / WT)", fontsize=12)
ax.set_ylabel("−log₁₀(padj gene, stageR)", fontsize=12)
ax.set_title("Differentially expressed genes — Illumina only (PerDKO vs WT, CT16-20, Liver)\nGSE130613 — DESeq2 + stageR, transcript-level", fontsize=11)
ax.set_xlim(-13, 12)
ax.set_ylim(0, Y_CAP + 0.3)
ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=9)
plt.tight_layout()

out_path = os.path.join(OUT_DIR, "step13_DE_volcano.png")
plt.savefig(out_path, dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved: {out_path}")

# MAPPING RATES — Dumbbell plot (before vs after trimming)
# Shows per-sample Salmon mapping rates before and after cutadapt trimming.
# Open circle = untrimmed, filled circle = trimmed.
# Colors match figure 13: WT = blue, PerDKO = orange.
print("\n=== Mapping rates dumbbell plot ===")

MAPPING = [
    # (label,      condition,  untrimmed, trimmed)
    ("WT a",      "WT",       56.71,     82.32),
    ("WT b",      "WT",       57.22,     82.28),
    ("WT c",      "WT",       70.17,     87.86),
    ("WT d",      "WT",       56.89,     82.22),
    ("PerDKO a",  "PerDKO",   47.26,     76.36),
    ("PerDKO b",  "PerDKO",   46.62,     75.35),
    ("PerDKO c",  "PerDKO",   43.82,     73.00),
    ("PerDKO d",  "PerDKO",   47.54,     76.43),
]

COLORS = {"WT": "#0072B2", "PerDKO": "#D55E00"}

fig, ax = plt.subplots(figsize=(8, 5))

n = len(MAPPING)
# Plot bottom to top: PerDKO first (lower), WT above, with a gap between groups
y_positions = list(range(n))   # 0..7, reversed below

for i, (label, cond, before, after) in enumerate(reversed(MAPPING)):
    y   = i
    col = COLORS[cond]
    # Add subtle gap between WT and PerDKO groups
    if i >= 4:
        y += 0.5

    # Connecting line
    ax.plot([before, after], [y, y], color=col, linewidth=1.8,
            alpha=0.7, zorder=1)
    # Untrimmed — open circle
    ax.scatter(before, y, color="white", edgecolors=col,
               s=80, linewidths=1.8, zorder=3)
    # Trimmed — filled circle
    ax.scatter(after, y, color=col, s=80, linewidths=0, zorder=3)

    # Improvement label
    delta = after - before
    ax.text(after + 0.8, y, f"+{delta:.1f}%",
            va="center", ha="left", fontsize=8.5, color=col)

    # Sample label on y-axis
    ax.text(-1.5, y, label, va="center", ha="right", fontsize=9,
            color="#333333")

ax.set_xlabel("Salmon mapping rate (%)", fontsize=11)
ax.set_title("Mapping rates before and after trimming\nGSE130613 — Salmon quasi-mapping", fontsize=11)
ax.set_xlim(30, 100)
ax.set_ylim(-0.6, n + 0.6)
ax.set_yticks([])   # labels drawn manually above

# Reference line at 60% and 80%
for ref in [60, 70, 80]:
    ax.axvline(ref, color="#CCCCCC", linewidth=0.8, linestyle=":", zorder=0)

# Legend
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0],[0], marker="o", color="w", markerfacecolor="white",
           markeredgecolor="#555555", markeredgewidth=1.5,
           markersize=8, label="Untrimmed"),
    Line2D([0],[0], marker="o", color="w", markerfacecolor="#555555",
           markersize=8, label="Trimmed"),
    Line2D([0],[0], marker="o", color="w", markerfacecolor=COLORS["WT"],
           markersize=8, label="WT"),
    Line2D([0],[0], marker="o", color="w", markerfacecolor=COLORS["PerDKO"],
           markersize=8, label="PerDKO"),
]
ax.legend(handles=legend_elements, fontsize=9,
          bbox_to_anchor=(1.02, 1), loc="upper left")

plt.tight_layout()
out_path = os.path.join(OUT_DIR, "mapping_rates.png")
plt.savefig(out_path, dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved: {out_path}")

print("\n=== DONE ===")
print(f"Figures saved in: {OUT_DIR}")
