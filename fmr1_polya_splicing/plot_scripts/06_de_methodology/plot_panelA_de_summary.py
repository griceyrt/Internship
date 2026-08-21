#!/usr/bin/env python3
# =============================================================================
# Panel A: Differentially expressed genes, WT vs Fmr1 KO, split by direction
# (up/down) and biotype (protein-coding vs non-coding). Style reference:
# Chikhaoui/Mamgain et al. PERIOD manuscript Fig 4A.
# Author: Gricey
#
# Inputs:
#   results/normalisation/deseq2_KO_vs_WT.tsv (gene_id, log2FoldChange, padj)
#   data/reference/gene_biotype_mm39.tsv (gene_id, gene_biotype -- produced
#     by extract_gene_biotype.sh from the same mm39 GTF used throughout
#     this project)
#
# Significance: padj < 0.05 (matches Overlap #1's own DESeq2 threshold --
# total significant count here should equal Overlap #1's n=142 as a
# consistency check).
#
# Biotype collapsed to a simple 2-way split (protein_coding vs everything
# else grouped as "non coding"), matching the reference figure's own
# 2-category split rather than the full GTF biotype list (lncRNA, TEC,
# pseudogene, miRNA etc. all fold into "non coding" here).
#
# Requires: pandas, matplotlib
# =============================================================================

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"]

BASE = "/Users/gricey/Desktop/Internship/fmr1_polya_splicing"
DESEQ_PATH = f"{BASE}/results/normalisation/deseq2_KO_vs_WT.tsv"
BIOTYPE_PATH = f"{BASE}/data/reference/gene_biotype_mm39.tsv"
OUT_DIR = f"{BASE}/results/panelAB_summary"

import os
os.makedirs(OUT_DIR, exist_ok=True)

deseq = pd.read_csv(DESEQ_PATH, sep="\t")
biotype = pd.read_csv(BIOTYPE_PATH, sep="\t")

sig = deseq[(deseq["padj"].notna()) & (deseq["padj"] < 0.05)].copy()
print(f"DESeq2-significant genes (padj<0.05): {len(sig)}")

sig = sig.merge(biotype, on="gene_id", how="left")
n_missing_biotype = sig["gene_biotype"].isna().sum()
if n_missing_biotype:
    print(f"WARNING: {n_missing_biotype} significant genes had no biotype match -- check gene_biotype_mm39.tsv coverage.")

sig["biotype_group"] = sig["gene_biotype"].apply(
    lambda b: "Protein coding" if b == "protein_coding" else "Non coding"
)
sig["direction"] = sig["log2FoldChange"].apply(lambda x: "Up" if x > 0 else "Down")

counts = sig.groupby(["direction", "biotype_group"]).size().unstack(fill_value=0)
print(counts)

fig, ax = plt.subplots(figsize=(6, 4))

groups = ["Protein coding", "Non coding"]
colors = {"Protein coding": "#2A9D8F", "Non coding": "#A8DADC"}
y_positions = [1, 0]

for y, group in zip(y_positions, groups):
    down_val = -counts.loc["Down", group] if "Down" in counts.index and group in counts.columns else 0
    up_val = counts.loc["Up", group] if "Up" in counts.index and group in counts.columns else 0
    ax.barh(y, down_val, color=colors[group], edgecolor="black", linewidth=0.6, label=group if y == 1 else None)
    ax.barh(y, up_val, color=colors[group], edgecolor="black", linewidth=0.6)
    ax.text(down_val - 2, y, str(abs(int(down_val))), va="center", ha="right", fontsize=9)
    ax.text(up_val + 2, y, str(int(up_val)), va="center", ha="left", fontsize=9)

ax.axvline(0, color="black", linewidth=0.8)
ax.set_yticks(y_positions)
ax.set_yticklabels(groups)
ax.set_xlabel("Down-regulated <-  -> Up-regulated  (number of genes)")
ax.set_title("Differentially expressed genes, WT vs Fmr1 KO\n(padj < 0.05)", fontsize=12, fontweight="bold")
for spine in ("top", "right", "left"):
    ax.spines[spine].set_visible(False)

fig.tight_layout()
png_path = f"{OUT_DIR}/panelA_de_summary.png"
svg_path = f"{OUT_DIR}/panelA_de_summary.svg"
fig.savefig(png_path, dpi=200, bbox_inches="tight")
fig.savefig(svg_path, bbox_inches="tight")
plt.close(fig)
print("Saved:", png_path)
print("Saved:", svg_path)
