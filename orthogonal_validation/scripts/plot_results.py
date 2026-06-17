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

# =============================================================================
# PATHS
# =============================================================================
BASE_DIR    = "/Users/gricey/Desktop/Internship/orthogonal_validation"
GTF_PATH    = "/Users/gricey/Desktop/Internship/boundary_analysis/data/transcriptome_productivity.gtf"
SUPPA_DATE  = "2026-06-17"
SUPPA_DIR   = os.path.join(BASE_DIR, "results", f"suppa_{SUPPA_DATE}", "diff")
DESEQ_FILE  = os.path.join(BASE_DIR, "results", "normalisation", "deseq2_KO_vs_WT.tsv")
OUT_DIR     = os.path.join(BASE_DIR, "results", "figures")
os.makedirs(OUT_DIR, exist_ok=True)

# =============================================================================
# STEP 11 — dPSI Volcano Plot
# Thresholds: |dPSI| > 0.1, p < 0.3
# =============================================================================
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
ax.set_title("Differential Splicing: PerDKO vs WT (CT16-20, Liver)\nGSE130613 — Illumina short-read", fontsize=11)
ax.legend(title="Event type", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=9)
plt.tight_layout()

out_path = os.path.join(OUT_DIR, "step11_dPSI_volcano.png")
plt.savefig(out_path, dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved: {out_path}")

# =============================================================================
# STEP 13 — DE Volcano Plot
# Thresholds: |log2FoldChange| > 0.5, padj < 0.05
# =============================================================================
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

de = pd.read_csv(DESEQ_FILE, sep="\t", index_col=0)
de = de.dropna(subset=["log2FoldChange", "padj"])
de = de[de["padj"] > 0]
de["neg_log10_padj"] = -np.log10(de["padj"])
de["gene_name"] = de.index.map(gene_map)

de["category"] = "not significant"
de.loc[(de["log2FoldChange"] >  LFC_THRESH) & (de["padj"] < PADJ_THRESH), "category"] = "up in PerDKO"
de.loc[(de["log2FoldChange"] < -LFC_THRESH) & (de["padj"] < PADJ_THRESH), "category"] = "down in PerDKO"

cat_colors = {
    "not significant": "lightgrey",
    "up in PerDKO":    "#D55E00",
    "down in PerDKO":  "#0072B2",
}

print(f"Total genes: {len(de)}")
for cat, count in de["category"].value_counts().items():
    print(f"  {cat}: {count}")

fig, ax = plt.subplots(figsize=(11, 8))
for cat, color in cat_colors.items():
    sub = de[de["category"] == cat]
    alpha = 0.4 if cat == "not significant" else 0.85
    size  = 10  if cat == "not significant" else 30
    ax.scatter(sub["log2FoldChange"], sub["neg_log10_padj"],
               color=color, alpha=alpha, s=size, linewidths=0,
               label=f"{cat} (n={len(sub)})", zorder=2 if cat != "not significant" else 1)

ax.axvline(x= LFC_THRESH,  color="black", linestyle="--", linewidth=0.8, alpha=0.6)
ax.axvline(x=-LFC_THRESH,  color="black", linestyle="--", linewidth=0.8, alpha=0.6)
ax.axhline(y=-np.log10(PADJ_THRESH), color="black", linestyle="--", linewidth=0.8, alpha=0.6)

# Label selected genes — keep original dot color, add white-box label + connector
texts = []
for gene in CIRCADIAN_GENES:
    gid = name_to_id.get(gene)
    if gid and gid in de.index:
        x = de.loc[gid, "log2FoldChange"]
        y = de.loc[gid, "neg_log10_padj"]
        dot_color = cat_colors.get(de.loc[gid, "category"], "lightgrey")
        ax.scatter(x, y, color=dot_color, s=70, zorder=5,
                   edgecolors="black", linewidths=0.8)
        texts.append(ax.text(
            x, y, gene,
            fontsize=9, fontweight="bold", zorder=6, clip_on=False,
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.75)
        ))

adjust_text(
    texts, ax=ax,
    expand=(2.2, 2.2),
    force_text=(1.0, 1.5),
    force_points=(0.4, 0.6),
    arrowprops=dict(arrowstyle="-", color="#333333", lw=0.8)
)

ax.set_xlabel("log₂ fold change (PerDKO / WT)", fontsize=12)
ax.set_ylabel("−log₁₀(adjusted p-value)", fontsize=12)
ax.set_title("Differential Expression: PerDKO vs WT (CT16-20, Liver)\nGSE130613 — Illumina short-read", fontsize=11)
ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=9)
x_min, x_max = ax.get_xlim()
y_min, y_max = ax.get_ylim()
ax.set_xlim(x_min - 0.8, x_max + 0.8)
ax.set_ylim(y_min, y_max + 1.5)
plt.tight_layout()

out_path = os.path.join(OUT_DIR, "step13_DE_volcano.png")
plt.savefig(out_path, dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved: {out_path}")

print("\n=== DONE ===")
print(f"Figures saved in: {OUT_DIR}")
