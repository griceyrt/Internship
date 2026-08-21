#!/usr/bin/env python3
# =============================================================================
# 2+3 replicate test: volcano + Panel-A-style bar chart for gene-level
# DESeq2 with WT_Rep1 excluded (2 WT + 3 KO, matching the paper's
# structure). Kept separate from de_2x2_test/ -- a different axis
# (replicate count), not the level x threshold x fold-change grid there.
#
# Same style as plot_de_2x2_volcano.py / plot_de_2x2_panelA.py, no log2FC
# cutoff, for direct comparison against tests 1 and 2 in that grid (same
# two thresholds on the full 3+3 dataset).
#
# Requires: pandas, numpy, matplotlib, adjustText
# =============================================================================

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from adjustText import adjust_text
import os

matplotlib.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 9,
    "axes.spines.top": False, "axes.spines.right": False, "axes.spines.left": True, "axes.spines.bottom": True,
    "axes.edgecolor": "#AAAAAA", "axes.linewidth": 0.7, "axes.grid": False, "figure.facecolor": "white",
    "axes.facecolor": "white", "xtick.color": "#444444", "ytick.color": "#444444",
    "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.labelsize": 8, "axes.labelcolor": "#222222",
    "xtick.major.size": 3, "ytick.major.size": 3, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
})

BASE = "/Users/gricey/Desktop/Internship/fmr1_polya_splicing"
DE_PATH = f"{BASE}/results/de_2plus3_test/deseq2_gene_level_2plus3_KO_vs_WT.tsv"
NAME_LOOKUP_PATH = f"{BASE}/results/isoform_proportions/entropy_shift_ranked_min10TPM.tsv"
BIOTYPE_PATH = f"{BASE}/data/reference/gene_biotype_mm39.tsv"
OUT_DIR = f"{BASE}/results/de_2plus3_test/figures"
os.makedirs(OUT_DIR, exist_ok=True)

df = pd.read_csv(DE_PATH, sep="\t")
name_lookup = pd.read_csv(NAME_LOOKUP_PATH, sep="\t")[["gene_id", "gene_name"]].drop_duplicates()
biotype = pd.read_csv(BIOTYPE_PATH, sep="\t")
df = df.merge(name_lookup, on="gene_id", how="left")

UP_COLOR, DOWN_COLOR, NS_COLOR = "#D55E00", "#0072B2", "#CCCCCC"


def make_volcano(sig_col, threshold_name, out_name):
    d = df[df[sig_col].notna() & df["log2FoldChange"].notna()].copy()
    d["category"] = "ns"
    d.loc[(d.log2FoldChange > 0) & (d[sig_col] < 0.05), "category"] = "up"
    d.loc[(d.log2FoldChange < 0) & (d[sig_col] < 0.05), "category"] = "down"
    d["neg_log10_sig"] = -np.log10(d[sig_col].clip(lower=1e-300))
    n_up, n_down = (d.category == "up").sum(), (d.category == "down").sum()

    fig, ax = plt.subplots(figsize=(4.6, 4.2))
    ns, up, down = d[d.category == "ns"], d[d.category == "up"], d[d.category == "down"]
    ax.scatter(ns.log2FoldChange, ns.neg_log10_sig, c=NS_COLOR, s=3, alpha=0.3, linewidths=0, zorder=1, rasterized=True)
    ax.scatter(up.log2FoldChange, up.neg_log10_sig, c=UP_COLOR, s=10, alpha=0.85, linewidths=0, zorder=2)
    ax.scatter(down.log2FoldChange, down.neg_log10_sig, c=DOWN_COLOR, s=10, alpha=0.85, linewidths=0, zorder=2)
    ax.axhline(-np.log10(0.05), color="#BBBBBB", ls="--", lw=0.7)
    ax.axvline(0, color="#BBBBBB", ls="-", lw=0.5)

    top_hits = pd.concat([up.nsmallest(5, sig_col), down.nsmallest(5, sig_col)])
    texts = []
    for _, row in top_hits.iterrows():
        ax.scatter(row.log2FoldChange, row.neg_log10_sig, s=25, facecolors="none",
                   edgecolors="black", linewidths=0.5, zorder=3)
        label = row["gene_name"] if pd.notna(row.get("gene_name", np.nan)) else row["gene_id"]
        texts.append(ax.text(row.log2FoldChange, row.neg_log10_sig, label, fontsize=6, fontstyle="italic"))
    if texts:
        adjust_text(texts, ax=ax, expand=(2.8, 3.2), force_text=(1.8, 2.8), force_points=(1.3, 1.3),
                    arrowprops=dict(arrowstyle="-", color="#888888", lw=0.5, shrinkA=3, shrinkB=3))

    sig_label = "padj" if sig_col == "padj" else "p"
    ax.set_xlabel("log$_2$FC (KO / WT)")
    ax.set_ylabel(f"$-$log$_{{10}}$({sig_label})")
    ax.set_title(f"Gene-level, {threshold_name}, 2WT+3KO\nup={n_up}, down={n_down}", fontsize=9)
    for ext in ["svg", "png"]:
        fig.savefig(f"{OUT_DIR}/{out_name}.{ext}", format=ext, bbox_inches="tight", dpi=300 if ext == "png" else None)
    plt.close(fig)
    print(f"{out_name}: up={n_up}, down={n_down}")
    return n_up, n_down


def make_panelA(sig_col, threshold_name, out_name):
    sig = df[(df[sig_col].notna()) & (df[sig_col] < 0.05)].copy()
    sig = sig.merge(biotype, on="gene_id", how="left")
    sig["biotype_group"] = sig["gene_biotype"].apply(lambda b: "Protein coding" if b == "protein_coding" else "Non coding")
    sig["direction"] = sig["log2FoldChange"].apply(lambda x: "Up" if x > 0 else "Down")
    counts = sig.groupby(["direction", "biotype_group"]).size().unstack(fill_value=0)

    fig, ax = plt.subplots(figsize=(6, 4))
    groups = ["Protein coding", "Non coding"]
    colors = {"Protein coding": "#2A9D8F", "Non coding": "#A8DADC"}
    for y, group in zip([1, 0], groups):
        down_val = -counts.loc["Down", group] if "Down" in counts.index and group in counts.columns else 0
        up_val = counts.loc["Up", group] if "Up" in counts.index and group in counts.columns else 0
        ax.barh(y, down_val, color=colors[group], edgecolor="black", linewidth=0.6)
        ax.barh(y, up_val, color=colors[group], edgecolor="black", linewidth=0.6)
        ax.text(down_val - (max(counts.values.max(), 1) * 0.02), y, str(abs(int(down_val))), va="center", ha="right", fontsize=9)
        ax.text(up_val + (max(counts.values.max(), 1) * 0.02), y, str(int(up_val)), va="center", ha="left", fontsize=9)

    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_yticks([1, 0])
    ax.set_yticklabels(groups)
    ax.set_xlabel("Down-regulated <-  -> Up-regulated  (count)")
    ax.set_title(f"Gene-level, {threshold_name}, 2WT+3KO (n={len(sig)} total)", fontsize=11, fontweight="bold")
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    for ext in ["svg", "png"]:
        fig.savefig(f"{OUT_DIR}/{out_name}.{ext}", dpi=200 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)
    print(f"{out_name}: n={len(sig)}\n{counts}\n")


rows = []
rows.append(("padj",) + make_volcano("padj", "padj<0.05", "volcano_gene_padj_2plus3"))
rows.append(("rawp",) + make_volcano("pvalue", "raw p<0.05", "volcano_gene_rawp_2plus3"))
make_panelA("padj", "padj<0.05", "panelA_gene_padj_2plus3")
make_panelA("pvalue", "raw p<0.05", "panelA_gene_rawp_2plus3")

summary = pd.DataFrame(rows, columns=["threshold", "n_up", "n_down"])
summary.to_csv(f"{BASE}/results/de_2plus3_test/2plus3_summary_counts.tsv", sep="\t", index=False)
print("\n", summary)
print("\nCompare against de_2x2_test's test 1 (gene/padj/3WT+3KO: 126 up, 16 down)")
print("and test 2 (gene/rawp/3WT+3KO: 593 up, 374 down).")
