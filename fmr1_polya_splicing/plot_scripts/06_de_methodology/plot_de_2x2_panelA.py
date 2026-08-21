#!/usr/bin/env python3
# =============================================================================
# Panel-A-style bar charts for the same 2x2(x2) grid as plot_de_2x2_volcano.py
# -- reuses plot_panelA_de_summary.py's diverging-bar-with-biotype-split
# design, generalized across all 8 combinations in one run.
#
# At isoform level, each transcript inherits its parent gene's biotype (via
# the gene_id column merged in by deseq2_isoform_level_test.R's tx2gene
# join) -- multiple transcripts of the same gene correctly share one
# biotype, not treated as independent observations.
#
# Requires: pandas, matplotlib
# =============================================================================

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"]

BASE = "/Users/gricey/Desktop/Internship/fmr1_polya_splicing"
GENE_PATH = f"{BASE}/results/normalisation/deseq2_KO_vs_WT.tsv"
ISO_PATH = f"{BASE}/results/de_2x2_test/deseq2_isoform_level_KO_vs_WT.tsv"
BIOTYPE_PATH = f"{BASE}/data/reference/gene_biotype_mm39.tsv"
OUT_DIR = f"{BASE}/results/de_2x2_test/figures"
os.makedirs(OUT_DIR, exist_ok=True)

gene_df = pd.read_csv(GENE_PATH, sep="\t")
iso_df = pd.read_csv(ISO_PATH, sep="\t")
biotype = pd.read_csv(BIOTYPE_PATH, sep="\t")


def make_panelA(df, sig_col, level_name, threshold_name, out_name, fc_cutoff=None):
    sig = df[(df[sig_col].notna()) & (df[sig_col] < 0.05)].copy()
    if fc_cutoff:
        sig = sig[sig.log2FoldChange.abs() > fc_cutoff]
    sig = sig.merge(biotype, on="gene_id", how="left")
    n_missing = sig["gene_biotype"].isna().sum()
    if n_missing:
        print(f"WARNING [{out_name}]: {n_missing} significant rows had no biotype match")
    sig["biotype_group"] = sig["gene_biotype"].apply(
        lambda b: "Protein coding" if b == "protein_coding" else "Non coding"
    )
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
        ax.text(down_val - (max(counts.values.max(), 1) * 0.02), y, str(abs(int(down_val))),
                va="center", ha="right", fontsize=9)
        ax.text(up_val + (max(counts.values.max(), 1) * 0.02), y, str(int(up_val)),
                va="center", ha="left", fontsize=9)

    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_yticks([1, 0])
    ax.set_yticklabels(groups)
    fc_label = f", |log2FC|>{fc_cutoff:.2f}" if fc_cutoff else ""
    ax.set_xlabel("Down-regulated <-  -> Up-regulated  (count)")
    ax.set_title(f"{level_name}, {threshold_name}{fc_label} (n={len(sig)} total)", fontsize=11, fontweight="bold")
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()

    for ext in ["svg", "png"]:
        fig.savefig(f"{OUT_DIR}/{out_name}.{ext}", dpi=200 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)
    print(f"{out_name}: n={len(sig)}\n{counts}\n")
    return counts


import math
FC_CUTOFF = math.log2(2.5)  # ~1.32, matches the paper's stated "2.5-fold" threshold

# Tests 1-4: significance threshold only (original 2x2)
make_panelA(gene_df, "padj", "Gene-level", "padj<0.05", "1_panelA_gene_padj")
make_panelA(gene_df, "pvalue", "Gene-level", "raw p<0.05", "2_panelA_gene_rawp")
make_panelA(iso_df, "padj", "Isoform-level", "padj<0.05", "3_panelA_isoform_padj")
make_panelA(iso_df, "pvalue", "Isoform-level", "raw p<0.05", "4_panelA_isoform_rawp")

# Tests 5-8: same 4, with |log2FC| > log2(2.5) added -- test 8 (isoform +
# rawp + FC) is directly comparable to the paper's "~100 RNAs, >=2.5-fold" claim
make_panelA(gene_df, "padj", "Gene-level", "padj<0.05", "5_panelA_gene_padj_fc", fc_cutoff=FC_CUTOFF)
make_panelA(gene_df, "pvalue", "Gene-level", "raw p<0.05", "6_panelA_gene_rawp_fc", fc_cutoff=FC_CUTOFF)
make_panelA(iso_df, "padj", "Isoform-level", "padj<0.05", "7_panelA_isoform_padj_fc", fc_cutoff=FC_CUTOFF)
make_panelA(iso_df, "pvalue", "Isoform-level", "raw p<0.05", "8_panelA_isoform_rawp_fc", fc_cutoff=FC_CUTOFF)
