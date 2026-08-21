#!/usr/bin/env python3
# =============================================================================
# 2x2(x2) comparison: gene-level vs isoform-level DESeq2, padj vs raw p
# threshold, with/without a 2.5-fold cutoff -- checks whether the
# asymmetric up/down split in Panel A (gene-level, padj<0.05) is explained
# by threshold type, analysis level, or both, against Shin/Chikhaoui et al.
# 2022's Fig 3F (isoform-level, raw-p-style). Holds mm39 + oarfish fixed
# throughout so only level/threshold/cutoff vary, isolating those factors
# from genome build and quantifier choice.
#
# Style matches GSE133398_make_final_figures.py (project 2 precedent):
# Okabe-Ito up/down colors, adjustText labels, Arial.
#
# Deliberately no log2FC cutoff on tests 1-4, matching the original Panel A
# convention -- tests 5-8 add |log2FC| > log2(2.5) on top of the same 4
# combinations, so the fold-change effect can be seen in isolation too.
#
# Inputs:
#   results/normalisation/deseq2_KO_vs_WT.tsv (gene-level, existing)
#   results/de_2x2_test/deseq2_isoform_level_KO_vs_WT.tsv (isoform-level, new)
#   results/isoform_proportions/entropy_shift_ranked_min10TPM.tsv (for
#     gene_id -> gene_name labels only, same source used for Overlap #1's
#     Gene_List tab)
#
# Requires: pandas, matplotlib, adjustText
#   pip install pandas matplotlib adjustText
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
GENE_PATH = f"{BASE}/results/normalisation/deseq2_KO_vs_WT.tsv"
ISO_PATH = f"{BASE}/results/de_2x2_test/deseq2_isoform_level_KO_vs_WT.tsv"
NAME_LOOKUP_PATH = f"{BASE}/results/isoform_proportions/entropy_shift_ranked_min10TPM.tsv"
OUT_DIR = f"{BASE}/results/de_2x2_test/figures"
os.makedirs(OUT_DIR, exist_ok=True)

gene_df = pd.read_csv(GENE_PATH, sep="\t")
iso_df = pd.read_csv(ISO_PATH, sep="\t")
name_lookup = pd.read_csv(NAME_LOOKUP_PATH, sep="\t")[["gene_id", "gene_name"]].drop_duplicates()

gene_df = gene_df.merge(name_lookup, on="gene_id", how="left")
iso_df = iso_df.merge(name_lookup, on="gene_id", how="left")

UP_COLOR, DOWN_COLOR, NS_COLOR = "#D55E00", "#0072B2", "#CCCCCC"


def make_volcano(df, sig_col, level_name, threshold_name, out_name, fc_cutoff=None):
    d = df[df[sig_col].notna() & df["log2FoldChange"].notna()].copy()
    d["category"] = "ns"
    passes_fc = d.log2FoldChange.abs() > fc_cutoff if fc_cutoff else True
    d.loc[(d.log2FoldChange > 0) & (d[sig_col] < 0.05) & passes_fc, "category"] = "up"
    d.loc[(d.log2FoldChange < 0) & (d[sig_col] < 0.05) & passes_fc, "category"] = "down"
    d["neg_log10_sig"] = -np.log10(d[sig_col].clip(lower=1e-300))

    n_up, n_down = (d.category == "up").sum(), (d.category == "down").sum()

    fig, ax = plt.subplots(figsize=(4.6, 4.2))
    ns, up, down = d[d.category == "ns"], d[d.category == "up"], d[d.category == "down"]

    ax.scatter(ns.log2FoldChange, ns.neg_log10_sig, c=NS_COLOR, s=3, alpha=0.3, linewidths=0, zorder=1, rasterized=True)
    ax.scatter(up.log2FoldChange, up.neg_log10_sig, c=UP_COLOR, s=10, alpha=0.85, linewidths=0, zorder=2)
    ax.scatter(down.log2FoldChange, down.neg_log10_sig, c=DOWN_COLOR, s=10, alpha=0.85, linewidths=0, zorder=2)

    ax.axhline(-np.log10(0.05), color="#BBBBBB", ls="--", lw=0.7)
    ax.axvline(0, color="#BBBBBB", ls="-", lw=0.5)
    if fc_cutoff:
        ax.axvline(fc_cutoff, color="#BBBBBB", ls="--", lw=0.7)
        ax.axvline(-fc_cutoff, color="#BBBBBB", ls="--", lw=0.7)

    top_hits = pd.concat([up.nsmallest(5, sig_col), down.nsmallest(5, sig_col)])
    texts = []
    for _, row in top_hits.iterrows():
        ax.scatter(row.log2FoldChange, row.neg_log10_sig, s=25, facecolors="none",
                   edgecolors="black", linewidths=0.5, zorder=3)
        label = row["gene_name"] if pd.notna(row.get("gene_name", np.nan)) else row.iloc[0]
        texts.append(ax.text(row.log2FoldChange, row.neg_log10_sig, label, fontsize=6, fontstyle="italic"))
    if texts:
        adjust_text(texts, ax=ax, expand=(2.8, 3.2), force_text=(1.8, 2.8), force_points=(1.3, 1.3),
                    arrowprops=dict(arrowstyle="-", color="#888888", lw=0.5, shrinkA=3, shrinkB=3))

    sig_label = "padj" if sig_col == "padj" else "p"
    fc_label = f", |log2FC|>{fc_cutoff:.2f}" if fc_cutoff else ""
    ax.set_xlabel("log$_2$FC (KO / WT)")
    ax.set_ylabel(f"$-$log$_{{10}}$({sig_label})")
    ax.set_title(f"{level_name}, {threshold_name}{fc_label}\nup={n_up}, down={n_down}", fontsize=9)

    for ext in ["svg", "png"]:
        fig.savefig(f"{OUT_DIR}/{out_name}.{ext}", format=ext, bbox_inches="tight", dpi=300 if ext == "png" else None)
    plt.close(fig)
    print(f"{out_name}: up={n_up}, down={n_down}")
    return n_up, n_down


import math
FC_CUTOFF = math.log2(2.5)  # ~1.32, matches the paper's stated "2.5-fold" threshold

rows = []
# Tests 1-4: significance threshold only, no fold-change cutoff (original 2x2)
rows.append(("gene", "padj", "none") + make_volcano(gene_df, "padj", "Gene-level", "padj<0.05", "1_volcano_gene_padj"))
rows.append(("gene", "rawp", "none") + make_volcano(gene_df, "pvalue", "Gene-level", "raw p<0.05", "2_volcano_gene_rawp"))
rows.append(("isoform", "padj", "none") + make_volcano(iso_df, "padj", "Isoform-level", "padj<0.05", "3_volcano_isoform_padj"))
rows.append(("isoform", "rawp", "none") + make_volcano(iso_df, "pvalue", "Isoform-level", "raw p<0.05", "4_volcano_isoform_rawp"))

# Tests 5-8: same 4 combinations, now with |log2FC| > log2(2.5) added on top --
# test 8 (isoform + rawp + FC) is the one directly comparable to the paper's
# "~100 RNAs altered by at least 2.5-fold" statement.
rows.append(("gene", "padj", "fc2.5") + make_volcano(gene_df, "padj", "Gene-level", "padj<0.05", "5_volcano_gene_padj_fc", fc_cutoff=FC_CUTOFF))
rows.append(("gene", "rawp", "fc2.5") + make_volcano(gene_df, "pvalue", "Gene-level", "raw p<0.05", "6_volcano_gene_rawp_fc", fc_cutoff=FC_CUTOFF))
rows.append(("isoform", "padj", "fc2.5") + make_volcano(iso_df, "padj", "Isoform-level", "padj<0.05", "7_volcano_isoform_padj_fc", fc_cutoff=FC_CUTOFF))
rows.append(("isoform", "rawp", "fc2.5") + make_volcano(iso_df, "pvalue", "Isoform-level", "raw p<0.05", "8_volcano_isoform_rawp_fc", fc_cutoff=FC_CUTOFF))

summary = pd.DataFrame(rows, columns=["level", "threshold", "fc_cutoff", "n_up", "n_down"])
summary.insert(0, "test", range(1, len(summary) + 1))
summary.to_csv(f"{BASE}/results/de_2x2_test/2x2_volcano_summary_counts.tsv", sep="\t", index=False)
print("\n", summary)
