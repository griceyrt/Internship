"""
Generates a volcano plot and an Illumina-vs-Nanopore concordance scatter for
every pipeline version tested in the MARS-seq spread investigation: GEO
official counts, 3'UTR-extended Salmon, STAR+htseq-count (extended GTF),
STAR+htseq-count (original GTF). All comparisons at GENE level (one row per
gene, representative = max |log2FC|), matched against Nanopore's DE_all
sheet. GEO's matrix is indexed by gene SYMBOL, matched on mgi_symbol;
everything else matches on gene_id. Threshold: padj<0.05, |log2FC|>0.5.
Author: Gricey

Run from: orthogonal_validation/
Outputs to: results/GSE130613/marsseq_sanity_check/investigation_figures/
"""
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT_DIR = "results/GSE130613/marsseq_sanity_check/investigation_figures"
os.makedirs(OUT_DIR, exist_ok=True)

LFC_THRESH = 0.5
PADJ_THRESH = 0.05

COLORS = {"not significant": "lightgrey", "up in PerDKO": "#D55E00", "down in PerDKO": "#0072B2"}


def gene_level_dedup(df, lfc_col, padj_col, id_col):
    """One row per gene: representative = transcript/row with max |log2FC|."""
    df = df.dropna(subset=[lfc_col, padj_col])
    df = df[df[padj_col] > 0]
    df = (df.sort_values(lfc_col, key=abs, ascending=False)
            .drop_duplicates(subset=id_col)
            .set_index(id_col))
    return df


def make_volcano(df, lfc_col, padj_col, title, out_name):
    d = df.copy()
    d["neglog10padj"] = -np.log10(d[padj_col].clip(lower=1e-300))
    d["category"] = "not significant"
    d.loc[(d[lfc_col] > LFC_THRESH) & (d[padj_col] < PADJ_THRESH), "category"] = "up in PerDKO"
    d.loc[(d[lfc_col] < -LFC_THRESH) & (d[padj_col] < PADJ_THRESH), "category"] = "down in PerDKO"

    fig, ax = plt.subplots(figsize=(7, 5.5))
    for cat, sub in d.groupby("category"):
        ax.scatter(sub[lfc_col], sub["neglog10padj"], s=8, c=COLORS[cat],
                   label=f"{cat} (n={len(sub)})", alpha=0.7, edgecolors="none")
    ax.axvline(LFC_THRESH, color="grey", lw=0.5, ls="--")
    ax.axvline(-LFC_THRESH, color="grey", lw=0.5, ls="--")
    ax.axhline(-np.log10(PADJ_THRESH), color="grey", lw=0.5, ls="--")
    ax.set_xlabel("log2FoldChange (PerDKO vs WT)")
    ax.set_ylabel("-log10(padj)")
    lfc_min, lfc_max = d[lfc_col].min(), d[lfc_col].max()
    ax.set_title(f"{title}\nlog2FC range: {lfc_min:.2f} to {lfc_max:.2f}", fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, fontsize=8)
    plt.tight_layout()
    path = os.path.join(OUT_DIR, out_name)
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"  Saved volcano: {path}")


def make_scatter(illumina_df, illumina_lfc_col, illumina_padj_col, nano_df, title, out_name):
    merged = illumina_df.join(nano_df, how="inner", lsuffix="_ill", rsuffix="_nano")
    if len(merged) == 0:
        print(f"  WARNING: no overlap for {title}, skipping scatter.")
        return

    ill_sig = (merged[illumina_padj_col] < PADJ_THRESH) & (merged[illumina_lfc_col].abs() > LFC_THRESH)
    nano_sig = (merged["padj_per_ko"] < PADJ_THRESH) & (merged["log2FC_perko"].abs() > LFC_THRESH)

    merged["category"] = "not significant"
    merged.loc[ill_sig & ~nano_sig, "category"] = "Illumina only"
    merged.loc[~ill_sig & nano_sig, "category"] = "Nanopore only"
    merged.loc[ill_sig & nano_sig, "category"] = "Both significant"

    colors = {"not significant": "lightgrey", "Illumina only": "#56B4E9",
              "Nanopore only": "#009E73", "Both significant": "#D55E00"}

    fig, ax = plt.subplots(figsize=(7, 6.5))
    for cat in ["not significant", "Illumina only", "Nanopore only", "Both significant"]:
        sub = merged[merged["category"] == cat]
        if len(sub) == 0:
            continue
        z = 1 if cat == "Both significant" else 0
        ax.scatter(sub[illumina_lfc_col], sub["log2FC_perko"], s=10, c=colors[cat],
                   label=f"{cat} (n={len(sub)})", alpha=0.7,
                   edgecolors="black" if cat == "Both significant" else "none",
                   linewidths=0.3, zorder=z)
    lim = max(merged[illumina_lfc_col].abs().max(), merged["log2FC_perko"].abs().max()) * 1.05
    ax.plot([-lim, lim], [-lim, lim], color="grey", lw=0.5, ls="--", zorder=-1)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("log2FC -- Illumina (PerDKO / WT)")
    ax.set_ylabel("log2FC -- Nanopore (PerDKO / WT)")
    ax.set_title(f"{title}\n(n={len(merged)} shared genes)", fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    plt.tight_layout()
    path = os.path.join(OUT_DIR, out_name)
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"  Saved scatter: {path}")


print("Loading Nanopore data...")
nano_raw = pd.read_excel("data/Table3-differential_expressiong_WTvsPerKO.xlsx", sheet_name="DE_all")
nano_gene = gene_level_dedup(nano_raw, "log2FC_perko", "padj_per_ko", "gene_id")
nano_gene_by_symbol = gene_level_dedup(nano_raw, "log2FC_perko", "padj_per_ko", "mgi_symbol")

print("\n=== GEO official counts ===")
geo = pd.read_csv("results/GSE130613/marsseq_sanity_check/deseq2_results_geo_counts.csv", index_col=0)
geo = geo.dropna(subset=["log2FoldChange", "padj"])
make_scatter(geo, "log2FoldChange", "padj", nano_gene_by_symbol,
             "GEO official counts vs Nanopore", "geo_vs_nanopore_scatter.png")
# volcano already exists: results/GSE130613/marsseq_sanity_check/sanity_check_volcano_geo_counts.png

print("\n=== 3'UTR-extended Salmon ===")
ext = pd.read_csv("results/GSE130613/normalisation_3utr_extended/deseq2_KO_vs_WT.tsv", sep="\t")
ext_gene = gene_level_dedup(ext, "log2FoldChange", "padj_gene_stageR", "gene_id")
make_volcano(ext_gene, "log2FoldChange", "padj_gene_stageR",
             "3'UTR-extended Salmon pipeline", "3utr_extended_volcano.png")
make_scatter(ext_gene, "log2FoldChange", "padj_gene_stageR", nano_gene,
             "3'UTR-extended Salmon vs Nanopore", "3utr_extended_vs_nanopore_scatter.png")

print("\n=== STAR+htseq (extended GTF) ===")
star_ext = pd.read_csv("results/GSE130613/normalisation_star_htseq/deseq2_star_htseq_KO_vs_WT.tsv", sep="\t")
star_ext = star_ext.dropna(subset=["log2FoldChange", "padj"]).set_index("gene_id")
make_volcano(star_ext, "log2FoldChange", "padj",
             "STAR + htseq-count (extended GTF)", "star_htseq_extended_volcano.png")
make_scatter(star_ext, "log2FoldChange", "padj", nano_gene,
             "STAR+htseq (extended GTF) vs Nanopore", "star_htseq_extended_vs_nanopore_scatter.png")

print("\n=== STAR+htseq (original GTF) ===")
star_orig = pd.read_csv("results/GSE130613/normalisation_star_htseq_original_gtf/deseq2_star_htseq_KO_vs_WT.tsv", sep="\t")
star_orig = star_orig.dropna(subset=["log2FoldChange", "padj"]).set_index("gene_id")
make_volcano(star_orig, "log2FoldChange", "padj",
             "STAR + htseq-count (original GTF)", "star_htseq_original_volcano.png")
make_scatter(star_orig, "log2FoldChange", "padj", nano_gene,
             "STAR+htseq (original GTF) vs Nanopore", "star_htseq_original_vs_nanopore_scatter.png")

print("\nAll figures saved to:", OUT_DIR)
