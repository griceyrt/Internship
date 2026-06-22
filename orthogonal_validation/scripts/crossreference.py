#!/usr/bin/env python3
"""
Orthogonal Validation — Steps 12 & 14
Author: Gricey

Step 12: Cross-reference Illumina vs Nanopore significant splicing events
Step 14: Venn diagram — DE genes vs Illumina differentially spliced genes

Significance thresholds:
  - Illumina splicing : |dPSI| > 0.1, p < 0.3
  - Nanopore splicing : |dPSI| > 0.1, p < 0.05
  - DE genes          : |log2FC| > 0.5, padj < 0.05

Run from: orthogonal_validation/
"""

import os
import glob
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# =============================================================================
# PATHS
# =============================================================================
BASE_DIR      = "/Users/gricey/Desktop/Internship/orthogonal_validation"
SUPPA_DATE    = "2026-06-17"
ILLUMINA_DIR  = os.path.join(BASE_DIR, "results", f"suppa_{SUPPA_DATE}", "diff")
NANOPORE_DIR  = os.path.join(BASE_DIR, "data", "nanopore_suppa")
OUT_DIR       = os.path.join(BASE_DIR, "results", "figures")
TABLES_DIR    = os.path.join(BASE_DIR, "results", "tables")
DESEQ_FILE    = os.path.join(BASE_DIR, "results", "normalisation", "deseq2_KO_vs_WT.tsv")
GTF_PATH      = "/Users/gricey/Desktop/Internship/boundary_analysis/data/transcriptome_productivity.gtf"
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(TABLES_DIR, exist_ok=True)

# =============================================================================
# THRESHOLDS
# =============================================================================
DPSI_THRESH        = 0.1
PVAL_ILLUMINA      = 0.3    # relaxed: Illumina has low power at n=3
PVAL_NANOPORE      = 0.05
LFC_THRESH         = 0.5
PADJ_THRESH        = 0.05

# =============================================================================
# LOAD SIGNIFICANT EVENTS
# =============================================================================

def build_gene_map(gtf_path):
    """Extract gene_id -> gene_name from GTF."""
    import re
    gene_map = {}
    with open(gtf_path) as fh:
        for line in fh:
            if "\tgene\t" not in line:
                continue
            gid   = re.search(r'gene_id "([^"]+)"', line)
            gname = re.search(r'gene_name "([^"]+)"', line)
            if gid and gname:
                gene_map[gid.group(1)] = gname.group(1)
    return gene_map


def load_significant_events(dpsi_dir, file_pattern, dpsi_thresh, pval_thresh):
    """Load all dpsi files and return (all_events set, significant_events set,
       sig_rows dict {event_id: {dPSI, pval, event_type}})."""
    all_events = []
    sig_events = []
    sig_rows   = {}
    for f in sorted(glob.glob(os.path.join(dpsi_dir, file_pattern))):
        df = pd.read_csv(f, sep="\t", index_col=0)
        df.columns = ["dPSI", "pval"]
        df = df.dropna(subset=["dPSI", "pval"])
        all_events.extend(df.index.tolist())
        sig = df[(abs(df["dPSI"]) > dpsi_thresh) & (df["pval"] < pval_thresh)]
        sig_events.extend(sig.index.tolist())
        for idx, row in sig.iterrows():
            sig_rows[idx] = {"dPSI": row["dPSI"], "pval": row["pval"]}
    return set(all_events), set(sig_events), sig_rows


print("Building gene name map from GTF...")
gene_map = build_gene_map(GTF_PATH)
print(f"  {len(gene_map)} genes loaded.")

print("Loading Illumina significant events...")
illumina_all, illumina_sig, illumina_rows = load_significant_events(
    ILLUMINA_DIR, "diff_*_strict.dpsi.temp.0", DPSI_THRESH, PVAL_ILLUMINA)

print("Loading Nanopore significant events...")
nanopore_all, nanopore_sig, nanopore_rows = load_significant_events(
    NANOPORE_DIR, "res_*.dpsi.temp.0", DPSI_THRESH, PVAL_NANOPORE)

# =============================================================================
# CROSS-REFERENCE
# =============================================================================
overlap       = illumina_sig & nanopore_sig
illumina_only = illumina_sig - nanopore_sig
nanopore_only = nanopore_sig - illumina_sig

print(f"\n{'='*50}")
print(f"Illumina significant events (p<{PVAL_ILLUMINA}) : {len(illumina_sig)}")
print(f"Nanopore significant events (p<{PVAL_NANOPORE}) : {len(nanopore_sig)}")
print(f"Overlap (both)                                   : {len(overlap)}")
print(f"Illumina only                                    : {len(illumina_only)}")
print(f"Nanopore only                                    : {len(nanopore_only)}")
print(f"{'='*50}")

if len(overlap) > 0:
    print("\nOverlapping event IDs:")
    for e in sorted(overlap):
        print(f"  {e}")
else:
    print("\nNo exact event ID overlap found.")

# =============================================================================
# STEP 12 — Excel: overlapping events
# =============================================================================
if len(overlap) > 0:
    rows = []
    for event_id in sorted(overlap):
        gene_id    = event_id.split(";")[0]
        event_type = event_id.split(";")[1].split(":")[0]
        rows.append({
            "gene_id":       gene_id,
            "gene_name":     gene_map.get(gene_id, "unannotated"),
            "event_type":    event_type,
            "event_id":      event_id,
            "dPSI_illumina": round(-illumina_rows[event_id]["dPSI"], 4),  # negate: KO-WT
            "pval_illumina": round( illumina_rows[event_id]["pval"], 4),
            "dPSI_nanopore": round( nanopore_rows[event_id]["dPSI"], 4),
            "pval_nanopore": round( nanopore_rows[event_id]["pval"], 4),
        })
    overlap_df = pd.DataFrame(rows)
    out = os.path.join(TABLES_DIR, "step12_overlap_events.xlsx")
    overlap_df.to_excel(out, index=False)
    print(f"\nStep 12 Excel saved: {out}")

# =============================================================================
# STEP 12 — Venn Diagram: Illumina vs Nanopore splicing
# =============================================================================
try:
    fig, ax = plt.subplots(figsize=(6, 5))
    venn2(
        subsets=(len(illumina_only), len(nanopore_only), len(overlap)),
        set_labels=(f"Illumina\n(n={len(illumina_sig)})",
                    f"Nanopore\n(n={len(nanopore_sig)})"),
        set_colors=("#56B4E9", "#E69F00"),
        alpha=0.6,
        ax=ax
    )
    ax.set_title(
        "Significant splicing events overlap\n"
        f"|dPSI| > {DPSI_THRESH} — Illumina p<{PVAL_ILLUMINA}, Nanopore p<{PVAL_NANOPORE}",
        fontsize=11
    )
    plt.tight_layout()
    out_path = os.path.join(OUT_DIR, "step12_overlap_venn.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Venn diagram saved: {out_path}")
except ImportError:
    print("matplotlib_venn not installed — pip install matplotlib-venn --break-system-packages")

# =============================================================================
# STEP 14 — Venn Diagram: DE genes vs Illumina spliced genes + Excel exports
# =============================================================================
print("\n=== Step 14: DE genes vs Spliced genes ===")

de = pd.read_csv(DESEQ_FILE, sep="\t", index_col=0)
de = de.dropna(subset=["log2FoldChange", "padj"])
de_full = de.copy()
de_sig  = set(de[(abs(de["log2FoldChange"]) > LFC_THRESH) & (de["padj"] < PADJ_THRESH)].index)

# Map event IDs to gene IDs (one row per unique gene)
illumina_sig_gene_events = {}   # gene_id -> list of (event_id, dPSI, pval)
for event_id, row in illumina_rows.items():
    gene_id = event_id.split(";")[0]
    illumina_sig_gene_events.setdefault(gene_id, []).append({
        "event_id":   event_id,
        "event_type": event_id.split(";")[1].split(":")[0],
        "dPSI":       round(-row["dPSI"], 4),   # negate: KO-WT
        "pval":       round( row["pval"], 4),
    })
illumina_sig_genes = set(illumina_sig_gene_events.keys())

de_splice_overlap = de_sig & illumina_sig_genes
de_only           = de_sig - illumina_sig_genes
splice_only       = illumina_sig_genes - de_sig

print(f"Significant DE genes       : {len(de_sig)}")
print(f"Significant spliced genes  : {len(illumina_sig_genes)}")
print(f"Both DE and spliced        : {len(de_splice_overlap)}")
print(f"DE only                    : {len(de_only)}")
print(f"Spliced only               : {len(splice_only)}")


def make_splice_rows(gene_set, label):
    """Build rows for Excel from a set of gene IDs (splicing side)."""
    rows = []
    for gene_id in sorted(gene_set):
        for ev in illumina_sig_gene_events.get(gene_id, []):
            row = {
                "gene_id":    gene_id,
                "gene_name":  gene_map.get(gene_id, "unannotated"),
                "event_type": ev["event_type"],
                "event_id":   ev["event_id"],
                "dPSI_illumina": ev["dPSI"],
                "pval_illumina": ev["pval"],
            }
            if gene_id in de_full.index:
                row["log2FoldChange"] = round(de_full.loc[gene_id, "log2FoldChange"], 4)
                row["padj"]           = round(de_full.loc[gene_id, "padj"], 6)
            rows.append(row)
    return rows


def make_de_only_rows(gene_set):
    """Build rows for Excel from DE-only genes (no splicing info)."""
    rows = []
    for gene_id in sorted(gene_set):
        rows.append({
            "gene_id":        gene_id,
            "gene_name":      gene_map.get(gene_id, "unannotated"),
            "log2FoldChange": round(de_full.loc[gene_id, "log2FoldChange"], 4),
            "padj":           round(de_full.loc[gene_id, "padj"], 6),
        })
    return rows


# Excel: DE only — sort named genes first, unannotated to bottom
de_only_df = pd.DataFrame(make_de_only_rows(de_only))
de_only_df = pd.concat([
    de_only_df[de_only_df["gene_name"] != "unannotated"].sort_values("gene_name"),
    de_only_df[de_only_df["gene_name"] == "unannotated"].sort_values("gene_id")
]).reset_index(drop=True)
de_only_df.to_excel(
    os.path.join(TABLES_DIR, "step14_DE_only.xlsx"), index=False)
print(f"Saved: step14_DE_only.xlsx  (n={len(de_only)})")

# Excel: overlap — both DE and spliced (18)
pd.DataFrame(make_splice_rows(de_splice_overlap, "overlap")).to_excel(
    os.path.join(TABLES_DIR, "step14_DE_and_spliced.xlsx"), index=False)
print(f"Saved: step14_DE_and_spliced.xlsx  (n={len(de_splice_overlap)})")

# Excel: spliced only (439)
pd.DataFrame(make_splice_rows(splice_only, "splice_only")).to_excel(
    os.path.join(TABLES_DIR, "step14_spliced_only.xlsx"), index=False)
print(f"Saved: step14_spliced_only.xlsx  (n={len(splice_only)})")

# Venn diagram
try:
    fig, ax = plt.subplots(figsize=(6, 5))
    venn2(
        subsets=(len(de_only), len(splice_only), len(de_splice_overlap)),
        set_labels=(
            f"DE genes\n(n={len(de_sig)})\n|log2FC|>{LFC_THRESH}, padj<{PADJ_THRESH}",
            f"Spliced genes\n(n={len(illumina_sig_genes)})\n|dPSI|>{DPSI_THRESH}, p<{PVAL_ILLUMINA}"
        ),
        set_colors=("#D55E00", "#56B4E9"),
        alpha=0.6,
        ax=ax
    )
    ax.set_title(
        "Genes with differential expression\nAND differential splicing (Illumina)",
        fontsize=11
    )
    plt.tight_layout()
    out_path = os.path.join(OUT_DIR, "step14_DE_vs_splicing_venn.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Venn diagram saved: {out_path}")
except ImportError:
    print("matplotlib_venn not installed — pip install matplotlib-venn --break-system-packages")

print("\n=== ALL DONE ===")
