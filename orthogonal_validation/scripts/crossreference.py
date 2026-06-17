#!/usr/bin/env python3
"""
Orthogonal Validation — Step 12: Cross-reference
Author: Gricey

Compare significant splicing events from:
  - Illumina short-read (GSE130613, this project)
  - Nanopore long-read (manuscript, Khushi's analysis)

Match criterion: exact event ID (both datasets used the same GTF)
Significance thresholds:
  - Illumina:  |dPSI| > 0.1 AND p < 0.05
  - Nanopore:  |dPSI| > 0.1 AND p < 0.05

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
os.makedirs(OUT_DIR, exist_ok=True)

# Significance thresholds
DPSI_THRESH = 0.1
PVAL_THRESH = 0.05

# =============================================================================
# LOAD SIGNIFICANT EVENTS
# =============================================================================

def load_significant_events(dpsi_dir, file_pattern, dpsi_thresh, pval_thresh):
    """Load all dpsi files and return set of significant event IDs."""
    all_events = []
    sig_events = []

    for f in sorted(glob.glob(os.path.join(dpsi_dir, file_pattern))):
        df = pd.read_csv(f, sep="\t", index_col=0)
        df.columns = ["dPSI", "pval"]
        df = df.dropna(subset=["dPSI", "pval"])
        all_events.extend(df.index.tolist())
        sig = df[(abs(df["dPSI"]) > dpsi_thresh) & (df["pval"] < pval_thresh)]
        sig_events.extend(sig.index.tolist())

    return set(all_events), set(sig_events)


print("Loading Illumina significant events...")
illumina_all, illumina_sig = load_significant_events(
    ILLUMINA_DIR, "diff_*_strict.dpsi.temp.0", DPSI_THRESH, PVAL_THRESH)

print("Loading Nanopore significant events...")
nanopore_all, nanopore_sig = load_significant_events(
    NANOPORE_DIR, "res_*.dpsi.temp.0", DPSI_THRESH, PVAL_THRESH)

# =============================================================================
# CROSS-REFERENCE
# =============================================================================
overlap = illumina_sig & nanopore_sig
illumina_only = illumina_sig - nanopore_sig
nanopore_only = nanopore_sig - illumina_sig

print(f"\n{'='*50}")
print(f"Illumina significant events : {len(illumina_sig)}")
print(f"Nanopore significant events : {len(nanopore_sig)}")
print(f"Overlap (both)              : {len(overlap)}")
print(f"Illumina only               : {len(illumina_only)}")
print(f"Nanopore only               : {len(nanopore_only)}")
print(f"{'='*50}")

if len(overlap) > 0:
    print("\nOverlapping event IDs:")
    for e in sorted(overlap):
        print(f"  {e}")
else:
    print("\nNo exact event ID overlap found.")
    print("Consider gene-level matching as a fallback.")

# =============================================================================
# SAVE OVERLAP TABLE
# =============================================================================

# Build a summary table of overlapping events
if len(overlap) > 0:
    rows = []
    for dpsi_file in sorted(glob.glob(os.path.join(ILLUMINA_DIR, "diff_*_strict.dpsi.temp.0"))):
        df = pd.read_csv(dpsi_file, sep="\t", index_col=0)
        df.columns = ["dPSI_illumina", "pval_illumina"]
        for event_id in overlap:
            if event_id in df.index:
                rows.append({
                    "event_id": event_id,
                    "dPSI_illumina": -df.loc[event_id, "dPSI_illumina"],  # negate: KO-WT
                    "pval_illumina": df.loc[event_id, "pval_illumina"]
                })

    # Add Nanopore values
    nano_data = {}
    for dpsi_file in sorted(glob.glob(os.path.join(NANOPORE_DIR, "res_*.dpsi.temp.0"))):
        df = pd.read_csv(dpsi_file, sep="\t", index_col=0)
        df.columns = ["dPSI_nanopore", "pval_nanopore"]
        for event_id in overlap:
            if event_id in df.index:
                nano_data[event_id] = {
                    "dPSI_nanopore": df.loc[event_id, "dPSI_nanopore"],
                    "pval_nanopore": df.loc[event_id, "pval_nanopore"]
                }

    overlap_df = pd.DataFrame(rows).set_index("event_id")
    for col in ["dPSI_nanopore", "pval_nanopore"]:
        overlap_df[col] = overlap_df.index.map(lambda x: nano_data.get(x, {}).get(col, np.nan))

    overlap_out = os.path.join(BASE_DIR, "results", "normalisation", "overlap_events.tsv")
    overlap_df.to_csv(overlap_out, sep="\t")
    print(f"\nOverlap table saved: {overlap_out}")

# =============================================================================
# VENN DIAGRAM
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
    ax.set_title("Significant splicing events overlap\n|dPSI| > 0.1, p < 0.05", fontsize=11)
    plt.tight_layout()
    out_path = os.path.join(OUT_DIR, "step12_overlap_venn.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Venn diagram saved: {out_path}")
except ImportError:
    print("matplotlib_venn not installed — skipping Venn diagram.")
    print("Install with: pip install matplotlib-venn --break-system-packages")

print("\n=== DONE ===")
