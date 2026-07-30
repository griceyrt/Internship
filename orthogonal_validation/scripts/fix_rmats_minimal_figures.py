#!/usr/bin/env python3
"""
Fixes for the two rMATS minimal volcano SVGs only (p0.3 and FDR0.05).
Does NOT touch GSE130613_STARhtseq_DE_volcano_minimal.svg or
GSE130613_STARhtseq_vs_Nanopore_DE_concordance_minimal.svg -- Gricey is
hand-editing those in Inkscape already, this script never writes those
filenames.

Fixes applied vs the original make_final_figures_GSE130613_new_pipeline.py:
  1. Removed rasterized=True from the background scatter -- that was
     embedding the grey dots as a bitmap, which is what looked pixelated
     when zoomed in Inkscape. Now fully vector, like everything else.
  2. Added explicit, sparse tick marks (was relying on matplotlib's auto
     ticks, which came out denser than the minimal style elsewhere).

Run from: orthogonal_validation/
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams.update({
    "font.family":       "sans-serif",
    "font.sans-serif":   ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":         8,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "axes.spines.left":  True,
    "axes.spines.bottom":True,
    "axes.edgecolor":    "#AAAAAA",
    "axes.linewidth":    0.7,
    "axes.grid":         False,
    "figure.facecolor":  "white",
    "axes.facecolor":    "white",
    "xtick.color":       "#444444",
    "ytick.color":       "#444444",
    "xtick.labelsize":   7,
    "ytick.labelsize":   7,
    "axes.labelsize":    8,
    "axes.labelcolor":   "#222222",
    "xtick.major.size":  3,
    "ytick.major.size":  3,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
})

BASE_DIR  = os.environ.get("OV_BASE_DIR", "/Users/gricey/Desktop/Internship/orthogonal_validation")
RMATS_DIR = os.path.join(BASE_DIR, "results", "SRP194523_rmats_2026-07-20")
OUT = os.path.join(BASE_DIR, "results", "figures", "figures_final_new_pipeline")

DPSI_THRESH = 0.1
PVAL_THRESH = 0.3
FDR_THRESH  = 0.05

EVENT_TYPE_MAP = {"SE": "SE", "A5SS": "A5", "A3SS": "A3", "MXE": "MX", "RI": "RI"}
EVENT_COLORS = {"A3": "#E69F00", "A5": "#56B4E9", "MX": "#0072B2", "RI": "#D55E00", "SE": "#CC79A7"}

def save_svg(fig, name):
    path = os.path.join(OUT, name)
    fig.savefig(path, format="svg", bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {name}")

print("Loading rMATS JC output...")
dfs = []
for rmats_type, code in EVENT_TYPE_MAP.items():
    d = pd.read_csv(os.path.join(RMATS_DIR, f"{rmats_type}.MATS.JC.txt"), sep="\t")
    d["GeneID"] = d["GeneID"].str.strip('"')
    d["geneSymbol"] = d["geneSymbol"].str.strip('"')
    d["event_type"] = code
    dfs.append(d[["event_type", "GeneID", "geneSymbol", "PValue", "FDR", "IncLevelDifference"]])
rmats = pd.concat(dfs, ignore_index=True)
rmats["dPSI"] = -rmats["IncLevelDifference"]

# =============================================================================
# FIGURE A — rMATS dPSI volcano, p<0.3 (fixed: vector dots, sparse ticks)
# =============================================================================
print("\n=== Figure A: rMATS dPSI volcano, p<0.3 ===")
d = rmats[rmats["PValue"] > 0].copy()
d["neg_log10_p"] = -np.log10(d["PValue"])
d["sig"] = (d["PValue"] < PVAL_THRESH) & (d["dPSI"].abs() > DPSI_THRESH)

fig, ax = plt.subplots(figsize=(3.5, 3.5))
ns = d[~d["sig"]]
ax.scatter(ns["dPSI"], ns["neg_log10_p"], c="#CCCCCC", s=4, alpha=0.4, linewidths=0, zorder=1)
for etype, col in EVENT_COLORS.items():
    sub = d[d["sig"] & (d["event_type"] == etype)]
    if len(sub):
        ax.scatter(sub["dPSI"], sub["neg_log10_p"], c=col, s=14, alpha=0.9, linewidths=0, zorder=2)
ax.axvline(DPSI_THRESH, color="#BBBBBB", ls="--", lw=0.7)
ax.axvline(-DPSI_THRESH, color="#BBBBBB", ls="--", lw=0.7)
ax.axhline(-np.log10(PVAL_THRESH), color="#BBBBBB", ls="--", lw=0.7)
ax.set_xlabel("dPSI (PerDKO − WT)")
ax.set_ylabel("−log₁₀(p)")
ax.set_xticks([-1, -0.5, 0, 0.5, 1])
ax.set_yticks([0, 2, 4, 6, 8])
save_svg(fig, "GSE130613_rmats_dPSI_volcano_p0.3_minimal.svg")

# =============================================================================
# FIGURE B — rMATS dPSI volcano, FDR<0.05 (fixed: vector dots, sparse ticks)
# =============================================================================
print("\n=== Figure B: rMATS dPSI volcano, FDR<0.05 ===")
d2 = rmats[rmats["FDR"] > 0].copy()
d2["neg_log10_fdr"] = -np.log10(d2["FDR"])
d2["sig"] = (d2["FDR"] < FDR_THRESH) & (d2["dPSI"].abs() > DPSI_THRESH)

fig, ax = plt.subplots(figsize=(3.5, 3.5))
ns2 = d2[~d2["sig"]]
ax.scatter(ns2["dPSI"], ns2["neg_log10_fdr"], c="#CCCCCC", s=4, alpha=0.4, linewidths=0, zorder=1)
for etype, col in EVENT_COLORS.items():
    sub = d2[d2["sig"] & (d2["event_type"] == etype)]
    if len(sub):
        ax.scatter(sub["dPSI"], sub["neg_log10_fdr"], c=col, s=14, alpha=0.9, linewidths=0, zorder=2)
ax.axvline(DPSI_THRESH, color="#BBBBBB", ls="--", lw=0.7)
ax.axvline(-DPSI_THRESH, color="#BBBBBB", ls="--", lw=0.7)
ax.axhline(-np.log10(FDR_THRESH), color="#BBBBBB", ls="--", lw=0.7)
ax.set_xlabel("dPSI (PerDKO − WT)")
ax.set_ylabel("−log₁₀(FDR)")
ax.set_xticks([-1, -0.5, 0, 0.5, 1])
ax.set_yticks([0, 2, 4, 6])
save_svg(fig, "GSE130613_rmats_dPSI_volcano_FDR0.05_minimal.svg")

print("\n=== DONE (only the 2 rMATS files were touched) ===")
