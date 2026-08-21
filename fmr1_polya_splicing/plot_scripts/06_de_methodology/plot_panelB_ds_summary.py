#!/usr/bin/env python3
# =============================================================================
# Panel B: Differentially spliced events, WT vs Fmr1 KO, up/down per SUPPA
# event type. Style reference: Chikhaoui/Mamgain et al. PERIOD manuscript
# Fig 4B.
# Author: Gricey
#
# Input: results/gse188840_suppa_2026-08-05/diff/diff_{EVENT}.dpsi.temp.0,
# all 7 event types.
#
# Significance: raw p < 0.05 (matches this project's SUPPA convention
# throughout, e.g. run_suppa.sh's own summary and Overlap #1).
#
# Direction: confirmed from run_suppa.sh's diffSplice call --
# `--psi psi/{EVENT}_WT.psi psi/{EVENT}_KO.psi`, i.e. PSI order is
# (WT, KO), so dPSI = PSI_KO - PSI_WT. Positive dpsi_value = increased in
# KO ("up"), negative = decreased in KO ("down") -- same KO-relative-to-WT
# direction as DESeq2's log2FoldChange, so Panel A and B are consistent.
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
DIFF_DIR = f"{BASE}/results/gse188840_suppa_2026-08-05/diff"
OUT_DIR = f"{BASE}/results/panelAB_summary"
os.makedirs(OUT_DIR, exist_ok=True)

EVENT_TYPES = ["SE", "MX", "RI", "A5", "A3", "AL", "AF"]  # order matches the reference figure

counts = {}
for event in EVENT_TYPES:
    path = f"{DIFF_DIR}/diff_{event}.dpsi.temp.0"
    df = pd.read_csv(path, sep="\t")
    df.columns = ["event_id", "dpsi_value", "pvalue"] + list(df.columns[3:])
    df["pvalue"] = pd.to_numeric(df["pvalue"], errors="coerce")
    df["dpsi_value"] = pd.to_numeric(df["dpsi_value"], errors="coerce")
    sig = df[df["pvalue"] < 0.05]
    n_up = int((sig["dpsi_value"] > 0).sum())
    n_down = int((sig["dpsi_value"] < 0).sum())
    counts[event] = {"up": n_up, "down": n_down}
    print(f"{event}: up={n_up}, down={n_down}, total={n_up + n_down}")

fig, ax = plt.subplots(figsize=(8, 5))
x = range(len(EVENT_TYPES))

up_vals = [counts[e]["up"] for e in EVENT_TYPES]
down_vals = [-counts[e]["down"] for e in EVENT_TYPES]

ax.bar(x, up_vals, color="lightskyblue", edgecolor="black", linewidth=0.6, label="up")
ax.bar(x, down_vals, color="coral", edgecolor="black", linewidth=0.6, label="down")

ax.axhline(0, color="black", linewidth=0.8)
ax.set_xticks(list(x))
ax.set_xticklabels(EVENT_TYPES)
ax.set_ylabel("Number of events")
ax.set_title("Differentially spliced events, WT vs Fmr1 KO\n(raw p < 0.05, per event type)",
              fontsize=12, fontweight="bold")
ax.legend(title=None, loc="upper right", frameon=False)
for spine in ("top", "right"):
    ax.spines[spine].set_visible(False)

fig.tight_layout()
png_path = f"{OUT_DIR}/panelB_ds_summary.png"
svg_path = f"{OUT_DIR}/panelB_ds_summary.svg"
fig.savefig(png_path, dpi=200, bbox_inches="tight")
fig.savefig(svg_path, bbox_inches="tight")
plt.close(fig)
print("Saved:", png_path)
print("Saved:", svg_path)
