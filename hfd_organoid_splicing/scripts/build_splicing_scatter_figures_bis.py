#!/usr/bin/env python3
"""
PSI scatter figures, "bis" minimal style (HFD organoid splicing project)
Author: Gricey

Same data/logic as build_splicing_scatter_figures.py, styled closer to
Chikhaoui_Mamgain_MS.pdf Fig 3D: only 3 axis ticks, a subsampled/lighter
background cloud (the real dataset has far more tested events than the
paper's, so matching their visual sparseness means subsampling, not just
shrinking points), simplified spines. Alternate version, doesn't replace
the main script's output.

Run locally, not on the cluster -- lightweight plotting.
USAGE: python3 scripts/build_splicing_scatter_figures_bis.py
Run from: hfd_organoid_splicing/
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")

BASE = "."
SPLIT_DIR = f"{BASE}/results/suppa/split_by_condition"
SIG_TABLE = f"{BASE}/results/suppa/significant/significant_events.tsv"
GENE_SYMBOLS = f"{BASE}/data/reference/gene_id_to_symbol.tsv"
FIGURES_DIR = f"{BASE}/figures"

EVENT_TYPES = ["SE", "MX", "RI", "A5", "A3", "AL", "AF"]
COMPARISONS = ["W8_WHFD", "W18_WHFD", "W24_WHFD", "W42_WHFD"]
PANEL_LAYOUT = [
    ("SE", 0, 0), ("MX", 0, 1), ("RI", 0, 2), ("A5", 0, 3),
    ("A3", 1, 0), ("AL", 1, 1), ("AF", 1, 2),
]

N_LABELS_PER_PANEL = 3
MAX_BACKGROUND_POINTS = 1500
RNG_SEED = 42

GOI_EXACT = {
    "Clock", "Arntl",
    "Per1", "Per2", "Per3",
    "Cry1", "Cry2",
    "Dbp",
    "Hnf1b",
    "Ndufa6",
    "Sun1", "Sun2",
}
GOI_PREFIXES = ("Ddx", "Srsf", "Cpeb", "Insig")

def is_gene_of_interest(symbol):
    if symbol in GOI_EXACT:
        return True
    return any(symbol.startswith(p) for p in GOI_PREFIXES)

gene_symbols = pd.read_csv(GENE_SYMBOLS, sep="\t").set_index("gene_id")["gene_name"].to_dict()

def symbol_for(gene_id):
    return gene_symbols.get(gene_id, gene_id)

sig_df = pd.read_csv(SIG_TABLE, sep="\t")
rng = np.random.default_rng(RNG_SEED)

def mean_psi(event_type, condition):
    path = f"{SPLIT_DIR}/events_{event_type}_strict_{condition}.tab"
    df = pd.read_csv(path, sep="\t", header=0, index_col=0)
    df.columns = [f"rep{i+1}" for i in range(df.shape[1])]
    return df.mean(axis=1).rename("mean_psi")

DPSI_IS_TIMEPOINT_MINUS_WT = True  # confirmed empirically in the main script

for comparison in COMPARISONS:
    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    fig.suptitle(f"{comparison} vs WT_SD — PSI comparison per event type",
                 fontsize=16, fontweight="bold")
    axes[1, 3].axis("off")

    for event_type, row, col in PANEL_LAYOUT:
        ax = axes[row, col]

        wt_mean = mean_psi(event_type, "WT_SD")
        tp_mean = mean_psi(event_type, comparison)
        merged = pd.concat([wt_mean.rename("wt"), tp_mean.rename("tp")], axis=1).dropna()

        comp_sig = sig_df[(sig_df["comparison"] == f"{comparison}_vs_WT_SD") &
                           (sig_df["event_type"] == event_type)]
        sig_lookup = comp_sig.set_index("event_id")[["dPSI", "p_value"]]

        merged = merged.join(sig_lookup, how="left")
        if not DPSI_IS_TIMEPOINT_MINUS_WT:
            merged["dPSI"] = -merged["dPSI"]

        is_sig = merged["dPSI"].notna()
        is_up = is_sig & (merged["dPSI"] > 0)
        is_down = is_sig & (merged["dPSI"] < 0)
        is_ns = ~is_sig

        # subsample not-significant background; keep all significant points
        ns_index = merged.index[is_ns]
        if len(ns_index) > MAX_BACKGROUND_POINTS:
            ns_index = rng.choice(ns_index, size=MAX_BACKGROUND_POINTS, replace=False)

        ax.scatter(merged.loc[ns_index, "wt"], merged.loc[ns_index, "tp"],
                   c="gainsboro", s=8, alpha=0.5, linewidths=0)
        ax.scatter(merged.loc[is_down, "wt"], merged.loc[is_down, "tp"],
                   c="firebrick", s=22, alpha=0.85, linewidths=0)
        ax.scatter(merged.loc[is_up, "wt"], merged.loc[is_up, "tp"],
                   c="forestgreen", s=22, alpha=0.85, linewidths=0)
        ax.plot([0, 1], [0, 1], color="black", linewidth=0.5, linestyle="--")

        gene_ids = pd.Series(merged.index, index=merged.index).str.split(";").str[0]
        gene_syms = gene_ids.map(symbol_for)
        is_goi = gene_syms.map(is_gene_of_interest) & is_sig
        if is_goi.any():
            ax.scatter(merged.loc[is_goi, "wt"], merged.loc[is_goi, "tp"],
                       facecolors="none", edgecolors="black", s=60,
                       linewidths=1.4, zorder=5)

        top_events = set(comp_sig.sort_values("p_value").head(N_LABELS_PER_PANEL)["event_id"])
        goi_events = set(merged.index[is_goi])
        for eid in top_events | goi_events:
            if eid in merged.index:
                x, y = merged.loc[eid, "wt"], merged.loc[eid, "tp"]
                gene_id = eid.split(";")[0]
                ax.annotate(symbol_for(gene_id), (x, y), fontsize=7,
                            fontweight="bold", xytext=(3, 3),
                            textcoords="offset points")

        # minimal styling: 3 ticks, no top/right spine, thicker remaining spines
        ax.set_title(event_type, fontsize=13, fontweight="bold")
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        ax.set_xticks([0, 0.5, 1])
        ax.set_yticks([0, 0.5, 1])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_linewidth(1.8)
        ax.spines["left"].set_linewidth(1.8)
        ax.tick_params(axis="both", labelsize=13, width=1.8, length=6)
        ax.set_xlabel("PSI (WT_SD)", fontsize=10, fontweight="bold")
        ax.set_ylabel(f"PSI ({comparison})", fontsize=10, fontweight="bold")

    handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="forestgreen", markersize=8, label="Significant up"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="firebrick", markersize=8, label="Significant down"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="lightgrey", markersize=8, label="Not significant"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="none", markeredgecolor="black", markeredgewidth=1.2, markersize=9, label="Circadian/splicing-factor gene of interest"),
    ]
    fig.legend(handles=handles, loc="lower right", bbox_to_anchor=(0.98, 0.08), fontsize=11)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    png_path = f"{FIGURES_DIR}/splicing_scatter_{comparison}_vs_WT_SD_bis.png"
    svg_path = f"{FIGURES_DIR}/splicing_scatter_{comparison}_vs_WT_SD_bis.svg"
    plt.savefig(png_path, dpi=150)
    plt.savefig(svg_path)
    plt.close(fig)
    print(f"Saved: {png_path}")
    print(f"Saved: {svg_path}")

print("\nDone.")
