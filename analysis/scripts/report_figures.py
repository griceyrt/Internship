"""
report_figures.py
Recreates the three figures for the UCBL internship report:
    fig1_central_dogma.png
    fig2_as_event_types.png
    fig3_psi_metric.png

Requirements: matplotlib, numpy
    pip install matplotlib numpy
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import matplotlib.patheffects as pe
import numpy as np

# ── Shared style ──────────────────────────────────────────────────────────────
BLUE   = "#185FA5"
BLUE_L = "#B5D4F4"
TEAL   = "#1D9E75"
TEAL_L = "#9FE1CB"
PURP   = "#534AB7"
PURP_L = "#CECBF6"
AMBER  = "#BA7517"
AMBER_L= "#FAC775"
CORAL  = "#993C1D"
CORAL_L= "#F5C4B3"
GRAY   = "#5F5E5A"
GRAY_L = "#D3D1C7"

FONT = "DejaVu Sans"
plt.rcParams.update({
    "font.family": FONT,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.spines.left": False,
    "axes.spines.bottom": False,
})

def exon_box(ax, x, y, w=0.08, h=0.06, color=BLUE_L, edge=BLUE,
             label="", fontsize=8):
    """Draw a single exon box."""
    box = FancyBboxPatch((x - w/2, y - h/2), w, h,
                         boxstyle="round,pad=0.005",
                         facecolor=color, edgecolor=edge, linewidth=0.8, zorder=3)
    ax.add_patch(box)
    if label:
        ax.text(x, y, label, ha="center", va="center",
                fontsize=fontsize, color=edge, fontweight="bold", zorder=4)

def intron_line(ax, x1, x2, y, color=GRAY, lw=1.0):
    ax.plot([x1, x2], [y, y], color=color, lw=lw, zorder=2)

def skip_arc(ax, x1, x2, y, color=GRAY, lw=1.0):
    mid = (x1 + x2) / 2
    arc_y = y + 0.055
    ax.annotate("", xy=(x2, y), xytext=(x1, y),
                arrowprops=dict(arrowstyle="-", color=color, lw=lw,
                                connectionstyle=f"arc3,rad=-0.5"),
                zorder=2)

def arrow(ax, x1, x2, y, color=GRAY, label="", fontsize=7):
    ax.annotate("", xy=(x2, y), xytext=(x1, y),
                arrowprops=dict(arrowstyle="-|>", color=color, lw=1.0,
                                mutation_scale=10),
                zorder=3)
    if label:
        ax.text((x1+x2)/2, y + 0.025, label, ha="center", va="bottom",
                fontsize=fontsize, color=GRAY)


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — Central Dogma
# ═══════════════════════════════════════════════════════════════════════════════
def fig_central_dogma(outpath="fig1_central_dogma.png"):
    fig, ax = plt.subplots(figsize=(10, 3.2))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # ── A: DNA ────────────────────────────────────────────────────────
    exon_box(ax, 0.09, 0.5, w=0.13, h=0.18, color=BLUE_L, edge=BLUE, label="A")
    ax.text(0.09, 0.34, "DNA\n(gene)", ha="center", va="top",
            fontsize=7.5, color=BLUE)

    # arrow: transcription
    arrow(ax, 0.16, 0.24, 0.5, color=GRAY, label="transcription")

    # ── B: pre-mRNA ───────────────────────────────────────────────────
    exon_box(ax, 0.30, 0.5, w=0.10, h=0.18, color=PURP_L, edge=PURP, label="B")
    ax.text(0.30, 0.34, "pre-mRNA", ha="center", va="top",
            fontsize=7.5, color=PURP)

    # splicing arrow (fans out to 3 isoforms)
    ax.annotate("", xy=(0.44, 0.74), xytext=(0.36, 0.54),
                arrowprops=dict(arrowstyle="-|>", color=GRAY, lw=0.9, mutation_scale=9))
    ax.annotate("", xy=(0.44, 0.50), xytext=(0.36, 0.50),
                arrowprops=dict(arrowstyle="-|>", color=GRAY, lw=0.9, mutation_scale=9))
    ax.annotate("", xy=(0.44, 0.26), xytext=(0.36, 0.46),
                arrowprops=dict(arrowstyle="-|>", color=GRAY, lw=0.9, mutation_scale=9))
    ax.text(0.37, 0.87, "splicing", ha="center", fontsize=7.5, color=GRAY)

    # ── C: 3 isoforms ─────────────────────────────────────────────────
    for i, (yc, label) in enumerate(zip([0.74, 0.50, 0.26],
                                        ["C1", "C2", "C3"])):
        exon_box(ax, 0.51, yc, w=0.12, h=0.15, color=TEAL_L, edge=TEAL, label=label)

    ax.text(0.51, 0.08, "mature mRNA\nisoforms", ha="center", va="top",
            fontsize=7.5, color=TEAL)

    # translation arrows
    for yc in [0.74, 0.50, 0.26]:
        arrow(ax, 0.58, 0.66, yc, color=GRAY)
    ax.text(0.62, 0.87, "translation", ha="center", fontsize=7.5, color=GRAY)

    # ── D: 3 proteins ─────────────────────────────────────────────────
    for i, (yc, label) in enumerate(zip([0.74, 0.50, 0.26],
                                        ["D1", "D2", "D3"])):
        exon_box(ax, 0.74, yc, w=0.12, h=0.15, color=CORAL_L, edge=CORAL, label=label)

    ax.text(0.74, 0.08, "proteins", ha="center", va="top",
            fontsize=7.5, color=CORAL)

    # Legend
    legend_text = (
        "A = DNA (gene)    B = pre-mRNA    "
        "C = mature mRNA isoform    D = protein\n"
        "One gene (A) can produce multiple isoforms (C1, C2, C3) via alternative splicing."
    )
    # Saving the legend for the word report since it doesn't fit well in the figure:
    # ax.text(0.5, 0.01, legend_text, ha="center", va="bottom",
    #        fontsize=6.8, color=GRAY, style="italic",
    #        transform=ax.transAxes)

    fig.tight_layout(pad=0.3)
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {outpath}")


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — AS Event Types
# ═══════════════════════════════════════════════════════════════════════════════
def draw_event_row(ax, y, event_name, event_desc,
                   color_l, color_e,
                   included_fn, excluded_fn):
    """Draw one AS event row: label | included isoform | excluded isoform."""
    # label
    ax.text(0.02, y, event_name, ha="left", va="center",
            fontsize=9, fontweight="bold", color=color_e)
    ax.text(0.02, y - 0.055, event_desc, ha="left", va="center",
            fontsize=7, color=GRAY)
    # isoforms
    included_fn(ax, 0.38, y, color_l, color_e)
    excluded_fn(ax, 0.73, y, color_l, color_e)


def se_included(ax, cx, y, cl, ce):
    for dx, lbl in [(-0.09, "E1"), (0, "E2"), (0.09, "E3")]:
        exon_box(ax, cx+dx, y, w=0.07, h=0.07, color=cl, edge=ce, label=lbl, fontsize=7)
    intron_line(ax, cx-0.055, cx-0.035, y, color=GRAY_L)
    intron_line(ax, cx+0.035, cx+0.055, y, color=GRAY_L)

def se_excluded(ax, cx, y, cl, ce):
    exon_box(ax, cx-0.07, y, w=0.07, h=0.07, color=cl, edge=ce, label="E1", fontsize=7)
    exon_box(ax, cx+0.07, y, w=0.07, h=0.07, color=cl, edge=ce, label="E3", fontsize=7)
    skip_arc(ax, cx-0.035, cx+0.035, y, color=GRAY)

def ri_included(ax, cx, y, cl, ce):
    exon_box(ax, cx-0.09, y, w=0.07, h=0.07, color=cl, edge=ce, label="E1", fontsize=7)
    # intron as thin rect
    rect = FancyBboxPatch((cx-0.045, y-0.015), 0.09, 0.03,
                          boxstyle="round,pad=0.003",
                          facecolor=GRAY_L, edgecolor=GRAY, linewidth=0.5, zorder=3)
    ax.add_patch(rect)
    ax.text(cx, y, "intron", ha="center", va="center", fontsize=6, color=GRAY)
    exon_box(ax, cx+0.09, y, w=0.07, h=0.07, color=cl, edge=ce, label="E2", fontsize=7)

def ri_excluded(ax, cx, y, cl, ce):
    exon_box(ax, cx-0.06, y, w=0.07, h=0.07, color=cl, edge=ce, label="E1", fontsize=7)
    exon_box(ax, cx+0.06, y, w=0.07, h=0.07, color=cl, edge=ce, label="E2", fontsize=7)
    intron_line(ax, cx-0.025, cx+0.025, y, color=GRAY_L)

def a3_included(ax, cx, y, cl, ce):
    exon_box(ax, cx-0.09, y, w=0.07, h=0.07, color=cl, edge=ce, label="E1", fontsize=7)
    exon_box(ax, cx+0.02, y, w=0.10, h=0.07, color=cl, edge=ce, label="E2 (long)", fontsize=6)
    intron_line(ax, cx-0.055, cx-0.03, y, color=GRAY_L)

def a3_excluded(ax, cx, y, cl, ce):
    exon_box(ax, cx-0.09, y, w=0.07, h=0.07, color=cl, edge=ce, label="E1", fontsize=7)
    exon_box(ax, cx+0.02, y, w=0.07, h=0.07, color=cl, edge=ce, label="E2 (short)", fontsize=5.5)
    intron_line(ax, cx-0.055, cx-0.015, y, color=GRAY_L)

def a5_included(ax, cx, y, cl, ce):
    exon_box(ax, cx-0.02, y, w=0.10, h=0.07, color=cl, edge=ce, label="E1 (long)", fontsize=6)
    exon_box(ax, cx+0.10, y, w=0.07, h=0.07, color=cl, edge=ce, label="E2", fontsize=7)
    intron_line(ax, cx+0.03, cx+0.065, y, color=GRAY_L)

def a5_excluded(ax, cx, y, cl, ce):
    exon_box(ax, cx-0.02, y, w=0.07, h=0.07, color=cl, edge=ce, label="E1 (short)", fontsize=5.5)
    exon_box(ax, cx+0.10, y, w=0.07, h=0.07, color=cl, edge=ce, label="E2", fontsize=7)
    intron_line(ax, cx+0.015, cx+0.065, y, color=GRAY_L)

def mx_included(ax, cx, y, cl, ce):
    exon_box(ax, cx-0.10, y, w=0.07, h=0.07, color=cl, edge=ce, label="E1", fontsize=7)
    exon_box(ax, cx,      y, w=0.07, h=0.07, color=cl, edge=ce, label="Ea", fontsize=7)
    exon_box(ax, cx+0.10, y, w=0.07, h=0.07, color=cl, edge=ce, label="E3", fontsize=7)
    intron_line(ax, cx-0.065, cx-0.035, y, color=GRAY_L)
    intron_line(ax, cx+0.035, cx+0.065, y, color=GRAY_L)

def mx_excluded(ax, cx, y, cl, ce):
    exon_box(ax, cx-0.10, y, w=0.07, h=0.07, color=cl, edge=ce, label="E1", fontsize=7)
    exon_box(ax, cx,      y, w=0.07, h=0.07, color=cl, edge=ce, label="Eb", fontsize=7)
    exon_box(ax, cx+0.10, y, w=0.07, h=0.07, color=cl, edge=ce, label="E3", fontsize=7)
    intron_line(ax, cx-0.065, cx-0.035, y, color=GRAY_L)
    intron_line(ax, cx+0.035, cx+0.065, y, color=GRAY_L)


def fig_as_events(outpath="fig2_as_event_types.png"):
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # Column headers
    ax.text(0.38, 0.97, "isoform: included", ha="center", va="top",
            fontsize=9, color=GRAY, fontweight="bold")
    ax.text(0.73, 0.97, "isoform: excluded", ha="center", va="top",
            fontsize=9, color=GRAY, fontweight="bold")
    ax.axhline(0.94, color=GRAY_L, lw=0.6)

    rows = [
        (0.84, "SE", "skipped exon",          BLUE_L,  BLUE,  se_included,  se_excluded),
        (0.67, "RI", "retained intron",        TEAL_L,  TEAL,  ri_included,  ri_excluded),
        (0.50, "A3", "alt. 3' splice site",    PURP_L,  PURP,  a3_included,  a3_excluded),
        (0.33, "A5", "alt. 5' splice site",    AMBER_L, AMBER, a5_included,  a5_excluded),
        (0.16, "MX", "mutually exclusive exons",CORAL_L, CORAL, mx_included, mx_excluded),
    ]

    for y, name, desc, cl, ce, inc_fn, exc_fn in rows:
        draw_event_row(ax, y, name, desc, cl, ce, inc_fn, exc_fn)
        ax.axhline(y - 0.10, color=GRAY_L, lw=0.4)

    # Legend
    # Saving the legend for the word report since it doesn't fit well in the figure:
    # ax.text(0.5, 0.02,
    #        "Colored boxes = exons.  Thin lines = introns.  Arc = skipped region.\n"
    #        "SUPPA2 detects and quantifies all five event types from RNA-seq data.",
    #        ha="center", va="bottom", fontsize=7, color=GRAY, style="italic")

    fig.tight_layout(pad=0.3)
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {outpath}")


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 3 — PSI Metric
# ═══════════════════════════════════════════════════════════════════════════════
def draw_psi_column(ax, cx, psi_val, n_included, n_total=4):
    """Draw n_total transcripts, n_included with exon, rest with arc."""
    y_start = 0.72
    dy = 0.13

    # title
    ax.text(cx, 0.93, f"PSI = {psi_val:.1f}", ha="center", va="center",
            fontsize=10, fontweight="bold",
            color=TEAL if psi_val == 1.0 else (BLUE if psi_val == 0.5 else GRAY))
    ax.text(cx, 0.83,
            ["event never\nincluded", "event in half\nof transcripts",
             "event always\nincluded"][int(psi_val * 2)],
            ha="center", va="center", fontsize=7.5, color=GRAY)

    for i in range(n_total):
        y = y_start - i * dy
        included = i < n_included

        if included:
            # 3 exons with middle one (event)
            for dx, lbl in [(-0.08, "E1"), (0, "E*"), (0.08, "E2")]:
                c_l = TEAL_L if psi_val == 1.0 else BLUE_L
                c_e = TEAL if psi_val == 1.0 else BLUE
                exon_box(ax, cx + dx, y, w=0.065, h=0.07,
                         color=c_l, edge=c_e, label=lbl, fontsize=6.5)
            intron_line(ax, cx-0.047, cx-0.033, y, color=GRAY_L)
            intron_line(ax, cx+0.033, cx+0.047, y, color=GRAY_L)
        else:
            # 2 flanking exons + arc (event skipped)
            exon_box(ax, cx-0.08, y, w=0.065, h=0.07,
                     color=GRAY_L, edge=GRAY, label="E1", fontsize=6.5)
            exon_box(ax, cx+0.08, y, w=0.065, h=0.07,
                     color=GRAY_L, edge=GRAY, label="E2", fontsize=6.5)
            skip_arc(ax, cx-0.047, cx+0.047, y, color=GRAY, lw=0.8)

    # count label
    ax.text(cx, 0.13,
            f"{n_included} / {n_total} transcripts\ninclude event",
            ha="center", va="center", fontsize=7.5, color=GRAY)


def fig_psi(outpath="fig3_psi_metric.png"):
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(0.5, 0.99, "PSI (percent spliced in)", ha="center", va="top",
            fontsize=12, fontweight="bold", color=GRAY)

    # vertical dividers
    ax.axvline(0.33, color=GRAY_L, lw=0.5, ymin=0.05, ymax=0.92)
    ax.axvline(0.66, color=GRAY_L, lw=0.5, ymin=0.05, ymax=0.92)

    draw_psi_column(ax, 0.165, 0.0, 0)
    draw_psi_column(ax, 0.495, 0.5, 2)
    draw_psi_column(ax, 0.825, 1.0, 4)

    # Saving the legend for the word report since it doesn't fit well in the figure:
    # ax.text(0.5, 0.06,
    #        "PSI = (transcripts including event) / (all transcripts at that locus)\n"
    #        "A change in PSI (dPSI) between conditions signals differential alternative splicing.",
    #        ha="center", va="center", fontsize=7.5, color=GRAY, style="italic")

    fig.tight_layout(pad=0.3)
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {outpath}")


# ═══════════════════════════════════════════════════════════════════════════════
# Outputs paths
if __name__ == "__main__":
    fig_central_dogma("figures/plots/fig1_central_dogma.png")
    fig_as_events("figures/plots/fig2_as_event_types.png")
    fig_psi("figures/plots/fig3_psi_metric.png")
