"""
report_figures.py
Recreates the figures for the UCBL internship report:
    fig1_central_dogma.png
    fig2_as_event_types.png
    fig3_psi_metric.png
    fig4_pipeline_overview.png
    fig5_sequencing_comparison.png

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
# FIGURE 4 — Pipeline Overview
# Add this function to report_figures.py, then add
# fig_pipeline("figures/plots/fig4_pipeline_overview.png")
# to the if __name__ == "__main__": block
# ═══════════════════════════════════════════════════════════════════════════════

def fig_pipeline(outpath="figures/plots/fig4_pipeline_overview.png"):

    def rounded_box(ax, x, y, w, h, facecolor, edgecolor, label, sublabel="", fontsize=8):
        box = FancyBboxPatch((x, y), w, h,
                             boxstyle="round,pad=0.012",
                             facecolor=facecolor, edgecolor=edgecolor,
                             linewidth=1.2, zorder=3)
        ax.add_patch(box)
        cx, cy = x + w/2, y + h/2
        if sublabel:
            ax.text(cx, cy + 0.022, label, ha="center", va="center",
                    fontsize=fontsize, fontweight="bold", color=edgecolor, zorder=4)
            ax.text(cx, cy - 0.022, sublabel, ha="center", va="center",
                    fontsize=fontsize - 1, color=GRAY, zorder=4)
        else:
            ax.text(cx, cy, label, ha="center", va="center",
                    fontsize=fontsize, fontweight="bold", color=edgecolor, zorder=4)

    def arr(ax, x1, x2, y, color=GRAY, lw=1.2):
        ax.annotate("", xy=(x2, y), xytext=(x1, y),
                    arrowprops=dict(arrowstyle="-|>", color=color,
                                    lw=lw, mutation_scale=10), zorder=3)

    fig, ax = plt.subplots(figsize=(13, 5.5))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    BW = 0.11
    BH = 0.13
    GAP = 0.015

    # Row labels
    ax.text(0.01, 0.78, "Project 1", ha="left", va="center",
            fontsize=9, fontweight="bold", color=BLUE, rotation=90)
    ax.text(0.01, 0.35, "Project 2", ha="left", va="center",
            fontsize=9, fontweight="bold", color=TEAL, rotation=90)

    # Track 1 — Nanopore
    track1_y = 0.62
    xs = [0.05, 0.21, 0.37, 0.53, 0.69, 0.85]
    labels1 = [
        ("Mouse liver\nsamples",    "CT8 / CT20 / PerDKO CT20"),
        ("Nanopore\nRNA-seq",       "Yutaka lab, Tokyo"),
        ("Salmon\nTPM",             "custom GTF 155k tx"),
        ("SUPPA2\ngenerateEvents",  "strict vs variable 1nt"),
        ("PSI comparison\n+ diffSplice", "WT CT20 vs PerDKO CT20"),
        ("IGV\nvalidation",         "candidate events"),
    ]
    for i, (x, (lbl, sub)) in enumerate(zip(xs, labels1)):
        rounded_box(ax, x, track1_y, BW, BH, BLUE_L, BLUE, lbl, sub, fontsize=7.5)
        if i < len(xs) - 1:
            arr(ax, x + BW + GAP, xs[i+1] - GAP, track1_y + BH/2, color=BLUE)

    # Track 2 — Illumina
    track2_y = 0.18
    xs2 = [0.05, 0.21, 0.37, 0.53, 0.69, 0.85]
    labels2 = [
        ("GEO GSE171975",           "Aviram et al. 2021"),
        ("Illumina\nRNA-seq",       "4x WT + 3x PerDKO CT20"),
        ("Salmon\nquantification",  "custom GTF 155k tx"),
        ("DESeq2 +\nbatch correction", "tximport + limma"),
        ("SUPPA2\npsiPerEvent",     "diffSplice WT vs PerDKO"),
        ("Overlap\nanalysis",       "vs Project 1 findings"),
    ]
    for i, (x, (lbl, sub)) in enumerate(zip(xs2, labels2)):
        rounded_box(ax, x, track2_y, BW, BH, TEAL_L, TEAL, lbl, sub, fontsize=7.5)
        if i < len(xs2) - 1:
            arr(ax, x + BW + GAP, xs2[i+1] - GAP, track2_y + BH/2, color=TEAL)

    # Vertical connector between last boxes
    cx_last = xs[-1] + BW/2
    ax.annotate("", xy=(cx_last, track1_y - 0.01),
                xytext=(cx_last, track2_y + BH),
                arrowprops=dict(arrowstyle="-|>", color=PURP,
                                lw=1.2, mutation_scale=9), zorder=3)
    ax.text(cx_last + 0.015, (track1_y + track2_y + BH) / 2,
            "orthogonal\nvalidation", ha="left", va="center",
            fontsize=7, color=PURP, style="italic")

    # Right bracket
    bx = 0.975
    ax.annotate("", xy=(bx, track1_y + BH/2), xytext=(bx, track2_y + BH/2),
                arrowprops=dict(arrowstyle="-", color=PURP, lw=1.5,
                                connectionstyle="bar,fraction=0.0"), zorder=3)
    ax.text(bx + 0.01, (track1_y + track2_y + BH) / 2,
            "Orthogonal\nValidation", ha="left", va="center",
            fontsize=8, fontweight="bold", color=PURP)

    # Shared GTF note
    mid_y = (track1_y + track2_y + BH) / 2
    ax.text(0.37 + BW/2, mid_y + 0.14,
            "same reference transcriptome (custom GTF)", ha="center", va="center",
            fontsize=7, color=GRAY, style="italic",
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                      edgecolor=GRAY_L, linewidth=0.6))

    # Title
    ax.text(0.5, 0.97, "Computational pipeline overview",
            ha="center", va="top", fontsize=11, fontweight="bold", color=GRAY)

    fig.tight_layout(pad=0.3)
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {outpath}")

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 5 — Short-read vs Long-read sequencing comparison
# Paste this function into report_figures.py after fig_pipeline()
# Then add this line to the if __name__ == "__main__": block:
#   fig_sequencing_comparison("figures/plots/fig5_sequencing_comparison.png")
# ═══════════════════════════════════════════════════════════════════════════════

def fig_sequencing_comparison(outpath="figures/plots/fig5_sequencing_comparison.png"):

    def draw_rna_molecule(ax, y=0.86):
        exon_colors = [BLUE, TEAL, BLUE, TEAL, BLUE]
        positions = [0.05, 0.18, 0.33, 0.46, 0.61, 0.74, 0.89]
        widths    = [0.10, 0.12, 0.10, 0.12, 0.10, 0.12, 0.08]
        ax.plot([0.05, 0.97], [y, y], color=GRAY_L, lw=0.5, zorder=2)
        for i, (x, w) in enumerate(zip(positions, widths)):
            color = exon_colors[i % len(exon_colors)] if i % 2 == 0 else GRAY_L
            box = FancyBboxPatch((x, y-0.015), w, 0.03,
                                 boxstyle="round,pad=0.003",
                                 facecolor=color, edgecolor="none", zorder=3)
            ax.add_patch(box)

    def arr(ax, x, y1, y2, color=GRAY):
        ax.annotate("", xy=(x, y2), xytext=(x, y1),
                    arrowprops=dict(arrowstyle="-|>", color=color,
                                    lw=1.2, mutation_scale=10), zorder=4)

    WHITE = "#FFFFFF"
    fig, axes = plt.subplots(1, 2, figsize=(13, 8))
    fig.patch.set_facecolor(WHITE)

    # ── LEFT: Illumina ────────────────────────────────────────────
    ax = axes[0]
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.add_patch(FancyBboxPatch((0, 0), 1, 1, boxstyle="round,pad=0.01",
                                facecolor="#FDF6EE", edgecolor=AMBER_L, linewidth=1.5))

    ax.text(0.5, 0.97, "Illumina", ha="center", va="top",
            fontsize=13, fontweight="bold", color=AMBER)
    ax.text(0.5, 0.93, "short-read sequencing (~150 nt)", ha="center", va="top",
            fontsize=8.5, color=GRAY)
    draw_rna_molecule(ax, y=0.86)
    ax.text(0.5, 0.82, "original RNA molecule", ha="center",
            fontsize=7.5, color=GRAY, style="italic")

    arr(ax, 0.5, 0.80, 0.77, AMBER)
    ax.text(0.5, 0.76, "1. Fragmentation into short pieces", ha="center",
            fontsize=8, color=GRAY, style="italic")

    # scattered fragments
    frag_xs = [0.05,0.14,0.23,0.33,0.43,0.52,0.60,0.69,0.77,0.85,0.42,0.62,0.22,0.55]
    frag_ys = [0.70,0.66,0.70,0.66,0.70,0.66,0.70,0.66,0.70,0.66,0.62,0.62,0.62,0.74]
    frag_ws = [0.08,0.08,0.08,0.08,0.08,0.07,0.08,0.07,0.07,0.08,0.08,0.08,0.08,0.08]
    frag_cs = [AMBER_L,BLUE_L,TEAL_L,AMBER_L,BLUE_L,TEAL_L,AMBER_L,BLUE_L,
               TEAL_L,AMBER_L,BLUE_L,TEAL_L,AMBER_L,BLUE_L]
    for x,y,w,c in zip(frag_xs,frag_ys,frag_ws,frag_cs):
        ax.add_patch(FancyBboxPatch((x,y-0.012),w,0.022,boxstyle="round,pad=0.002",
                                    facecolor=c,edgecolor="white",linewidth=0.5,zorder=3))
    ax.text(0.5, 0.58, "millions of short fragments", ha="center",
            fontsize=7.5, color=GRAY, style="italic")

    arr(ax, 0.5, 0.57, 0.54, AMBER)
    ax.text(0.5, 0.53, "2. Computational alignment and assembly", ha="center",
            fontsize=8, color=GRAY, style="italic")

    # stacked short reads with gaps
    read_rows = [
        [(0.05,0.09),(0.17,0.09),(0.32,0.09),(0.47,0.09),(0.62,0.09),(0.78,0.09)],
        [(0.08,0.09),(0.22,0.09),(0.40,0.09),(0.55,0.09),(0.70,0.09),(0.85,0.09)],
        [(0.06,0.09),(0.26,0.09),(0.50,0.09),(0.74,0.09),(0.88,0.07)],
        [(0.10,0.09),(0.38,0.09),(0.60,0.09),(0.80,0.09)],
    ]
    for row, ry in zip(read_rows, [0.47,0.44,0.41,0.38]):
        for (rx,rw) in row:
            ax.add_patch(FancyBboxPatch((rx,ry-0.01),rw,0.018,
                                        boxstyle="round,pad=0.002",
                                        facecolor=AMBER_L,edgecolor=AMBER,
                                        linewidth=0.4,zorder=3))

    # gap annotations
    for gx, gl in [(0.22, "gap?"), (0.66, "gap?")]:
        ax.annotate("", xy=(gx+0.06, 0.34), xytext=(gx, 0.34),
                    arrowprops=dict(arrowstyle="<->", color=CORAL, lw=1.2))
        ax.text(gx+0.03, 0.32, gl, ha="center", fontsize=7, color=CORAL)

    arr(ax, 0.5, 0.31, 0.28, AMBER)
    ax.text(0.5, 0.27, "3. Reconstructed sequence (ambiguous)", ha="center",
            fontsize=8, color=GRAY, style="italic")

    positions = [0.05,0.16,0.28,0.40,0.52,0.64,0.76,0.87]
    widths_r  = [0.09,0.08,0.09,0.09,0.09,0.09,0.08,0.09]
    for i,(x,w) in enumerate(zip(positions,widths_r)):
        if i in [2,5]:
            ax.add_patch(FancyBboxPatch((x,0.21),w,0.025,boxstyle="round,pad=0.002",
                                        facecolor=CORAL_L,edgecolor=CORAL,
                                        linewidth=0.8,linestyle="dashed",zorder=3))
            ax.text(x+w/2,0.2225,"?",ha="center",va="center",
                    fontsize=7,color=CORAL,fontweight="bold")
        else:
            ax.add_patch(FancyBboxPatch((x,0.21),w,0.025,boxstyle="round,pad=0.002",
                                        facecolor=AMBER_L if i%2==0 else TEAL_L,
                                        edgecolor="none",zorder=3))

    ax.add_patch(FancyBboxPatch((0.05,0.04),0.90,0.12,boxstyle="round,pad=0.01",
                                facecolor=CORAL_L,edgecolor=CORAL,linewidth=1.2))
    ax.text(0.5,0.135,"Ambiguous isoform reconstruction",ha="center",
            fontsize=8.5,color=CORAL,fontweight="bold")
    ax.text(0.5,0.085,"Short reads cannot distinguish\nbetween similar isoforms",
            ha="center",fontsize=7.5,color=GRAY)

    # ── RIGHT: Nanopore ───────────────────────────────────────────
    ax = axes[1]
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.add_patch(FancyBboxPatch((0,0),1,1,boxstyle="round,pad=0.01",
                                facecolor="#EEF7F3",edgecolor=TEAL_L,linewidth=1.5))

    ax.text(0.5,0.97,"Nanopore",ha="center",va="top",
            fontsize=13,fontweight="bold",color=TEAL)
    ax.text(0.5,0.93,"long-read sequencing (up to ~50,000 nt)",ha="center",va="top",
            fontsize=8.5,color=GRAY)
    draw_rna_molecule(ax, y=0.86)
    ax.text(0.5,0.82,"original RNA molecule",ha="center",
            fontsize=7.5,color=GRAY,style="italic")

    arr(ax, 0.5, 0.80, 0.77, TEAL)
    ax.text(0.5,0.76,"1. Direct reading of full-length molecules",ha="center",
            fontsize=8,color=GRAY,style="italic")

    long_reads = [
        (0.05,0.92,TEAL_L,TEAL),(0.05,0.92,BLUE_L,BLUE),
        (0.05,0.92,TEAL_L,TEAL),(0.05,0.70,BLUE_L,BLUE),
        (0.05,0.92,TEAL_L,TEAL),
    ]
    for (x,xe,fc,ec),ry in zip(long_reads,[0.72,0.68,0.64,0.60,0.56]):
        ax.add_patch(FancyBboxPatch((x,ry-0.012),xe-x,0.022,
                                    boxstyle="round,pad=0.002",
                                    facecolor=fc,edgecolor=ec,linewidth=0.6,zorder=3))

    ax.annotate("alternative\nisoform!", xy=(0.72,0.60),xytext=(0.80,0.67),
                fontsize=7,color=BLUE,
                arrowprops=dict(arrowstyle="->",color=BLUE,lw=0.8))
    ax.text(0.5,0.52,"Each read = one complete molecule",ha="center",
            fontsize=7.5,color=GRAY,style="italic")

    arr(ax, 0.5, 0.51, 0.48, TEAL)
    ax.text(0.5,0.47,"2. No assembly required",ha="center",
            fontsize=8,color=GRAY,style="italic")

    for x,w in zip(positions,widths_r):
        ax.add_patch(FancyBboxPatch((x,0.40),w,0.025,boxstyle="round,pad=0.002",
                                    facecolor=TEAL_L,edgecolor=TEAL,
                                    linewidth=0.5,zorder=3))

    arr(ax, 0.5, 0.39, 0.36, TEAL)
    ax.text(0.5,0.35,"3. Two isoforms clearly identified",ha="center",
            fontsize=8,color=GRAY,style="italic")

    for i,(ry,exon_set,label,col,edg) in enumerate([
        (0.29,[0,1,2,3,4,5,6,7],"known isoform",TEAL_L,TEAL),
        (0.22,[0,1,3,5,6,7],"novel isoform (skipped exon)",BLUE_L,BLUE),
    ]):
        for j,(x,w) in enumerate(zip(positions,widths_r)):
            if j in exon_set:
                ax.add_patch(FancyBboxPatch((x,ry-0.012),w,0.020,
                                            boxstyle="round,pad=0.002",
                                            facecolor=col,edgecolor=edg,
                                            linewidth=0.5,zorder=3))
        ax.text(0.5,ry-0.026,label,ha="center",fontsize=7,color=edg)

    ax.add_patch(FancyBboxPatch((0.05,0.04),0.90,0.12,boxstyle="round,pad=0.01",
                                facecolor=TEAL_L,edgecolor=TEAL,linewidth=1.2))
    ax.text(0.5,0.135,"Direct and unambiguous isoform identification",
            ha="center",fontsize=8.5,color=TEAL,fontweight="bold")
    ax.text(0.5,0.085,"Full-length reads reveal novel transcripts\nnot detectable with short reads",
            ha="center",fontsize=7.5,color=GRAY)

    fig.suptitle("Short-read vs long-read RNA sequencing: reconstruction challenge",
                 fontsize=11,fontweight="bold",color=GRAY,y=0.999)
    plt.tight_layout(rect=[0,0,1,0.998])
    fig.savefig(outpath, dpi=200, bbox_inches="tight", facecolor=WHITE)
    plt.close(fig)
    print(f"Saved: {outpath}")


# ═══════════════════════════════════════════════════════════════════════════════
# Outputs paths
if __name__ == "__main__":
    fig_central_dogma("figures/plots/fig1_central_dogma.png")
    fig_as_events("figures/plots/fig2_as_event_types.png")
    fig_psi("figures/plots/fig3_psi_metric.png")
    fig_pipeline("figures/plots/fig4_pipeline_overview.png")
    fig_sequencing_comparison("figures/plots/fig5_sequencing_comparison.png")