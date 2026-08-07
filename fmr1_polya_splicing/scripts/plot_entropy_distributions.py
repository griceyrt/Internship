#!/usr/bin/env python3
# =============================================================================
# Entropy distribution comparison, WT vs Fmr1 KO -- Kiran's stats request
# Author: Gricey
#
# Kiran wants to see the per-gene Shannon entropy of isoform usage (already
# computed in isoform_proportions_step3.R -- entropy_WT and entropy_KO
# columns) plotted as two overlaid distributions, to check by eye whether
# the KO distribution is shifted relative to WT (isoform usage globally
# becoming more/less diverse under FMRP loss) or whether the two bell
# curves stay centered on each other.
#
# Statistical backing: entropy_WT and entropy_KO are PAIRED per gene (same
# gene contributes one value to each), so the right test for "is there a
# systematic shift" is a paired test, not an unpaired one. Entropy values
# are bounded and not expected to be normally distributed, so the
# nonparametric paired option -- Wilcoxon signed-rank -- is used instead
# of a paired t-test.
#
# Run once per TPM threshold (10/20/50), matching the existing family of
# outputs, since Gricey is still comparing all three even though Kiran
# leans toward the most stringent (50).
#
# Inputs: results/isoform_proportions/entropy_shift_ranked_min{10,20,50}TPM.tsv
# Outputs (per threshold): entropy_distribution_min{T}TPM.png/.svg
# Output (combined): entropy_distribution_stats_summary.tsv
#
# This only needs the small TSV files, not the cluster -- run it locally.
# Requires: pandas, matplotlib, scipy (pip install pandas matplotlib scipy)
# =============================================================================

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon

matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"]

# Input tables (entropy_shift_ranked_min*TPM.tsv) are still produced by the
# R pipeline (step1-4) and shared with plot_isoform_stacked_bars.py, so they
# stay put in isoform_proportions/. Outputs of THIS script get their own
# folder, entropy_distribution/, since it's a separate deliverable from the
# stacked bars (2026-08-06 reorg, per Gricey's request to keep the two
# figure types in clearly separate folders).
BASE = "/Users/gricey/Desktop/Internship/fmr1_polya_splicing/results/isoform_proportions"
OUT_DIR = "/Users/gricey/Desktop/Internship/fmr1_polya_splicing/results/entropy_distribution"

# Deliberately a DIFFERENT palette family from the stacked bars (Gricey's
# call, 2026-08-06): the stacked bars use hue to encode isoform identity
# within one gene, but this figure is a two-group comparison (WT vs KO
# across all genes), which is a different visual task -- red vs blue reads
# as "two conditions" more clearly than two shades of the same blue would.
# Modeled on a reference overlapping-histogram style (coral vs slate blue,
# semi-transparent so the overlap region blends to purple).
WT_COLOR = "#F08080"   # light coral
KO_COLOR = "#6A5ACD"   # slate blue


def analyze_threshold(threshold, bins=40):
    df = pd.read_csv(f"{BASE}/entropy_shift_ranked_min{threshold}TPM.tsv", sep="\t")
    wt = df["entropy_WT"].to_numpy()
    ko = df["entropy_KO"].to_numpy()

    # Paired Wilcoxon signed-rank: tests whether entropy_KO - entropy_WT
    # (the same entropy_diff column used for ranking) is systematically
    # non-zero across genes, i.e. whether the KO distribution is shifted
    # relative to WT rather than just individually noisy per gene.
    stat, pval = wilcoxon(wt, ko)
    mean_wt, mean_ko = float(np.mean(wt)), float(np.mean(ko))
    median_wt, median_ko = float(np.median(wt)), float(np.median(ko))

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.hist(wt, bins=bins, alpha=0.55, color=WT_COLOR,
            label=f"WT (median={median_wt:.2f})", density=True)
    ax.hist(ko, bins=bins, alpha=0.55, color=KO_COLOR,
            label=f"KO (median={median_ko:.2f})", density=True)
    ax.axvline(median_wt, color=WT_COLOR, linestyle="--", linewidth=1.2)
    ax.axvline(median_ko, color=KO_COLOR, linestyle="--", linewidth=1.2)

    ax.set_xlabel("Shannon entropy of isoform usage (bits)")
    ax.set_ylabel("Density")
    ax.set_title(f"Entropy distribution, WT vs KO (min TPM > {threshold}, n={len(df)} genes)",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=9, frameon=False, loc="upper right")

    if median_ko < median_wt:
        direction = "KO shifted lower (less diverse isoform usage on average)"
    elif median_ko > median_wt:
        direction = "KO shifted higher (more diverse isoform usage on average)"
    else:
        direction = "no median shift"

    # Stacked directly under the legend (also upper-right) rather than the
    # opposite corner, per Gricey's request -- keeps both annotations in
    # one visual cluster instead of pulling the eye to two corners.
    ax.text(
        0.98, 0.80,
        f"Wilcoxon signed-rank (paired):\nW={stat:.0f}, p={pval:.2e}\n{direction}",
        transform=ax.transAxes, fontsize=8, va="top", ha="right",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.85, edgecolor="0.7"),
    )

    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()

    png_path = f"{OUT_DIR}/entropy_distribution_min{threshold}TPM.png"
    svg_path = f"{OUT_DIR}/entropy_distribution_min{threshold}TPM.svg"
    fig.savefig(png_path, dpi=200, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)

    print(f"--- min{threshold}TPM (n={len(df)} genes) ---")
    print(f"  WT: mean={mean_wt:.4f}, median={median_wt:.4f}")
    print(f"  KO: mean={mean_ko:.4f}, median={median_ko:.4f}")
    print(f"  Wilcoxon signed-rank: W={stat:.1f}, p={pval:.4e}  ({direction})")
    print(f"  Saved: {png_path}")
    print(f"  Saved: {svg_path}\n")

    return {
        "threshold": threshold, "n_genes": len(df),
        "mean_WT": mean_wt, "mean_KO": mean_ko,
        "median_WT": median_wt, "median_KO": median_ko,
        "wilcoxon_W": stat, "wilcoxon_p": pval, "direction": direction,
    }


if __name__ == "__main__":
    results = [analyze_threshold(t) for t in (10, 20, 50)]
    summary = pd.DataFrame(results)
    summary_path = f"{OUT_DIR}/entropy_distribution_stats_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    print("Saved summary table:", summary_path)
