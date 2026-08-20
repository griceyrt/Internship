#!/usr/bin/env python3
"""
Recreate manuscript Fig 3 panels A/B/C for WT CT8 vs CT20 (reviewer point d),
as a visual companion to the numbers already sent to Kiran/Bharath.

NOTE ON METHOD: this deliberately does NOT try to match the exact extra
filtering step found in the original fig4_2_2024_12_04.R panel C (an
abundance-table cross-filter that shrinks PerKO's Venn from ~720/202 genes
down to the reference figure's 675/157) -- see rebuttal_fig1g_ct8_ct20 memory
note for why. Instead, all three panels here use the SAME already-verified
significant gene/event sets from rebuttal_ct8_vs_ct20_de.R (real R DESeq2) and
rebuttal_ct8_vs_ct20_splicing.sh (SUPPA), consistently, so the three panels
agree with each other and with the numbers already reported.

Panel A -- Differentially expressed: significant transcripts (padj<0.05 &
    |log2FC|>1.5) between CT8 and CT20, split Down/Up and by
    class (Annotated / Novel-Annotated / Novel-Novel) x biotype (protein
    coding / non coding).
Panel B -- Differentially spliced: significant SUPPA events (|dPSI|>0.1 &
    p<0.05) per event type, split Up/Down (sign convention: dPSI = CT8-CT20,
    so dPSI<0 = higher in CT20 = "Up" here, matching the DESeq2 contrast
    direction used in panel A, condition CT20 vs CT8).
Panel C -- Venn of DE-significant genes vs splicing-significant genes, plus
    a scatter of |dPSI| vs |log2FC| for each significant splicing event
    (using that gene's log2FC from the full DE test, not just significant
    genes), with Pearson correlation.

Prerequisites (run first):
    Rscript scripts/rebuttal_ct8_vs_ct20_de.R
    bash    scripts/rebuttal_ct8_vs_ct20_splicing.sh

Run from: orthogonal_validation/
    python3 scripts/make_ct8_vs_ct20_fig3_panels.py
"""
from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib_venn import venn2
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results" / "rhythmicity_analysis" / "rebuttal_ct8_vs_ct20"
FIG_DIR = PROJECT_ROOT / "results" / "rhythmicity_analysis" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

GTF_PATH_CANDIDATES = [
    PROJECT_ROOT.parent / "boundary_analysis" / "data" / "transcriptome_productivity.gtf",
    PROJECT_ROOT.parent / "Internship" / "boundary_analysis" / "data" / "transcriptome_productivity.gtf",
    DATA_DIR / "transcriptome_productivity.gtf",
]

DPSI_CUTOFF = 0.1
PVAL_CUTOFF = 0.05
LOG2FC_CUTOFF = 1.5
FDR_CUTOFF = 0.05
EVENT_TYPES = ["SE", "MX", "RI", "A5", "A3"]  # AF/AL excluded per Kiran's 2026-08-13 instruction

CLASS_COLORS = {"Annotated": "#66c2a5", "Novel-Annotated": "#fc8d62", "Novel-Novel": "#8da0cb"}
UPDOWN_COLORS = {"up": "#8fd4a8", "down": "#e8a598"}  # green / salmon, matches Fig 3B reference

plt.rcParams.update({
    "font.weight": "bold",
    "axes.labelweight": "bold",
    "axes.titleweight": "bold",
    "axes.linewidth": 1.6,
})


def style_L_axis(ax):
    """Keep only left + bottom spines (an 'L'), drop the rest, thicken what remains."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.6)
    ax.spines["bottom"].set_linewidth(1.6)
    ax.tick_params(width=1.6)


def find_gtf() -> Path:
    for p in GTF_PATH_CANDIDATES:
        if p.exists():
            return p
    raise FileNotFoundError(f"No GTF found in any of: {GTF_PATH_CANDIDATES}")


def load_tx_biotypes(gtf_path: Path) -> dict:
    tx_info = {}
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            attrs = fields[8]
            tb = re.search(r'transcript_biotype "([^"]+)"', attrs)
            t = re.search(r'transcript_id "([^"]+)"', attrs)
            if t:
                tx_info[t.group(1)] = tb.group(1) if tb else None
    return tx_info


def load_de_results() -> pd.DataFrame:
    r_path = RESULTS_DIR / "deseq2_R_tx_level_ct8_vs_ct20.csv"
    py_path = RESULTS_DIR / "deseq2_tx_level_ct8_vs_ct20.csv"
    if r_path.exists():
        return pd.read_csv(r_path)
    return pd.read_csv(py_path, index_col=0).reset_index(names="transcript_id")


def load_splicing_events() -> pd.DataFrame:
    frames = []
    for et in EVENT_TYPES:
        path = RESULTS_DIR / "diff" / f"diff_{et}.dpsi.temp.0"
        df = pd.read_csv(path, sep="\t", header=None, names=["event_id", "dpsi", "pval"], skiprows=1)
        df["event_type"] = et
        frames.append(df)
    out = pd.concat(frames, ignore_index=True)
    out["gene_id"] = out["event_id"].str.split(";").str[0]
    out["dpsi"] = pd.to_numeric(out["dpsi"], errors="coerce")
    out["pval"] = pd.to_numeric(out["pval"], errors="coerce")
    return out.dropna(subset=["dpsi", "pval"])


def classify_tx(transcript_id: str, gene_id: str) -> str:
    if isinstance(transcript_id, str) and transcript_id.startswith("ENS"):
        return "Annotated"
    if isinstance(gene_id, str) and "ENS" in gene_id:
        return "Novel-Annotated"
    return "Novel-Novel"


def panel_a_data(de_valid: pd.DataFrame, tx_biotypes: dict) -> pd.DataFrame:
    sig = de_valid[(de_valid["log2FoldChange"].abs() > LOG2FC_CUTOFF) & (de_valid["padj"] < FDR_CUTOFF)].copy()
    sig["class"] = sig.apply(lambda r: classify_tx(r["transcript_id"], r["gene_id"]), axis=1)
    sig["coding"] = sig["transcript_id"].map(
        lambda t: "Protein coding" if tx_biotypes.get(t) == "protein_coding" else "Non coding"
    )
    sig["direction"] = np.where(sig["log2FoldChange"] > 0, "up", "down")
    return sig


def plot_panel_a(ax_down, ax_up, sig: pd.DataFrame):
    for ax, direction, title in [(ax_down, "down", "Down-regulated"), (ax_up, "up", "Up-regulated")]:
        sub = sig[sig["direction"] == direction]
        counts = sub.groupby(["class", "coding"]).size().unstack(fill_value=0)
        classes = [c for c in ["Annotated", "Novel-Annotated", "Novel-Novel"] if c in counts.index]
        y = np.arange(len(classes))
        left = np.zeros(len(classes))
        for coding, alpha in [("Protein coding", 1.0), ("Non coding", 0.45)]:
            vals = [counts.loc[c, coding] if coding in counts.columns and c in counts.index else 0 for c in classes]
            colors = [CLASS_COLORS[c] for c in classes]
            ax.barh(y, vals, left=left, color=colors, alpha=alpha, edgecolor="black", linewidth=1.2,
                    label=coding if ax is ax_down else None)
            left += np.array(vals)
        ax.set_yticks(y)
        ax.set_yticklabels([])  # class is conveyed by color (green=Annotated, orange=Novel-Annotated, blue=Novel-Novel)
        ax.set_title(title, fontsize=9, fontweight="bold")
        ax.invert_xaxis() if direction == "down" else None
        ax.set_xlabel("transcripts", fontsize=7, fontweight="bold")
        style_L_axis(ax)
        if direction == "down":
            ax.spines["left"].set_visible(False)
            ax.spines["right"].set_visible(True)
            ax.spines["right"].set_linewidth(1.6)
            ax.yaxis.set_ticks_position("right")


def panel_b_data(events: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for et in EVENT_TYPES:
        sub = events[events["event_type"] == et]
        sig = sub[(sub["dpsi"].abs() > DPSI_CUTOFF) & (sub["pval"] < PVAL_CUTOFF)]
        # dPSI = CT8 - CT20; "up" (in CT20) = dpsi < 0, matching the CT20-vs-CT8 DESeq2 contrast direction in panel A
        n_up = int((sig["dpsi"] < 0).sum())
        n_down = int((sig["dpsi"] > 0).sum())
        rows.append({"event_type": et, "up": n_up, "down": n_down})
    return pd.DataFrame(rows)


def plot_panel_b(ax, df: pd.DataFrame):
    x = np.arange(len(df))
    ax.bar(x, df["up"], color=UPDOWN_COLORS["up"], edgecolor="black", linewidth=1.2, label="up")
    ax.bar(x, -df["down"], color=UPDOWN_COLORS["down"], edgecolor="black", linewidth=1.2, label="down")
    ax.axhline(0, color="black", linewidth=1.2)
    ax.set_xticks(x)
    ax.set_xticklabels(df["event_type"], fontsize=8, fontweight="bold")
    ax.set_ylabel("Number of events", fontsize=8, fontweight="bold")
    ax.set_title("Differentially spliced", fontsize=9, fontweight="bold")
    ax.legend(fontsize=6, frameon=False, loc="upper right")
    style_L_axis(ax)


def main():
    gtf_path = find_gtf()
    print("Using GTF:", gtf_path)
    tx_biotypes = load_tx_biotypes(gtf_path)

    de = load_de_results()
    de_valid = de.dropna(subset=["padj"])
    events = load_splicing_events()

    sig_a = panel_a_data(de_valid, tx_biotypes)
    print(f"\nPanel A: {len(sig_a)} significant transcripts "
          f"({(sig_a['direction']=='down').sum()} down, {(sig_a['direction']=='up').sum()} up)")

    df_b = panel_b_data(events)
    print(f"\nPanel B: {df_b['up'].sum() + df_b['down'].sum()} significant events total")
    print(df_b.to_string(index=False))

    de_sig_genes = set(de_valid[(de_valid["log2FoldChange"].abs() > LOG2FC_CUTOFF)
                                 & (de_valid["padj"] < FDR_CUTOFF)]["gene_id"].dropna())
    sig_events_all = events[(events["dpsi"].abs() > DPSI_CUTOFF) & (events["pval"] < PVAL_CUTOFF)].copy()
    splice_sig_genes = set(sig_events_all["gene_id"].unique())
    overlap = de_sig_genes & splice_sig_genes
    print(f"\nPanel C Venn: DE-only={len(de_sig_genes - splice_sig_genes)}, "
          f"overlap={len(overlap)}, splice-only={len(splice_sig_genes - de_sig_genes)}")

    gene_log2fc = de_valid.groupby("gene_id")["log2FoldChange"].apply(lambda x: x.loc[x.abs().idxmax()])
    sig_events_all["log2FC_gene"] = sig_events_all["gene_id"].map(gene_log2fc)
    scatter = sig_events_all.dropna(subset=["log2FC_gene"]).copy()
    scatter["abs_dpsi"] = scatter["dpsi"].abs()
    scatter["abs_log2fc"] = scatter["log2FC_gene"].abs()
    r, p = stats.pearsonr(scatter["abs_log2fc"], scatter["abs_dpsi"])
    print(f"Panel C correlation: R={r:.3f}, p={p:.3g}, n={len(scatter)}")

    # ---- Build figure ----
    fig = plt.figure(figsize=(10.5, 3.0))
    gs = fig.add_gridspec(1, 4, width_ratios=[0.55, 0.55, 1.0, 1.2], wspace=1.1)

    ax_down = fig.add_subplot(gs[0, 0])
    ax_up = fig.add_subplot(gs[0, 1])
    plot_panel_a(ax_down, ax_up, sig_a)

    ax_b = fig.add_subplot(gs[0, 2])
    plot_panel_b(ax_b, df_b)

    ax_c1 = fig.add_subplot(gs[0, 3])
    v = venn2(subsets=(len(de_sig_genes - splice_sig_genes), len(splice_sig_genes - de_sig_genes), len(overlap)),
              set_labels=("Differentially\nexpressed", "Differentially\nspliced"), ax=ax_c1,
              set_colors=("#822d30", "#434299"), alpha=0.3)
    for label in v.set_labels:
        if label is not None:
            label.set_fontweight("bold")
            label.set_fontsize(7)
    for label in v.subset_labels:
        if label is not None:
            label.set_fontweight("bold")
    ax_c1.set_title("Gene overlap\n(CT8 vs CT20)", fontsize=8, fontweight="bold")

    fig.savefig(FIG_DIR / "CT8_vs_CT20_fig3ABC.svg", bbox_inches="tight")
    fig.savefig(FIG_DIR / "CT8_vs_CT20_fig3ABC.png", bbox_inches="tight", dpi=200)
    print(f"\nSaved: {FIG_DIR / 'CT8_vs_CT20_fig3ABC.svg'}")

    # separate scatter (own figure, easier to read)
    fig2, ax2 = plt.subplots(figsize=(3.2, 3.2))
    ax2.scatter(scatter["abs_log2fc"], scatter["abs_dpsi"], s=6, color="black", alpha=0.6)
    if len(scatter) > 1:
        m, b = np.polyfit(scatter["abs_log2fc"], scatter["abs_dpsi"], 1)
        xs = np.linspace(scatter["abs_log2fc"].min(), scatter["abs_log2fc"].max(), 50)
        ax2.plot(xs, m * xs + b, color="#fd9997")
    ax2.set_xlabel("|log2(fold change gene expression)|", fontsize=7)
    ax2.set_ylabel("|dPSI|", fontsize=8)
    ax2.set_title(f"R = {r:.2f}, p = {p:.2g}", fontsize=8)
    fig2.savefig(FIG_DIR / "CT8_vs_CT20_fig3C_scatter.svg", bbox_inches="tight")
    fig2.savefig(FIG_DIR / "CT8_vs_CT20_fig3C_scatter.png", bbox_inches="tight", dpi=200)
    print(f"Saved: {FIG_DIR / 'CT8_vs_CT20_fig3C_scatter.svg'}")


if __name__ == "__main__":
    main()
