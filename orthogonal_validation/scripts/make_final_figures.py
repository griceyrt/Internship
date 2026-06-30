#!/usr/bin/env python3
"""
Final publication-ready figures for Orthogonal Validation.
Output: results/figures/figures_final/
Format: SVG (vector), Arial font, no title, no legend, bare-minimum axes.
A companion legend.txt is generated explaining colors and thresholds.

Figures produced:
  Illumina_events_dPSI_volcano.svg
  Illumina_Nanopore_events_scatter.svg
  Illumina_DE_genes_volcano.svg
  Illumina_Nanopore_DE_genes_concordance.svg
  Illumina_Nanopore_DE_transcript_concordance.svg
"""

import os, glob, re, json, openpyxl
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from adjustText import adjust_text

# =============================================================================
# PATHS
# =============================================================================
BASE      = "/Users/gricey/Desktop/Internship/orthogonal_validation"
GTF       = "/Users/gricey/Desktop/Internship/boundary_analysis/data/transcriptome_productivity.gtf"
SUPPA_DIR = os.path.join(BASE, "results", "suppa_2026-06-20", "diff")
NAN_DIR   = os.path.join(BASE, "data", "nanopore_suppa")
DESEQ     = os.path.join(BASE, "results", "normalisation", "deseq2_KO_vs_WT.tsv")
KHUSHI    = os.path.join(BASE, "data", "Table3-differential_expressiong_WTvsPerKO.xlsx")
OUT       = os.path.join(BASE, "results", "figures", "figures_final")
os.makedirs(OUT, exist_ok=True)

# =============================================================================
# GLOBAL STYLE — Arial, no title, no legend, open axes, no grid
# =============================================================================
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

THRESHOLDS = {
    "dpsi":      0.1,
    "pval_ill":  0.3,
    "pval_nan":  0.05,
    "lfc":       0.5,
    "lfc_15":    1.5,
    "padj":      0.05,
}

# =============================================================================
# HELPERS
# =============================================================================
def build_gene_map(gtf):
    gmap = {}
    with open(gtf) as fh:
        for line in fh:
            if "\tgene\t" not in line: continue
            gid   = re.search(r'gene_id "([^"]+)"', line)
            gname = re.search(r'gene_name "([^"]+)"', line)
            if gid and gname:
                gmap[gid.group(1)] = gname.group(1)
    return gmap

def load_dpsi(dpsi_dir, pattern):
    dfs = []
    for f in sorted(glob.glob(os.path.join(dpsi_dir, pattern))):
        etype = os.path.basename(f).replace("diff_","").replace("_strict.dpsi.temp.0","").replace("res_","").replace(".dpsi.temp.0","")
        df = pd.read_csv(f, sep="\t", index_col=0, header=0)
        df.columns = ["dPSI","pval"]
        df["event_type"] = etype
        df = df.dropna().query("pval > 0")
        df["dPSI"] = -df["dPSI"]   # negate: file is WT-KO → KO-WT
        dfs.append(df)
    return pd.concat(dfs) if dfs else pd.DataFrame()

def save_svg(fig, name):
    path = os.path.join(OUT, name)
    fig.savefig(path, format="svg", bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {name}")

print("Building gene map...")
gene_map = build_gene_map(GTF)

# =============================================================================
# FIGURE 1 — Illumina_events_dPSI_volcano.svg
# =============================================================================
print("\n=== Figure 1: dPSI volcano ===")

EVENT_COLORS = {
    "A3": "#E69F00", "A5": "#56B4E9", "AF": "#009E73",
    "AL": "#F0E442", "MX": "#0072B2", "RI": "#D55E00", "SE": "#CC79A7",
}

dpsi = load_dpsi(SUPPA_DIR, "diff_*_strict.dpsi.temp.0")
dpsi["neg_log10_p"] = -np.log10(dpsi["pval"])
dpsi["sig"] = (dpsi["dPSI"].abs() > THRESHOLDS["dpsi"]) & (dpsi["pval"] < THRESHOLDS["pval_ill"])
print(f"  Significant events: {dpsi['sig'].sum()}")

fig, ax = plt.subplots(figsize=(3.5, 3.5))

# Background
ns = dpsi[~dpsi["sig"]]
ax.scatter(ns["dPSI"], ns["neg_log10_p"], c="#CCCCCC", s=4, alpha=0.4,
           linewidths=0, zorder=1, rasterized=True)

# Significant by event type
for etype, col in EVENT_COLORS.items():
    sub = dpsi[dpsi["sig"] & (dpsi["event_type"] == etype)]
    if len(sub):
        ax.scatter(sub["dPSI"], sub["neg_log10_p"], c=col, s=14, alpha=0.9,
                   linewidths=0, zorder=2, rasterized=True)

# Threshold lines
ax.axvline( THRESHOLDS["dpsi"], color="#BBBBBB", ls="--", lw=0.7)
ax.axvline(-THRESHOLDS["dpsi"], color="#BBBBBB", ls="--", lw=0.7)
ax.axhline(-np.log10(THRESHOLDS["pval_ill"]), color="#BBBBBB", ls="--", lw=0.7)

ax.set_xlabel("dPSI (PerDKO − WT)")
ax.set_ylabel("−log₁₀(p)")
ax.set_xticks([-1, -0.5, 0, 0.5, 1])
ax.set_yticks([0, 1, 2, 3])

save_svg(fig, "Illumina_events_dPSI_volcano.svg")

# =============================================================================
# FIGURE 2 — Illumina_Nanopore_events_scatter.svg
# =============================================================================
print("\n=== Figure 2: dPSI scatter ===")

ill = load_dpsi(SUPPA_DIR, "diff_*_strict.dpsi.temp.0")
nan = load_dpsi(NAN_DIR,   "res_*.dpsi.temp.0")

ill_sig = set(ill[(ill.dPSI.abs() > THRESHOLDS["dpsi"]) & (ill.pval < THRESHOLDS["pval_ill"])].index)
nan_sig = set(nan[(nan.dPSI.abs() > THRESHOLDS["dpsi"]) & (nan.pval < THRESHOLDS["pval_nan"])].index)
common  = set(ill.index) & set(nan.index)
overlap = ill_sig & nan_sig & common

fig, ax = plt.subplots(figsize=(3.5, 3.5))
ax.grid(False)

data = []
for eid in common:
    x = ill.loc[eid, "dPSI"]
    y = nan.loc[eid, "dPSI"]
    if eid in overlap:         col, sz, zo = "#E69F00", 30, 4
    elif eid in ill_sig:       col, sz, zo = "#56B4E9", 10, 2
    elif eid in nan_sig:       col, sz, zo = "#009E73", 10, 3
    else:                      col, sz, zo = "#CCCCCC",  3, 1
    data.append((x, y, col, sz, zo, eid))

x_arr = np.array([d[0] for d in data])
y_arr = np.array([d[1] for d in data])
for zo_val in [1, 2, 3, 4]:
    mask = np.array([d[4] == zo_val for d in data])
    c    = [d[2] for d in data]
    s    = [d[3] for d in data]
    ras  = zo_val == 1
    ax.scatter(x_arr[mask], y_arr[mask],
               c=np.array(c)[mask], s=np.array(s)[mask],
               linewidths=0.5 if zo_val == 4 else 0,
               edgecolors="black" if zo_val == 4 else "none",
               alpha=0.6 if zo_val == 1 else 1.0,
               zorder=zo_val, rasterized=ras)

# Labels on orange events
gene_best = {}
for x, y, col, sz, zo, eid in data:
    if col != "#E69F00": continue
    if max(abs(x), abs(y)) <= 0.2: continue
    gid   = eid.split(";")[0]
    gname = gene_map.get(gid, gid)
    if gname not in gene_best or abs(x) > abs(gene_best[gname][0]):
        gene_best[gname] = (x, y)

texts, gene_order = [], []
for gname, (dx, dy) in gene_best.items():
    ox = 0.04 if dx >= 0 else -0.04
    oy = 0.03 if dy >= 0 else -0.03
    texts.append(ax.text(dx+ox, dy+oy, gname, fontsize=6, fontstyle="italic"))
    gene_order.append((gname, dx, dy))

adjust_text(texts, ax=ax, expand=(3,3), force_text=(2,2), force_points=(1.5,1.5), lim=500)

for txt, (gname, dot_x, dot_y) in zip(texts, gene_order):
    lx, ly = txt.get_position(); txt.remove()
    ax.annotate(gname, xy=(dot_x,dot_y), xytext=(lx,ly),
                fontsize=6, fontstyle="italic", ha="center", va="center",
                arrowprops=dict(arrowstyle="-", color="#888888", lw=0.5, shrinkA=3, shrinkB=3))

lim = 1.05
ax.plot([-lim, lim], [-lim, lim], color="#BBBBBB", lw=0.7, ls="--", alpha=0.6)
ax.axhline(0, color="#BBBBBB", lw=0.5); ax.axvline(0, color="#BBBBBB", lw=0.5)
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
ax.set_xlabel("dPSI Illumina (PerDKO − WT)")
ax.set_ylabel("dPSI Nanopore (PerDKO − WT)")
ax.set_xticks([-1, -0.5, 0, 0.5, 1])
ax.set_yticks([-1, -0.5, 0, 0.5, 1])

save_svg(fig, "Illumina_Nanopore_events_scatter.svg")

# =============================================================================
# FIGURE 3 — Illumina_DE_genes_volcano.svg
# =============================================================================
print("\n=== Figure 3: DE volcano ===")

CIRCADIAN_GENES = [
    "Ren1","Orm2","Prtn3","Cyp2b9","Dbp","Cxcl1","Saa1","Ighg2b",
    "Celf4","Egr1","Ciart","Nr1d1","Nr1d2","Btg2","Per3","Ephx1",
    "Gstm2","Clock","Cyp2a5","Arntl","Cry1","Kmt2d","Gstm3",
    "Prkca","Cyp7a1","Fancm","Cspg5"
]

de_tx = pd.read_csv(DESEQ, sep="\t").dropna(subset=["log2FoldChange","padj_gene_stageR"])
de_tx = de_tx[de_tx["padj_gene_stageR"] > 0]
de = (de_tx.sort_values("log2FoldChange", key=abs, ascending=False)
           .drop_duplicates("gene_id").set_index("gene_id"))
de["neg_log10_padj"] = -np.log10(de["padj_gene_stageR"])

Y_CAP = 13.5
de["y_plot"]   = de["neg_log10_padj"].clip(upper=Y_CAP)
de["clipped"]  = de["neg_log10_padj"] > Y_CAP
de["category"] = "ns"
de.loc[(de.log2FoldChange >  THRESHOLDS["lfc"]) & (de.padj_gene_stageR < THRESHOLDS["padj"]), "category"] = "up"
de.loc[(de.log2FoldChange < -THRESHOLDS["lfc"]) & (de.padj_gene_stageR < THRESHOLDS["padj"]), "category"] = "down"

name_to_id = {row["gene_name"]: gid for gid, row in de.iterrows() if pd.notna(row.get("gene_name"))}

fig, ax = plt.subplots(figsize=(4, 3.5))

# Background
ns_mask = de["category"] == "ns"
ax.scatter(de.loc[ns_mask, "log2FoldChange"], de.loc[ns_mask, "y_plot"],
           c="#CCCCCC", s=3, alpha=0.3, linewidths=0, zorder=1, rasterized=True)
# Up
up_mask = (de["category"] == "up") & ~de["clipped"]
ax.scatter(de.loc[up_mask, "log2FoldChange"], de.loc[up_mask, "y_plot"],
           c="#D55E00", s=10, alpha=0.85, linewidths=0, zorder=2)
up_clip = (de["category"] == "up") & de["clipped"]
if up_clip.any():
    ax.scatter(de.loc[up_clip, "log2FoldChange"], [Y_CAP-0.2]*up_clip.sum(),
               c="#D55E00", s=14, alpha=0.85, marker="^", linewidths=0, zorder=3)
# Down
dn_mask = (de["category"] == "down") & ~de["clipped"]
ax.scatter(de.loc[dn_mask, "log2FoldChange"], de.loc[dn_mask, "y_plot"],
           c="#0072B2", s=10, alpha=0.85, linewidths=0, zorder=2)

# Threshold lines
ax.axvline( THRESHOLDS["lfc"], color="#BBBBBB", ls="--", lw=0.7)
ax.axvline(-THRESHOLDS["lfc"], color="#BBBBBB", ls="--", lw=0.7)
ax.axhline(-np.log10(THRESHOLDS["padj"]), color="#BBBBBB", ls="--", lw=0.7)

# Labels — selected genes
blue_genes, orange_texts = [], []
for gene in CIRCADIAN_GENES:
    gid = name_to_id.get(gene)
    if not gid or gid not in de.index or de.loc[gid,"category"] == "ns": continue
    x = de.loc[gid,"log2FoldChange"]; y = de.loc[gid,"y_plot"]
    if y < 1.35: continue
    col = "#D55E00" if de.loc[gid,"category"] == "up" else "#0072B2"
    ax.scatter(x, y, c=col, s=25, zorder=5, edgecolors="black", linewidths=0.5)
    if de.loc[gid,"category"] == "down":
        blue_genes.append((gene, x, y))
    else:
        orange_texts.append(ax.text(x+1.0, y, gene, fontsize=6, fontstyle="italic"))

adjust_text(orange_texts, ax=ax, expand=(2,2), force_text=(1.2,1.2), force_points=(1,1), lim=400)
orange_gene_list = [(g,de.loc[name_to_id[g],"log2FoldChange"],de.loc[name_to_id[g],"y_plot"])
                    for g in CIRCADIAN_GENES
                    if name_to_id.get(g) and name_to_id[g] in de.index
                    and de.loc[name_to_id[g],"category"]=="up"
                    and de.loc[name_to_id[g],"y_plot"]>=1.35]

for txt, (gene, dot_x, dot_y) in zip(orange_texts, orange_gene_list):
    lx, ly = txt.get_position(); txt.remove()
    ax.annotate(gene, xy=(dot_x,dot_y), xytext=(lx,ly), fontsize=6, fontstyle="italic",
                ha=txt.get_ha(), va="center",
                arrowprops=dict(arrowstyle="-", color="#888888", lw=0.5, shrinkA=3, shrinkB=3))

blue_sorted = sorted(blue_genes, key=lambda t: t[2], reverse=True)
label_pos = []
for gene, dot_x, dot_y in blue_sorted:
    lx = dot_x - 1.5; ly = dot_y
    for _, _, prev_ly in label_pos:
        if abs(ly - prev_ly) < 0.6: ly = prev_ly - 0.6
    label_pos.append((gene, lx, ly))

for (gene, dot_x, dot_y), (_, lx, ly) in zip(blue_sorted, label_pos):
    ax.annotate(gene, xy=(dot_x,dot_y), xytext=(lx,ly), fontsize=6, fontstyle="italic",
                ha="right", va="center",
                arrowprops=dict(arrowstyle="-", color="#888888", lw=0.5, shrinkA=0, shrinkB=0))

ax.set_xlabel("log₂FC (PerDKO / WT)")
ax.set_ylabel("−log₁₀(padj)")
ax.set_xlim(-13, 12); ax.set_ylim(0, Y_CAP+0.3)
ax.set_xticks([-10, -5, 0, 5, 10])
ax.set_yticks([0, 5, 10])

save_svg(fig, "Illumina_DE_genes_volcano.svg")

# =============================================================================
# FIGURE 4 — Illumina_Nanopore_DE_genes_concordance.svg  (step14, gene level)
# =============================================================================
print("\n=== Figure 4: DE gene concordance ===")

wb = openpyxl.load_workbook(KHUSHI, read_only=True, data_only=True)

ws_sig = wb["perKO_sig"]
rows_sig = list(ws_sig.iter_rows(values_only=True))
header_sig = [str(c) if c else "idx" for c in rows_sig[0]]
nano_sig_gid = {}
for r in rows_sig[1:]:
    row = dict(zip(header_sig, r))
    gid = row.get("gene_id")
    if gid and gid not in nano_sig_gid:
        try: nano_sig_gid[gid] = float(row["log2FC_perko"])
        except: pass

ws_all = wb["DE_all"]
rows_all = list(ws_all.iter_rows(values_only=True))
header_all = list(rows_all[0])
nano_all_gid = {}
for r in rows_all[1:]:
    row = dict(zip(header_all, r))
    gid = row.get("gene_id")
    if not gid: continue
    try: lfc = float(row["log2FC_perko"])
    except: continue
    if gid not in nano_all_gid or abs(lfc) > abs(nano_all_gid[gid]): nano_all_gid[gid] = lfc

de14 = pd.read_csv(DESEQ, sep="\t").dropna(subset=["log2FoldChange","padj_gene_stageR"])
de14["abs_lfc"] = de14["log2FoldChange"].abs()
de14_gene = de14.sort_values("abs_lfc",ascending=False).drop_duplicates("gene_id").set_index("gene_id")
ill_sig14 = set(de14_gene[(de14_gene.padj_gene_stageR < THRESHOLDS["padj"]) &
                           (de14_gene.abs_lfc > THRESHOLDS["lfc_15"])].index)
nano_sig14 = set(nano_sig_gid.keys())
common14   = set(de14_gene.index) & set(nano_all_gid.keys())
overlap14  = ill_sig14 & nano_sig14 & common14
print(f"  Common genes: {len(common14)}, overlap: {len(overlap14)}")

genes14 = list(common14)
x14 = np.array([float(de14_gene.loc[g,"log2FoldChange"]) for g in genes14])
y14 = np.array([float(nano_all_gid[g]) for g in genes14])

mask14_both = np.array([g in overlap14 for g in genes14])
mask14_ill  = np.array([g in ill_sig14 and g not in overlap14 for g in genes14])
mask14_nano = np.array([g in nano_sig14 and g not in overlap14 for g in genes14])
mask14_grey = ~(mask14_both | mask14_ill | mask14_nano)

MANUAL_POS = {
    "Orm1":     (-25,  9,  4, 4),
    "N4bp2l1":  (-18,  9,  4, 4),
    "Cyp4a12a": ( -7,  4,  0, 0),
    "Hsd17b6":  (  1, -1.63, 4, 4),
    "Hsp90aa1": (-6.36, 9, 4, 4),
}

fig, ax = plt.subplots(figsize=(4, 4))
ax.grid(False)
ax.scatter(x14[mask14_grey], y14[mask14_grey], c="#CCCCCC", s=3, alpha=0.35,
           linewidths=0, zorder=1, rasterized=True)
ax.scatter(x14[mask14_ill],  y14[mask14_ill],  c="#56B4E9", s=12, alpha=0.85,
           linewidths=0, zorder=2)
ax.scatter(x14[mask14_nano], y14[mask14_nano], c="#009E73", s=12, alpha=0.85,
           linewidths=0, zorder=3)
ax.scatter(x14[mask14_both], y14[mask14_both], c="#E69F00", s=30, alpha=0.95,
           edgecolors="black", linewidths=0.5, zorder=4)

lim14 = max(abs(x14).max(), abs(y14).max()) * 1.05
ax.plot([-lim14,lim14],[-lim14,lim14], color="#BBBBBB", lw=0.7, ls="--", alpha=0.6)
ax.axhline(0, color="#BBBBBB", lw=0.5); ax.axvline(0, color="#BBBBBB", lw=0.5)

texts14, gene_order14 = [], []
for i, g in enumerate(genes14):
    if not mask14_both[i]: continue
    gname = de14_gene.loc[g,"gene_name"] if "gene_name" in de14_gene.columns else gene_map.get(g,g)
    if gname in MANUAL_POS: continue
    ox = 1.0 if x14[i] >= 0 else -1.0
    oy = 0.6 if y14[i] >= 0 else -0.6
    texts14.append(ax.text(x14[i]+ox, y14[i]+oy, gname, fontsize=6, fontstyle="italic"))
    gene_order14.append((gname, x14[i], y14[i]))

adjust_text(texts14, ax=ax, expand=(3,3), force_text=(2,2), force_points=(1.5,1.5), lim=600)
for txt,(gname,dot_x,dot_y) in zip(texts14,gene_order14):
    lx,ly=txt.get_position(); txt.remove()
    ax.annotate(gname,xy=(dot_x,dot_y),xytext=(lx,ly),fontsize=6,fontstyle="italic",
                ha="center",va="center",
                arrowprops=dict(arrowstyle="-",color="#888888",lw=0.5,shrinkA=3,shrinkB=3))

for i,g in enumerate(genes14):
    if not mask14_both[i]: continue
    gname = de14_gene.loc[g,"gene_name"] if "gene_name" in de14_gene.columns else gene_map.get(g,g)
    if gname not in MANUAL_POS: continue
    lx,ly,sA,sB = MANUAL_POS[gname]
    ax.annotate(gname,xy=(x14[i],y14[i]),xytext=(lx,ly),fontsize=6,fontstyle="italic",
                ha="center",va="center",
                arrowprops=dict(arrowstyle="-",color="#888888",lw=0.5,shrinkA=sA,shrinkB=sB))

ax.set_xlabel("log₂FC Illumina (PerDKO / WT)")
ax.set_ylabel("log₂FC Nanopore (PerDKO / WT)")
ax.set_xticks([-20,-10,0,10,20])
ax.set_yticks([-20,-10,0,10,20])

save_svg(fig, "Illumina_Nanopore_DE_genes_concordance.svg")

# =============================================================================
# FIGURE 5 — Illumina_Nanopore_DE_transcript_concordance.svg  (step15, tx level)
# =============================================================================
print("\n=== Figure 5: DE transcript concordance ===")

ws_all2 = wb["DE_all"]
rows_all2 = list(ws_all2.iter_rows(values_only=True))
header_all2 = list(rows_all2[0])
khushi_all_tx = {}
for r in rows_all2[1:]:
    row = dict(zip(header_all2, r))
    txid = str(row.get("transcript_id","")).strip()
    if not txid or txid=="None": continue
    try: lfc=float(row["log2FC_perko"]); padj=float(row["padj_per_ko"])
    except: continue
    khushi_all_tx[txid] = {"log2FC":lfc,"padj":padj}

ws_sig2 = wb["perKO_sig"]
rows_sig2 = list(ws_sig2.iter_rows(values_only=True))
header_sig2 = [str(c) if c else "idx" for c in rows_sig2[0]]
khushi_sig_tx = set()
for r in rows_sig2[1:]:
    row = dict(zip(header_sig2,r))
    txid = str(row.get("transcript_id","")).strip()
    if txid and txid!="None": khushi_sig_tx.add(txid)

our_tx = pd.read_csv(DESEQ, sep="\t").dropna(subset=["log2FoldChange","padj_gene_stageR"])
our_tx = our_tx.set_index("transcript_id")
our_sig_tx = set(our_tx[(our_tx.padj_gene_stageR < THRESHOLDS["padj"]) &
                         (our_tx.log2FoldChange.abs() > THRESHOLDS["lfc_15"])].index)
tx_to_name = our_tx["gene_name"].to_dict()

common_tx = set(our_tx.index) & set(khushi_all_tx.keys())
overlap_tx = our_sig_tx & khushi_sig_tx & common_tx
print(f"  Common transcripts: {len(common_tx)}, overlap: {len(overlap_tx)}")

txids = list(common_tx)
x15 = np.array([float(our_tx.loc[t,"log2FoldChange"]) for t in txids])
y15 = np.array([float(khushi_all_tx[t]["log2FC"]) for t in txids])

mask15_both = np.array([t in overlap_tx for t in txids])
mask15_ill  = np.array([t in our_sig_tx and t not in overlap_tx for t in txids])
mask15_nano = np.array([t in khushi_sig_tx and t not in overlap_tx for t in txids])
mask15_grey = ~(mask15_both | mask15_ill | mask15_nano)

CLIP = 15
x15p = np.clip(x15, -CLIP, CLIP)
y15p = np.clip(y15, -CLIP, CLIP)

fig, ax = plt.subplots(figsize=(4, 4))
ax.grid(False)
ax.scatter(x15p[mask15_grey], y15p[mask15_grey], c="#CCCCCC", s=2, alpha=0.3,
           linewidths=0, zorder=1, rasterized=True)
ax.scatter(x15p[mask15_ill],  y15p[mask15_ill],  c="#56B4E9", s=10, alpha=0.85,
           linewidths=0, zorder=2)
ax.scatter(x15p[mask15_nano], y15p[mask15_nano], c="#009E73", s=10, alpha=0.85,
           linewidths=0, zorder=3)
ax.scatter(x15p[mask15_both], y15p[mask15_both], c="#E69F00", s=28, alpha=0.95,
           edgecolors="black", linewidths=0.5, zorder=4)

ax.plot([-CLIP,CLIP],[-CLIP,CLIP], color="#BBBBBB", lw=0.7, ls="--", alpha=0.6)
ax.axhline(0, color="#BBBBBB", lw=0.5); ax.axvline(0, color="#BBBBBB", lw=0.5)
ax.set_xlim(-CLIP*1.05, CLIP*1.05); ax.set_ylim(-CLIP*1.05, CLIP*1.05)

# Labels on orange
texts15, gene_order15 = [], []
for i,t in enumerate(txids):
    if not mask15_both[i]: continue
    gname = tx_to_name.get(t, t)
    ox = 0.8 if x15p[i] >= 0 else -0.8
    oy = 0.5 if y15p[i] >= 0 else -0.5
    texts15.append(ax.text(x15p[i]+ox, y15p[i]+oy, gname, fontsize=6, fontstyle="italic"))
    gene_order15.append((gname, x15p[i], y15p[i]))

adjust_text(texts15, ax=ax, expand=(3,3), force_text=(2,2), force_points=(1.5,1.5), lim=600)
for txt,(gname,dot_x,dot_y) in zip(texts15,gene_order15):
    lx,ly=txt.get_position(); txt.remove()
    ax.annotate(gname,xy=(dot_x,dot_y),xytext=(lx,ly),fontsize=6,fontstyle="italic",
                ha="center",va="center",
                arrowprops=dict(arrowstyle="-",color="#888888",lw=0.5,shrinkA=3,shrinkB=3))

ax.text(0.01,0.99,"▲ axes clipped ±15",transform=ax.transAxes,fontsize=5,color="#AAAAAA",va="top")
ax.set_xlabel("log₂FC Illumina (PerDKO / WT)")
ax.set_ylabel("log₂FC Nanopore (PerDKO / WT)")
ax.set_xticks([-10, -5, 0, 5, 10])
ax.set_yticks([-10, -5, 0, 5, 10])

save_svg(fig, "Illumina_Nanopore_DE_transcript_concordance.svg")

# =============================================================================
# LEGEND TXT
# =============================================================================
legend_text = """FIGURE LEGEND — Orthogonal Validation (GSE130613 vs Nanopore, PerDKO vs WT)
Generated: 2026-06-30
Dataset: Liver, CT16-20, constant darkness, n=4 per condition

=== Illumina_events_dPSI_volcano.svg ===
dPSI volcano plot — Illumina SUPPA2 diffSplice results (PerDKO vs WT).
X-axis: dPSI (PerDKO − WT). Y-axis: −log10(p-value).
Threshold lines: |dPSI| = 0.1 (vertical), p = 0.3 (horizontal).
Grey dots: not significant.
Coloured dots (by event type, significant in Illumina, |dPSI|>0.1, p<0.3):
  SE (skipping exon)       = #CC79A7 (pink/magenta)
  MX (mutually exclusive)  = #0072B2 (dark blue)
  RI (retained intron)     = #D55E00 (orange-red)
  A3 (alternative 3' ss)   = #E69F00 (amber)
  A5 (alternative 5' ss)   = #56B4E9 (light blue)
  AF (alternative first)   = #009E73 (teal/green)
  AL (alternative last)    = #F0E442 (yellow)
Total significant events: 514

=== Illumina_Nanopore_events_scatter.svg ===
dPSI concordance scatter — Illumina vs Nanopore splicing (PerDKO vs WT).
X-axis: dPSI Illumina (PerDKO − WT). Y-axis: dPSI Nanopore (PerDKO − WT).
Dashed diagonal: y=x line (perfect agreement).
Thresholds: Illumina |dPSI|>0.1, p<0.3 | Nanopore |dPSI|>0.1, p<0.05.
Grey dots:   not significant in either platform (n=22)
Blue dots:   significant in Illumina only (n=481)
Green dots:  significant in Nanopore only (n=888)
Orange dots: significant in BOTH platforms — orthogonal validation (n=32)
Total events in both datasets: 1,423. Labels on orange dots only (gene names).

=== Illumina_DE_genes_volcano.svg ===
DE volcano plot — Illumina DESeq2 + stageR results (PerDKO vs WT), gene level.
X-axis: log2 fold change (PerDKO / WT). Y-axis: −log10(padj, stageR stage 1).
Threshold lines: |log2FC| = 0.5 (vertical), padj = 0.05 (horizontal).
Grey dots:   not significant
Orange/red dots (#D55E00): upregulated in PerDKO (padj<0.05, log2FC>0.5)
Blue dots (#0072B2):       downregulated in PerDKO (padj<0.05, log2FC<-0.5)
Triangles (▲): points clipped at y=13.5 (extreme outliers).
Labels: circadian and selected genes (Kiran's list).

=== Illumina_Nanopore_DE_genes_concordance.svg ===
DE concordance scatter — Illumina vs Nanopore, GENE level (PerDKO vs WT).
X-axis: log2FC Illumina. Y-axis: log2FC Nanopore.
One dot per gene (representative transcript = highest |log2FC|).
Threshold (both sides): padj<0.05 AND |log2FC|>1.5.
Grey dots:   present in both datasets, not significant in either (n=9,519)
Blue dots:   significant in Illumina only (n=188)
Green dots:  significant in Nanopore only (n=116)
Orange dots: significant in BOTH platforms (n=27) — labeled with gene names
Total genes in both datasets: 9,850.
See: step14_step15_gene_transcript_mapping.xlsx

=== Illumina_Nanopore_DE_transcript_concordance.svg ===
DE concordance scatter — Illumina vs Nanopore, TRANSCRIPT level (PerDKO vs WT).
X-axis: log2FC Illumina. Y-axis: log2FC Nanopore. Axes clipped at ±15.
One dot per transcript (includes annotated ENSMUST and novel UUID isoforms).
Both datasets used the same custom transcriptome (transcriptome_ext.fa).
Threshold (both sides): padj<0.05 AND |log2FC|>1.5.
Grey dots:   present in both datasets, not significant in either (n=12,602)
Blue dots:   significant in Illumina only (n=256)
Green dots:  significant in Nanopore only (n=92)
Orange dots: significant in BOTH platforms (n=16) — labeled with gene names
              (11 annotated ENSMUST + 5 novel UUID transcripts)
              11/16 show consistent direction of change across platforms
Total transcripts in both datasets: 12,966.
See: step15_DE_overlap_Illumina_vs_Nanopore.xlsx
     step15_khushi_267_transcript_breakdown.xlsx

Note: Khushi's 267 Nanopore-significant transcripts are accounted for as:
  16 orange + 92 green + 18 NA padj + 141 zero reads in Illumina = 267
"""

with open(os.path.join(OUT, "figure_legend.txt"), "w") as f:
    f.write(legend_text)
print("\nSaved: figure_legend.txt")
print("\n=== ALL DONE ===")
