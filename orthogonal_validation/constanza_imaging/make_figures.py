"""
Build colocalization figures from tidy_colocalization.csv.
Metric shown: Pearson's R. Strip/dot plot, one dot per cell/nucleus, mean +/- SD line.

Fig 1: three separate figures, one per pair (PER2/SRSF3, PER2/SC35, SC35/SRSF3),
       each showing HepG2 endogenous vs HA-SRSF3 overexpression.
Fig 2: mouse liver nuclei, CT20, Z3 (Zoom 3) only -- single-column strip plot.
Fig 3: not needed.
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_theme(style="ticks", context="paper", font_scale=1.1)
df = pd.read_csv("tidy_colocalization.csv")
pear = df[df.metric == "Pearson's R"].copy()


def save_fig(fig, fname_png):
    fig.savefig(fname_png, dpi=300, bbox_inches="tight")
    fig.savefig(fname_png.replace(".png", ".pdf"), bbox_inches="tight")
    fig.savefig(fname_png.replace(".png", ".svg"), bbox_inches="tight")


def strip_with_mean(ax, data, x, order, color=None, palette=None, jitter=0.15):
    sns.stripplot(data=data, x=x, y="value", order=order, hue=x, ax=ax, size=6,
                  jitter=jitter, alpha=0.75, linewidth=0.4, edgecolor="white",
                  color=color, palette=palette, legend=False)
    for i, cat in enumerate(order):
        vals = data.loc[data[x] == cat, "value"]
        if len(vals) == 0:
            continue
        m, sd = vals.mean(), vals.std()
        ax.plot([i - 0.25, i + 0.25], [m, m], color="black", lw=1.8, zorder=5)
        ax.errorbar(i, m, yerr=sd, color="black", capsize=3, lw=1, zorder=4)
    ax.set_ylim(-0.05, 1.05)
    ax.set_ylabel("Pearson's R")
    ax.set_xlabel("")


# ============================================================
# FIGURE 1 — one figure per pair, HepG2 endogenous vs overexpression
# ============================================================
pairs = [
    ("PER2 vs SRSF3", "fig1a_PER2vsSRSF3.png", "#2b6cb0"),
    ("PER2 vs SC35", "fig1b_PER2vsSC35.png", "#a0aec0"),
    ("SC35 vs SRSF3", "fig1c_SC35vsSRSF3.png", "#718096"),
]
cond_order = ["endogenous", "HA-SRSF3 overexpression"]

for comparison, fname, color in pairs:
    fig, ax = plt.subplots(figsize=(3.5, 4.5))
    sub = pear[(pear.cell_line == "HepG2") & (pear.comparison == comparison)]
    strip_with_mean(ax, sub, "condition", cond_order, color=color)
    ns = [sub[sub.condition == c].shape[0] for c in cond_order]
    ax.set_xticks(range(len(cond_order)))
    ax.set_xticklabels([f"Endogenous\n(n={ns[0]})", f"HA-SRSF3\noverexpression\n(n={ns[1]})"], fontsize=8)
    ax.set_title(f"Fig 1 — {comparison}\n(HepG2)", fontsize=11)
    fig.tight_layout()
    save_fig(fig, fname)
    plt.close(fig)

# ============================================================
# FIGURE 2 — mouse liver nuclei, CT20, Z3 (Zoom 3) only
# ============================================================
fig, ax = plt.subplots(figsize=(3, 4.5))
d_mouse = pear[(pear.cell_line == "mouse liver (in vivo)") & (pear.condition == "CT20 Z=3")].copy()
d_mouse["label"] = "PER2 / SRSF3"
strip_with_mean(ax, d_mouse, "label", ["PER2 / SRSF3"], color="#6b46c1")
n = d_mouse.shape[0]
ax.set_xticks([0])
ax.set_xticklabels([f"PER2 / SRSF3\n(n={n})"])
ax.set_title("Fig 2 — PER2/SRSF3 colocalization\nmouse liver nuclei, CT20 (Zoom 3)", fontsize=11)
fig.tight_layout()
save_fig(fig, "fig2_mouse_liver_nuclei_Z3.png")
plt.close(fig)

# ============================================================
# FIGURE 1 (endogenous-only versions) — same 3 pairs, single column, no overexpression
# ============================================================
pairs_endog = [
    ("PER2 vs SRSF3", "fig1a_PER2vsSRSF3_endogenous.png", "#2b6cb0"),
    ("PER2 vs SC35", "fig1b_PER2vsSC35_endogenous.png", "#a0aec0"),
    ("SC35 vs SRSF3", "fig1c_SC35vsSRSF3_endogenous.png", "#718096"),
]

for comparison, fname, color in pairs_endog:
    fig, ax = plt.subplots(figsize=(2.4, 4.5))
    sub = pear[(pear.cell_line == "HepG2") & (pear.comparison == comparison) & (pear.condition == "endogenous")].copy()
    sub["label"] = comparison
    strip_with_mean(ax, sub, "label", [comparison], color=color)
    n = sub.shape[0]
    ax.set_xticks([0])
    ax.set_xticklabels([f"{comparison}\n(n={n})"], fontsize=9)
    ax.set_title(f"Fig 1 — {comparison}\n(HepG2, endogenous)", fontsize=11)
    fig.tight_layout()
    save_fig(fig, fname)
    plt.close(fig)

# ============================================================
# FIGURE 1 (combined) — all 3 pairs, endogenous only, side by side in one figure
# ============================================================
fig, ax = plt.subplots(figsize=(6, 4.5))
comparison_order = ["PER2 vs SRSF3", "PER2 vs SC35", "SC35 vs SRSF3"]
palette = {"PER2 vs SRSF3": "#2b6cb0", "PER2 vs SC35": "#a0aec0", "SC35 vs SRSF3": "#718096"}
d_endo_all = pear[(pear.cell_line == "HepG2") & (pear.condition == "endogenous")]
strip_with_mean(ax, d_endo_all, "comparison", comparison_order, palette=palette)
ns = [d_endo_all[d_endo_all.comparison == c].shape[0] for c in comparison_order]
ax.set_xticks(range(len(comparison_order)))
ax.set_xticklabels([f"{c}\n(n={n})" for c, n in zip(comparison_order, ns)], fontsize=9)
ax.set_title("Fig 1 — Colocalization specificity in HepG2 cells\n(endogenous SRSF3)", fontsize=11)
fig.tight_layout()
save_fig(fig, "fig1_combined_endogenous.png")
plt.close(fig)

print("Saved fig1a/b/c_*.{png,pdf,svg} (endog vs overexpression)")
print("Saved fig1a/b/c_*_endogenous.{png,pdf,svg} (endogenous only, single column)")
print("Saved fig1_combined_endogenous.{png,pdf,svg} (all 3 pairs, endogenous only)")
print("Saved fig2_mouse_liver_nuclei_Z3.{png,pdf,svg}")
