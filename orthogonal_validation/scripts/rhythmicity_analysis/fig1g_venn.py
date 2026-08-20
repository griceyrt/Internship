#!/usr/bin/env python3
"""
Fig 1G redo -- overlap between the 773 rhythmic genes (gene_to_tx_rhythmicity.R)
and the 228 rhythmic AS genes (data/tables/Table5-rhythmic_AS_events.xlsx).
Author: Gricey

Run from: orthogonal_validation/
    python3 scripts/fig1g_venn.py
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

fig, ax = plt.subplots(figsize=(6.5, 3.8))

v = venn2(subsets=(727, 182, 46), set_labels=("Rhythmic genes", "Rhythmic AS events"), ax=ax,
          set_colors=("#b3b3e6", "#e8b3c2"), alpha=0.7)

for text in v.set_labels:
    if text:
        text.set_fontsize(13)
        text.set_fontweight("bold")
for text in v.subset_labels:
    if text:
        text.set_fontsize(13)

for patch_id in ("10", "01", "11"):
    patch = v.get_patch_by_id(patch_id)
    if patch:
        patch.set_edgecolor("black")
        patch.set_linewidth(1.5)

plt.tight_layout()
plt.savefig("results/rhythmicity_analysis/gene_rhythmicity/figures/fig1g_venn_773.svg")
