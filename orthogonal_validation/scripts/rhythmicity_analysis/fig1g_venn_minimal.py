#!/usr/bin/env python3
"""
Fig 1G minimal -- same overlap as fig1g_venn.py (773 rhythmic genes and the
228 rhythmic AS genes), circles only, no counts/labels, matching the
manuscript's thicker black outline. Kiran adds the numbers/text via textbox
in pptx.
Author: Gricey

Run from: orthogonal_validation/
    python3 scripts/rhythmicity_analysis/fig1g_venn_minimal.py
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

fig, ax = plt.subplots(figsize=(6.5, 3.8))

v = venn2(subsets=(727, 182, 46), set_labels=("", ""), ax=ax,
          set_colors=("#b3b3e6", "#e8b3c2"), alpha=0.85)

for text in v.subset_labels:
    if text:
        text.set_text("")

# matplotlib_venn's auto-blend for the overlap patch washes out to near-white
# with these pastel colors -- set it explicitly to the purple/pink blend seen
# in the manuscript figure instead.
overlap = v.get_patch_by_id("11")
if overlap:
    overlap.set_facecolor("#ceb3d4")
    overlap.set_alpha(0.85)

for patch_id in ("10", "01", "11"):
    patch = v.get_patch_by_id(patch_id)
    if patch:
        patch.set_edgecolor("black")
        patch.set_linewidth(2.5)

plt.tight_layout()
plt.savefig("results/rhythmicity_analysis/gene_rhythmicity/figures/fig1g_venn_773_minimal.svg")
