# Colocalization Analysis — PER2 / SRSF3 Confocal Data
**Context:** Constanza (wet lab, IGFL) is analyzing confocal images from U2OS-PER2-GFP-HA-SRSF3-WT cells using Fiji (Coloc2 plugin). The goal is to demonstrate colocalization between PER2 and SRSF3 across 25 cells, one condition.

This is directly relevant to the Chikhaoui & Mamgain et al. manuscript, which shows via IP (Fig. 4d) that endogenous PER2 physically interacts with SRSF3 at CT20. The confocal colocalization is the spatial/visual confirmation of that interaction.

---

## Metrics output by Fiji / Coloc2

| Metric | Range | Interpretation |
|---|---|---|
| **Pearson's R** | -1 to +1 | Overall pixel intensity correlation; >0.5 generally accepted as colocalization |
| **Spearman's ρ** | -1 to +1 | Rank-based correlation; more robust to outliers |
| **Li's ICQ** | -0.5 to +0.5 | >0 = colocalized; 0 = random; <0 = exclusion |
| **Costes P-value** | 0–1 | Statistical significance of colocalization above chance; want >0.95 |
| **tM1** | 0–1 | Fraction of PER2 signal overlapping with SRSF3 (Manders' thresholded) |
| **tM2** | 0–1 | Fraction of SRSF3 signal overlapping with PER2 (Manders' thresholded) |

**tM1 and tM2** (Manders' coefficients) are the most biologically interpretable and most requested by reviewers.

---

## What Fiji is already generating

### Cytofluorogram (scatter plot)
- X-axis: PER2 channel intensity per pixel
- Y-axis: SRSF3 channel intensity per pixel
- Each dot = one pixel
- The tighter and more diagonal the cloud → the higher the colocalization
- Pearson's R is the correlation coefficient of this scatter

### Li's ICQ plot (what's on Constanza's screen)
- X-axis: Channel 1 intensity minus its mean (centered at 0)
- Y-axis: Channel 2 intensity minus its mean (centered at 0)
- The **bowtie/hourglass shape** = good colocalization signature
- Color = pixel density (orange/red = high density core, purple = outliers)
- ICQ value = proportion of pixels where both channels are simultaneously above OR below their mean

Both of these are **exploratory/representative** plots — good for supplementary or methods, not the main quantification figure.

---

## Recommended plot for publication

### Single-column strip/dot plot (one condition, 25 cells)

**Setup:**
- Y-axis: Pearson's R (0 to 1)
- X-axis: one label only — "PER2 / SRSF3"
- Each dot = one cell (25 dots total)
- Horizontal line at the **mean** (or median)
- Optional: thin box or violin behind the dots to show distribution

**Why this works:**
- No groups to compare → no need for bars or multiple columns
- 25 cells is enough to show the distribution is consistently high
- Pearson's R is the most universally understood metric
- Costes P-value can be reported in the figure legend as statistical validation

**Conceptual layout:**
```
1.0 |        •  •
    |     •  •  •  •
0.8 |   •  •  •  •  •   ——— mean
    |      •  •  •
0.6 |   •     •
    |
    ——————————————————
       PER2 / SRSF3
```

---

## Suggested figure structure (main figure panel)

| Panel | Content |
|---|---|
| A | Merged confocal images (PER2 channel, SRSF3 channel, merge) — representative cell |
| B | Cytofluorogram or Li ICQ plot — representative cell |
| C | Strip/dot plot of Pearson's R across 25 cells with mean line |

Report in figure legend: mean ± SD, n=25 cells, Costes P-value.

---

## Important technical note

Pearson's R is sensitive to background pixels — a large dark background artificially inflates the correlation. Coloc2's **Costes automatic thresholding** handles this. Confirm it is enabled before reporting R values.

---

## Tools to make the final plot

- **Prism** — easiest for Constanza if she already uses it; "scatter with mean" plot type
- **R (ggplot2)** — `geom_jitter()` + `stat_summary()` for mean line
- **Python (matplotlib/seaborn)** — `stripplot()` + `axhline()` for mean

Data input = one column of 25 Pearson's R values exported from Fiji/Coloc2.
