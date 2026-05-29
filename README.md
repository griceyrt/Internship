# Circadian Rhythm & Alternative Splicing

**Internship project at the Institut de Génomique Fonctionnelle de Lyon (ENS Lyon)**  
Supervisor: Khushi Mamgain · Kiran Padmanabhan  

---

## Overview

This project investigates **alternative splicing** in the context of circadian rhythms, using mouse RNA-seq data (GRCm39/mm39) across three biological conditions:

| Condition | Samples | Description |
|---|---|---|
| CT8 | KM1, KM2, KM3 | Circadian time 8 |
| CT20 | KM4, KM5, KM6 | Circadian time 20 |
| PerKO | KM7, KM8, KM9 | Period gene knockout |

The central question: **how do PERIOD proteins influence alternative splicing, and how does the choice of splice site boundary mode affect what we detect?**

This work builds on: *"PERIOD proteins dictate targeting of splicing factors to chromatin, shaping circadian transcriptome diversity"* — Chikhaoui & Mamgain et al.

---

## Repository Structure

```
Internship/
├── analysis/
│   ├── notebooks/
│   │   ├── 00_numbers.ipynb          # Event counts across boundary modes
│   │   ├── 01_filter_events.ipynb    # Unique events: variable vs strict
│   │   ├── 02_barplots.ipynb         # IOE population bar plots
│   │   └── 03_psi_comparison.ipynb   # PSI comparison: strict vs variable 1nt
│   └── 01_strict_vs_1nt.key          # Keynote presentation slides
├── scripts/
│   ├── run_suppa.sh                  # Full SUPPA pipeline (generateEvents → clusterEvents)
│   ├── compare_runs.sh               # Compare significant events across boundary modes
│   ├── compare_event_populations.sh  # Q1/Q2/Q3 gene population analysis
│   ├── get_igv_candidates.sh         # Filter candidates for IGV visualisation
│   └── igv_filter.sh                 # IGV region filtering by event type
├── figures/
│   └── plots/                        # Generated tables and charts
└── data/                             # (gitignored) SUPPA outputs, Salmon TPM, GTF files
```

---

## Methods

**Tool:** [SUPPA2](https://github.com/comprna/SUPPA)  
**Reference genome:** Mouse GRCm39/mm39  
**Annotation:** `transcriptome_productivity.gtf`  
**Quantification:** Salmon (TPM)

### Pipeline

1. `generateEvents` — detect splicing events from GTF in strict or variable boundary mode
2. `psiPerEvent` — quantify splicing usage (PSI) per event using Salmon TPM files
3. `diffSplice` — identify differentially spliced events between conditions
4. `clusterEvents` — cluster events by PSI profile (SE events)

### Boundary modes compared

| Mode | Description |
|---|---|
| Strict | Exact splice site match required |
| Variable 1nt | ±1 nucleotide tolerance at splice boundaries |
| Variable 3–20nt | Broader tolerances for sensitivity analysis |

---

## Key Findings

### IOE-level: what does relaxing the boundary unlock?

| Event | Q1 gained | Q2 new (0→events) | Q3 novel transcripts |
|---|---|---|---|
| A3 | 2,711 genes | 3,423 genes | 2,881 events |
| A5 | 2,163 genes | 3,097 genes | 2,530 events |
| RI | 1,618 genes | 4,405 genes | 2,286 events |

- **Q1** — genes already in strict that gain more events in variable
- **Q2** — genes completely absent in strict, newly detected in variable
- **Q3** — events involving unannotated (UUID-format) transcript IDs

The majority of the increase comes from **Q2**: genes that strict missed entirely, not just genes gaining a few extra events.

### PSI-level: are these new events biologically expressed?

Of the variable-only events with ENSMUSG gene IDs:
- ~65% have mean PSI > 0 across samples — genuinely expressed
- The remainder have no TPM support (annotation artefacts)

---

## Usage

All scripts are run from the `Internship/` root. Paths to `data/` are hardcoded as absolute paths inside each script.

```bash
# Run full SUPPA pipeline (edit BOUNDARY and THRESHOLD at top of script first)
bash scripts/run_suppa.sh

# Compare event populations across boundary modes
bash scripts/compare_event_populations.sh

# Get IGV candidates for a given event type
bash scripts/get_igv_candidates.sh

# Filter IGV candidates interactively
bash scripts/igv_filter.sh RI 3 1 4
```

---

## Environment

```bash
conda activate suppa_env   # Python 3.8.20
```

Dependencies: SUPPA2, Python (pandas, numpy, matplotlib, seaborn), Salmon, IGV
