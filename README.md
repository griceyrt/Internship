# Internship — Alternative Splicing Analysis of the Murine Circadian Transcriptome

**BSc in Bioinformatics internship** at the Institut de Genomique Fonctionnelle de Lyon (IGFL), ENS Lyon.  
Supervised by Dr. Kiran Padmanabhan. April 27 – August 28, 2026.

---

## Projects

### Project 1 · SUPPA2 Boundary Parameter Exploration

This project explores how SUPPA2's boundary parameter (strict vs. variable mode, 1–20 nt) affects alternative splicing event detection and quantification in Nanopore long-read direct RNA-seq data from mouse liver across three circadian conditions (CT8, CT20, PerKO). Rather than benchmarking event-calling accuracy, the analysis characterises how the choice of boundary mode changes detected event counts and PSI values across the seven SUPPA2 event types (SE, RI, A3, A5, MX, AF, AL), confirming variable mode as a superset of strict mode at both the IOE and PSI levels and identifying candidate splicing events with notable PSI shifts between conditions (e.g. the C8g A3 event).

### Project 2 · Orthogonal Validation (in progress)

Orthogonal validation means confirming a result using a completely different technology. The lab previously found differential splicing events in PerDKO vs. WT mouse liver using Nanopore long-read sequencing. This project attempts to confirm the same events using public Illumina short-read RNA-seq data (GSE171975, Asher lab, Weizmann Institute) — if the same splicing changes appear independently in both datasets, the finding is much more robust.

The pipeline runs: SRA download → Salmon transcript quantification (aligned to the lab's custom transcriptome) → tximport into DESeq2 (normalised expression) → SUPPA2 alternative splicing analysis → comparison of dPSI values against the Nanopore results.

---

## Repository Structure

```
Internship/
│
├── boundary_analysis/              # Project 1: SUPPA2 boundary parameter exploration
│   ├── data/                       # Raw and processed data (gitignored)
│   ├── scripts/                    # Bash pipeline scripts (SUPPA2, IGV filtering)
│   ├── notebooks/                  # Jupyter notebooks for analysis and figures
│   ├── figures/
│   │   ├── plots/                  # Output figures (barplots, Venns, density plots)
│   │   └── igv/                    # IGV screenshots per gene and event type
│   └── presentations/              # Keynote slides for lab meetings
│
├── orthogonal_validation/          # Project 2: Illumina short-read validation (in progress)
│   ├── data/                       # Raw and processed data (gitignored)
│   ├── scripts/                    # Pipeline scripts (SRA download, Salmon, DESeq2, SUPPA2)
│   ├── notebooks/                  # Jupyter notebooks (in progress)
│   └── presentations/              # Keynote slides for lab meetings
│
├── report/                         # UCBL internship report
│   ├── figures/                    # Report figures (fig1-fig5) + report_figures.py
│   └── docs/                       # Report PDF (when finalized)
│
├── papers/                         # Reference PDFs (gitignored)
├── IGV/                            # IGV application and saved sessions
├── SUPPA/                          # SUPPA2 tool
└── README.md
```

A visual overview of how the Project 1 scripts, notebooks, and data folders connect is available in [`ucbl_report/figures/web/suppa_project_overview.svg`](ucbl_report/figures/web/suppa_project_overview.svg).

---

## Tools and Languages

- **SUPPA2** v2.3 — alternative splicing event detection and PSI quantification
- **Salmon** v1.5.2 — transcript quantification (TPM)
- **IGV** v2.19.7 — genome browser for visual inspection
- **Python** 3.8.20 (`suppa_env`) — data analysis and figure generation
- **R** — tximport / DESeq2 for normalised expression (Project 2)
- **Bash** — pipeline automation

---

## Setup

1. Clone the repository:
   ```bash
   git clone https://github.com/griceyrt/Internship.git
   cd Internship
   ```

2. Create and activate the Python environment:
   ```bash
   python3 -m venv suppa_env
   source suppa_env/bin/activate
   pip install pandas matplotlib matplotlib-venn jupyter
   ```

3. Install SUPPA2 by following the instructions at:
   [https://github.com/comprna/SUPPA](https://github.com/comprna/SUPPA)
