# Internship — Alternative Splicing Analysis of the Murine Circadian Transcriptome

**BSc in Bioinformatics internship** at the Institut de Genomique Fonctionnelle de Lyon (IGFL), ENS Lyon.  
Supervised by Dr. Kiran Padmanabhan. April 27 – August 28, 2026.

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
├── IGV/                            # IGV application and saved sessions
├── SUPPA/                          # SUPPA2 tool
└── README.md
```

## Tools and Languages

- **SUPPA2** v2.3 — alternative splicing event detection and PSI quantification
- **Salmon** v1.5.2 — transcript quantification (TPM)
- **IGV** v2.19.7 — genome browser for visual inspection
- **Python** 3.8.20 (`suppa_env`) — data analysis and figure generation
- **Bash** — pipeline automation

## Notes

- All large data files (GTF, BAM, FASTQ, SUPPA outputs) are excluded via `.gitignore`.
- Run `report_figures.py` from the `Internship/` root to regenerate all report figures.
- Bash scripts use dynamic paths relative to their location — no hardcoded absolute paths.
