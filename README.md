# Internship — Alternative Splicing Analysis of the Murine Circadian Transcriptome

**BSc in Bioinformatics internship** at the Institut de Génomique Fonctionnelle de Lyon (IGFL), ENS Lyon.
Supervised by Dr. Kiran Padmanabhan, with Bharath, Khushi Mamgain, and Lies. April 27 – August 28, 2026.

**Handoff note (2026-08-21, last internship day):** four projects were carried out over the internship, at various stages of completion. Each project folder has its own detailed `README.md` — this file is just the index. Start there for anything specific; this file is for orienting a newcomer to the whole internship's scope.

---

## Projects

### Project 1 · SUPPA2 Boundary Parameter Exploration (`boundary_analysis/`) — DONE

Explores how SUPPA2's boundary parameter (strict vs. variable mode, 1–20 nt) affects alternative splicing event detection and PSI quantification in Nanopore long-read direct RNA-seq data from mouse liver across three circadian conditions (CT8, CT20, PerKO). Confirmed variable mode as a superset of strict mode at both the IOE and PSI level, and identified candidate splicing events with notable PSI shifts between conditions. No handoff README written for this one (predates the handoff convention) — see `notebooks/` and `ucbl_report/` for the write-up.

### Project 2 · Orthogonal Validation (`orthogonal_validation/`) — see its own README

Validates splicing/DE events from the lab's manuscript (PER2/circadian clock, in revision) against independent public Illumina short-read datasets, since the lab's own primary data is Nanopore. Several parallel workstreams: GSE130613 (closed), GSE171975 (active, awaiting reply), GSE133398/SRSF3 (paused on an unresolved event-count gap), E-MTAB-9701 (dropped, insufficient replicates), plus two closed-out manuscript rebuttal figure reworks and a wet-lab imaging figure subproject. **Read `orthogonal_validation/README.md` for the full status table before picking anything up.**

### Project 3 · Fmr1 KO Poly(A)/Splicing (`fmr1_polya_splicing/`) — see its own README

GSE188840 — WT vs Fmr1 KO mouse brain cortex, ONT direct RNA-seq. Covers isoform usage/entropy analysis, GO enrichment (per-GO-term paired Wilcoxon method, not standard hypergeometric), an Overlap #1 (DESeq2 vs SUPPA significant genes) figure set, and a poly(A) tail length analysis (Nanopolish) that was **still running when the internship ended** — the last job submitted was not confirmed complete. **Read `fmr1_polya_splicing/README.md`** — it documents exactly what's built vs. what's mid-run vs. what was never started (custom transcriptome/IsoQuant, Overlap #2).

### Project 4 · HFD Organoid Splicing Landscape (`hfd_organoid_splicing/`) — see its own README

Adds a splicing-landscape layer on top of Desmond's already-completed transcriptomics/proteomics on a High Fat Diet mouse liver organoid time course (5 timepoints × 3 replicates, Illumina). **Core deliverables complete and already emailed to Kiran/Desmond.** One real open item: Desmond flagged a possible threshold mismatch in the DE-gene-list comparison used for one figure — **read `hfd_organoid_splicing/README.md`'s "what's genuinely unfinished" section** before touching this again.

---

## Repository structure

```
Internship/
│
├── boundary_analysis/          Project 1 — SUPPA2 boundary parameter exploration
├── orthogonal_validation/      Project 2 — Illumina short-read validation + rebuttal figures (own README)
├── fmr1_polya_splicing/        Project 3 — Fmr1 KO ONT poly(A)/splicing (own README)
├── hfd_organoid_splicing/      Project 4 — HFD organoid splicing landscape (own README)
├── ucbl_report/                UCBL internship report (figures + docs)
├── papers/                     Reference PDFs (gitignored)
├── IGV/                        IGV application and saved sessions
├── SUPPA/                      SUPPA2 tool (own README)
└── README.md                   this file
```

Each project folder generally follows the same internal convention: `scripts/` (numbered by pipeline order where practical), `results/` or `data/`, `figures/`, `notes/`, `meta/`. Where a project's own README documents a different convention (e.g. Project 3's `scripts/` numbering is build-order, not the numbering used in its workflow diagram), trust that project's README over this general pattern.

## Cluster access (applies to all 4 projects)

PSMN cluster, ENS Lyon, username `grangel`. Full details (SSH jump-host chain, storage tiers, `biotools` conda env, known transfer gotchas) are documented once per project's own README/notes where relevant — the fullest version is in `orthogonal_validation/notes/` (project 2 set up the cluster workflow first) and repeated/extended in each later project's README as new gotchas were found.

## Tools and languages

- **Salmon**, **STAR**, **rMATS-turbo**, **SUPPA2**, **DESeq2**/**tximport**/**stageR**, **Nanopolish**, **minimap2**, **oarfish** — quantification, alignment, and splicing/DE tools, used selectively per project depending on data type (short-read Illumina vs. ONT long-read) and library prep (standard vs. MARS-seq 3')
- **clusterProfiler**, **GO.db**/**org.Mm.eg.db** — GO enrichment, both standard hypergeometric and (project 3) a custom paired-Wilcoxon method
- **IGV** — genome browser, visual validation of novel transcripts/splicing events
- **Python** (pandas, matplotlib, openpyxl) and **R** (Bioconductor) — figure and table generation
- **Bash/SLURM** — pipeline automation on the PSMN cluster

## Setup

1. Clone the repository:
   ```bash
   git clone https://github.com/griceyrt/Internship.git
   cd Internship
   ```
2. Each project's own README lists its specific conda/venv environment(s) — there is no single shared environment across all 4 projects (they diverged as tool needs diverged, e.g. project 3/4's GO enrichment work uses a cloned `go_entropy`-family env kept separate from `biotools` to avoid cross-project breakage).
3. SUPPA2 install instructions: [https://github.com/comprna/SUPPA](https://github.com/comprna/SUPPA)
