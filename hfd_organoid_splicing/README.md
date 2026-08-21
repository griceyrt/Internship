# HFD Organoid Splicing Landscape (Project 4) — handoff README

**Author:** Gricey

**Last updated:** 2026-08-21, my last internship day. **Status: core deliverables complete and already emailed to Kiran/Desmond** (2026-08-11/12). One real open follow-up from Desmond's reply is described below — pick that up first if continuing this project.

## Goal

Desmond (lab member) already completed transcriptomics (DE) and proteomics on a High Fat Diet (HFD) mouse liver organoid time course. My job was to add the missing layer: **does alternative splicing change over the HFD time course, and if so, what biological pathways does it hit (GO enrichment)?** This is a genuinely distinct question from Desmond's DE analysis — a gene can be spliced differently with no expression change at all (confirmed empirically: only 17-23% overlap between splicing-significant and DE-significant gene sets at every timepoint).

Biological framing from Kiran: HFD disrupts the liver's normal circadian rhythm (same broader theme as Project 2's PerDKO work) progressively over time, eventually producing damaged tissue that can progress toward carcinoma. Early damage is reversible; beyond some threshold it isn't. GO enrichment across timepoints was meant to trace that progression.

## Dataset

15 single-end Illumina fastq files, 5 timepoints × 3 replicates, mouse liver organoids:

| Timepoint | Label | Files |
|---|---|---|
| Baseline | WT_SD | Lib1-3 |
| 8 weeks HFD | W8_WHFD | Lib4-6 |
| 18 weeks HFD | W18_WHFD | Lib7-9 |
| 24 weeks HFD | W24_WHFD | Lib10-12 |
| 42 weeks HFD | W42_WHFD | Lib13-15 |

Reference: Ensembl mouse **GRCm39, release 116** (not 115 — Desmond's own script had a stale reference in a comment; he confirmed 116 by email after checking what he actually downloaded). ~490M raw reads total, 97.78% retained after trimming, 77.17% overall Salmon mapping rate, no outlier samples.

## Pipeline

All comparisons are timepoint-vs-baseline (W8/W18/W24/W42, each vs WT_SD).

1. `00_download_reference.sh` — genome + cDNA + GTF, Ensembl release 116
2. `01_raw_qc.sh` — FastQC + MultiQC (raw)
3. `02_fastp_trim.sh` — poly-G trim, quality sliding-window, TruSeq adapter auto-detect, <35bp discard
4. `03_salmon_index.sh` — decoy-aware Salmon index (genome as decoy + cDNA)
5. `04_salmon_quant.sh` — `salmon quant -l SF --seqBias --gcBias --validateMappings` (single-end forward stranded, confirmed via RSeQC)
6. `05_suppa_generate_events.sh` — all 7 SUPPA2 event types (SE, A3, A5, MX, RI, AF, AL)
7. `06_build_tpm_matrix.sh` — combines 15 quant.sf into one TPM table (own script, Desmond's pipeline never needed this)
8. `07_suppa_psi.sh` — PSI per event (777,087 rows)
9. `08_suppa_diffsplice.sh` — 28 comparisons (4 timepoints × 7 event types) vs WT_SD, empirical method
10. `09_extract_significant_events.sh` — |dPSI|>0.1, p<0.05 (same thresholds as the lab's own PERIOD paper)
11. `10_extract_gene_symbols.sh` — GTF-derived ID→symbol mapping
12. `11_go_enrichment.R` / `.sh` — clusterProfiler GO BP enrichment on genes with significant splicing changes, per timepoint
13. `12_trimmed_qc.sh` — FastQC + MultiQC (trimmed) — numbered 12 not 03 since scripts are numbered by write-order, added after the initial pipeline run, not before
14. `13_comprehensive_qc.sh` — aggregates raw + trimmed + salmon into one MultiQC report (raw/trimmed rows side by side for direct before/after comparison)

Figure/table-building scripts (`build_*.py`) run locally, not on the cluster — lightweight, faster iteration. All scripts follow the lab's minimal-commenting convention (cf. `orthogonal_validation/scripts/reference_scripts/scripts_from_github`) — short header, no debugging narrative in the comments.

**Environments:** cluster pipeline uses `biotools`. GO enrichment uses a cloned `go_entropy_hfd` (cloned from Project 3's `go_entropy` env specifically to avoid risking that project's reproducibility — never install directly into a shared env another project depends on).

## Key results

**W18 (18 weeks HFD) stands out consistently across every layer** — highest significant-gene count (1,155 vs 840-873 elsewhere), highest event count (driven by SE and A3 events specifically), and by far the strongest GO signal (496 significant terms vs 43-145 at other timepoints). Worth leading with this if presenting.

**GO theme progression across timepoints** (pattern-level observation, not a formal trajectory analysis, but lines up suggestively with Kiran's verbal description of HFD progression):
- Present at every timepoint: ribosome biogenesis/translation, nuclear-cytoplasmic transport/RNA export (the "core" splicing signature of HFD exposure)
- W8: DNA damage response, negative regulation of hepatocyte proliferation
- W18 (peak): nucleotide/purine biosynthesis, mitotic nuclear division (active proliferation)
- W24: viral/innate immune response
- W42: respiration/mitochondrial/hypoxia terms

**Splicing vs. Desmond's DE genes**: only 17-23% overlap at any timepoint — most splicing-significant genes show no significant expression change, confirming splicing captures a mostly distinct regulatory layer. (Two valid ways to express an overlap %, denominator matters — see script docstrings in `build_de_overlap_figures.py` before quoting a number: overlap÷splicing-sig gives 16.8-23.4%, overlap÷DE-sig gives 6.7-7.8%, both true, different questions.)

## Deliverables

- `figures/` — 9 figure families, PNG+SVG: 4 splicing scatter grids (+ "bis" minimal-style variants), event-type bar chart, gene-count-over-time dot plot, 5 GO dot plots (4 timepoints + union), DE-vs-splicing Venn (4-panel)
- `tables/` — 6 paired xlsx workbooks (`significant_events_by_type`, `significant_gene_count_over_time`, `splicing_scatter_significant_events`, `go_enrichment_significant_terms`, `de_splicing_overlap`, `sequencing_qc_summary`), each opens with a README sheet (thresholds, column defs, which figure it pairs with), scoped to significant/summarized results only (not raw PSI, which would be hundreds of thousands of rows)
- `notes/Splicing_HFD_Organoids.pptx` — 14-slide deck I built myself in PowerPoint, already emailed to Kiran and Desmond

## What's genuinely unfinished — read this before doing anything else here

**Desmond's reply (2026-08-12) to the emailed deck raised a methods question about slide 13/14 (the DE-vs-splicing Venn) — largely resolved 2026-08-21 by checking Kiran's team share in person over ethernet:**

1. ~~The DE gene lists used may be an earlier, more permissive version of his DE calls~~ — **resolved**: the share (`organoid/`) has a separate `Genes - Proteomics DEGs - logFC above 3/` folder distinct from `Genes - RNA Seq DEGs - logFC above 1/`. Desmond's "I think I raised it to 3" almost certainly refers to the proteomics threshold, not RNA-seq — the RNA-seq threshold used throughout is 1, matching what this project already used.
2. ~~A fuller DE dataframe may exist~~ — **checked, doesn't exist on the share**: every folder in `organoid/` on the share was opened and inspected — `Genes - RNA Seq DEGs - logFC above 1/` (same 7 filenames, same gene+cluster-only columns as what's already local, confirmed by content, not just row counts), `Genes - RNAseq + proteomics/` (gene-symbol intersect lists only, no logFC/p-value), and `raw data/` (`proteome/` = raw proteomics matrices + xlsx experiment files; `transcriptome/` and `transcriptome.zip` = identical, both just the raw fastqs). No file anywhere on the share has adj-p-value + directional logFC for RNA-seq. If it exists, Desmond has it locally and never exported it — would need to ask him directly, not worth another share-drive search.
3. Desmond offered to add Ensembl IDs to his future exports (currently symbol-only, which required a symbol-to-symbol match on my side instead of the more robust ID-based join). Still open, low priority.

**Next steps if resumed:**
- The threshold question (point 1) is resolved enough to not need a reply — the drafted email to Desmond can likely be dropped, or sent as a quick FYI/confirmation rather than a real ask.
- If the fuller dataframe is still wanted (for the dPSI-vs-log2FC correlation scatter, part of the reference paper's Fig 3C), ask Desmond directly — it is not sitting anywhere in this project or the team share.
- If a fuller dataframe does arrive, `build_de_overlap_figures.py` and the Venn figure/table/slide would need to be regenerated with real thresholds/direction.
- **Deprioritized, never started**, only pick up if there's spare time: entropy-based figures (paper's Fig 3E/F, Fig 2C-E — the "delta entropy" metric Project 3 used isn't reproduced here), Fig 2A-style UpSet plot.

## Known gotchas (PSMN-specific, apply to all cluster projects)

- `$HOME` quota is shared lab-group-wide — keep only scripts there. All data (raw fastq, reference, results) goes on `/scratch/Lake/$USER/`.
- Never transfer straight to `ssh.psmn.ens-lyon.fr` (the outer gateway) — it silently succeeds but lands in a ~5MB stub invisible from any real node. Two stray leftover scripts from this exact project were found there 2026-08-06 (`01_raw_qc.sh`, `00_download_reference.sh`) — confirms this isn't theoretical.
- MultiQC renders blank plots in the `biotools` env unless `python-kaleido` is installed via **conda-forge specifically** (`pip install kaleido` is unreliable on headless HPC).
- SUPPA2 (git-cloned into `biotools`) needs `scikit-learn`, `scipy<1.15` (newer scipy has a real regression breaking `scipy.special`), and `statsmodels` all installed together before it will run at all — see the `reference_biotools_env_gotchas` memory file for the full story if this needs re-debugging.
