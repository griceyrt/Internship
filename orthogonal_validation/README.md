# Orthogonal Validation (Project 2) — handoff README

**Last updated:** 2026-08-21, Gricey's last internship day. This project has several parallel, semi-independent workstreams — read the "Status by workstream" section below to find what's active vs. closed before doing anything new.

## Goal

Validate splicing/DE events from the lab's own manuscript (Chikhaoui, Mamgain et al. — PER2/circadian clock, in revision) against independent public short-read Illumina RNA-seq datasets, since the lab's own primary data is Nanopore long-read. "Orthogonal validation" = does an independent platform (Illumina) + independent dataset see the same genes/splicing events as the lab's Nanopore data. Deliverable scope was deliberately kept simple per Gricey's own call: **figures + xlsx tables only**, no tracking database/Slides system.

## Key contacts

- **Kiran Padmanabhan** — PI (INSERM, IGFL/ENS Lyon). Sets priorities, corresponds with reviewers and dataset owners. kiran.padmanabhan@ens-lyon.fr
- **Bharath** — supervisor, pipeline architecture. Provided the original bash scripts under `scripts/reference_scripts/bharaths_code/`. Teaching most weeks, prefers email.
- **Khushi** — PhD student. Source of the Nanopore data/tables, GitHub analysis code, and rebuttal-writing coordination.
- **Lies** (he/him) — manuscript co-author, made the original PerKO isoform-usage figure.
- **Gad Asher** (Weizmann) — GSE130613/GSE171975 data generator.
- **Constanza** — wet lab, IGFL. Source of the IF colocalization imaging data (`Constanzas_files/`).

## Cluster access

PSMN cluster (ENS Lyon), username `grangel`. Full connection details, storage rules, and known gotchas are in `notes/` — but the authoritative, most up-to-date version of this lives in Claude's memory file `reference_cluster_connection.md` (applies to all 4 internship projects, not just this one). Quick version:
```
ssh-add ~/.ssh/id_ed25519
ssh -A -J grangel@ssh.psmn.ens-lyon.fr grangel@allo-psmn.psmn.ens-lyon.fr
ssh -A cl5218comp1
```
Never SSH/rsync straight to `ssh.psmn.ens-lyon.fr` — it has no real storage (~5MB), files silently vanish. Conda env on cluster: `biotools`.

## Status by workstream (read this first)

| Workstream | Dataset | Status | Where |
|---|---|---|---|
| GSE130613 pipeline | Gad Asher, WT vs PerDKO liver, MARS-seq 3' | **CLOSED 2026-07-20**, wrap-up confirmed by Kiran/Bharath/Khushi | `scripts/GSE130613/`, `results/GSE130613/` |
| GSE130613 second-dataset check | GSE171975 (CT20 window) | **Sent to team 2026-08-17**, awaiting Bharath's reply | `scripts/GSE171975/`, `results/GSE171975/` |
| GSE133398 (SRSF3 KO) | Nick Webster, TruSeq PE | **PAUSED 2026-08-03** by Kiran — 498 (paper) vs 84 (ours) sig splicing events gap unresolved | `scripts/GSE133398/`, `results/GSE133398/` |
| E-MTAB-9701 (Per2-KO HCT116) | Relógio lab, human colon cancer | **DROPPED 2026-07-20** — only 1 replicate/condition, not statistically usable | never downloaded past initial checks |
| Rebuttal Fig 1G / CT8 vs CT20 | Reviewer point (d) | **DONE, sent 2026-08-13** — Kiran satisfied with the result | `scripts/rhythmicity_analysis/`, `results/rhythmicity_analysis/rebuttal_ct8_vs_ct20/` |
| Reviewer 3 pt.1 — gene-level rhythmicity | Fig 1B replacement (501→773 genes) | **DONE, sent 2026-08-20**, composite figure fixed | `scripts/rhythmicity_analysis/`, `results/rhythmicity_analysis/` |
| Constanza IF colocalization | HepG2/mouse liver imaging | **DONE 2026-07-07** | `Constanzas_files/constanza_imaging/` |

**If picking this back up:** the two live threads waiting on a reply are GSE171975 (Bharath) and GSE133398 (Kiran — the 498 vs 84 gap). Everything else is closed unless a new reviewer question reopens it.

---

## Datasets

### GSE130613 (primary, CLOSED)
Gad Asher lab, PNAS 2020 (Manella & Asher). WT vs PerDKO (Per1/Per2 double KO) mouse liver, constant darkness, CT16-20 window. **MARS-seq 3'-end** — this drove the whole pipeline rework story below. n=4 WT (SRR9002567-570) vs n=4 PerDKO (SRR9002583-586), SRP194523.

**Final adopted pipeline (as of 2026-07-19, after a long Salmon→STAR investigation):** STAR (splice-aware, mm39, `transcriptome_productivity.gtf` — original/non-extended) → htseq-count (gene-level, deterministic) → DESeq2 (plain Wald test, no tximport/stageR needed since htseq-count is already gene-level). **390 significant genes** (238 up, 152 down in PerDKO), padj<0.05 & |log2FC|>0.5.

Why not Salmon: Salmon's EM-based transcript-level quantification produced an implausibly wide log2FC spread (±28) that didn't reproduce in Asher's own GEO-deposited official counts (±4-6) or the lab's Nanopore data (±7-8). Root-caused to Salmon's ambiguous-read handling on this large custom transcriptome, not UMI/PCR-duplication (that was investigated and ruled out — this bulk MARS-seq variant uses T7 IVT linear amplification, not heavy PCR, so no UMI was ever needed).

Splicing: **rMATS-turbo** replaced SUPPA2/Salmon for the same reason (43 sig events at p<0.3, 8 at FDR<0.05, vs SUPPA's original 513 — expected big drop, method AND data both changed). **The 3'-bias caveat is permanent and unresolved**: MARS-seq only sequences the 3' end, so internal splice junctions are structurally undercaptured no matter which tool is used. Kiran/Bharath's explicit instruction: splicing results from this dataset should never be over-interpreted.

Full investigation trail (Salmon→STAR, UMI dead-end, rMATS tool selection) is documented in Claude's memory file `gse130613_marsseq_rework.md` if the reasoning needs to be reconstructed later.

### GSE171975 (second Asher dataset, ACTIVE)
Same Asher lab, PLOS Biology 2021, same WT/PerDKO/constant-darkness design but a longer 2-day time course — only the CT20 window was used (matching GSE130613). Requested by Kiran 2026-08-17 as a second independent check. WT n=4, PerDKO n=3 (unbalanced, confirmed fine by Kiran — DESeq2 handles it natively). Same STAR+htseq+DESeq2 pipeline as GSE130613.

**Result: 1,759 significant genes (770 up, 989 down)** — 10.9% of tested genes vs GSE130613's 2.4%, but Nanopore validation rate is *lower* (8.3% vs GSE130613's 13.7%). Kiran's reading: this supports the reviewer's original concern that Illumina over-calls DE genes relative to Nanopore, likely via amplification artifacts. Sent to the team, **awaiting Bharath's reply**.

### GSE133398 (SRSF3 KO, PAUSED)
Nick Webster lab, standard TruSeq stranded PE (not MARS-seq — no 3'-bias caveat applies). CRE-KO (SRSF3) vs GFP control hepatocytes, n=4 each. Ran the full standard Salmon (mm39) → DESeq2 → SUPPA2 pipeline (unlike GSE130613, Salmon was fine here since it's standard PE data). Two DE options built (Option A gene-level: 1,595 sig genes; Option B transcript-level+stageR: 998 genes/950 transcripts) since Bharath wanted both to choose from.

**Paused 2026-08-03**: the original paper reports 498 SRSF3-dependent splicing events; this pipeline found only 84 (p<0.05, SUPPA2). Kiran is thinking about how to present/reconcile this gap — likely candidates are a different splicing tool in the original paper (rMATS vs SUPPA2 disagree substantially) or a looser significance threshold there. **No resolution yet — check with Kiran before resuming.**

### E-MTAB-9701 (DROPPED)
Reviewer originally asked for this specific dataset (Per2-KO HCT116 human colon cancer cells) instead of a mouse liver dataset — considered weak science by Kiran (wrong organism/tissue) but investigated per Bharath's directive to "at least try." Found only 1 biological replicate per condition (the "2 replicates" Kiran remembered were actually R1/R2 paired-end file rows, not real replicates) — statistically unusable for DESeq2 (no dispersion estimate possible at n=1). Both Bharath and Kiran confirmed dropping it, using the 1-replicate finding itself as the reviewer response. **Never downloaded past a few test files.**

---

## Rebuttal workstreams (reviewer-driven, separate from the dataset validation above)

These live in `scripts/rhythmicity_analysis/` and use the lab's own Nanopore/Salmon data (`data/salmon_2024-09-25/`, `data/tables/`), not the public GEO datasets above.

### Reviewer point (d) — CT8 vs CT20 (DONE, sent 2026-08-13)
Tests whether "more splicing-significant than DE-significant genes" is a genotype-clock-disruption-specific pattern (as claimed in the manuscript) or just generic to any two-state comparison. Ran the same WT day-vs-night (CT8 vs CT20) comparison and computed the same ratio. Result: CT8-vs-CT20 shows an *even bigger* splicing:DE disparity (10.12x) than the WT-vs-PerKO comparison (3.64x) — initially looked like it undercut the manuscript's argument, but **Kiran's actual take was positive**: low absolute counts (23 DE / 19 splicing genes) from a same-genotype comparison support that the much larger PerKO numbers reflect real clock-disruption biology, not a generic artifact. Fully resolved, no action needed.

### Reviewer 3 point 1 — gene-level rhythmicity (DONE, sent 2026-08-20)
Manuscript's original Fig 1B tested rhythmicity at transcript-level (501 genes via cosinor+DESeq2 LRT+stageR). Reviewer asked for a proper gene-level collapse. Rebuilt at gene level (no stageR needed — one test per gene): **773 significant rhythmic genes** (7.1% of testable genes, up from 3.9%), still below literature's 15-25% but a real improvement. 458 of the original 501 genes carried over, 315 newly gained, 43 lost — the "isoform fragmentation dilutes power" mechanism is real, not just a looser threshold (all 7 core clock genes present with strong significance as a sanity check).

Kiran's hand-drawn "proposed composite" figure had a factual issue (showed 70 "masked" genes — isoform-significant but not gene-significant — as a nested subset of the 773-gene panel; they're actually the complementary set, not nested). Flagged and fixed — final figure uses two stacked, non-nested blocks. Sent 2026-08-20.

Key files: `gene_rhythmicity_significant.xlsx` (773 genes, 4 tabs incl. lost/masked genes), `gene_rhythmicity_heatmap.svg`, `proposed_composite_v2.png`. The DTU (differential transcript usage) list Kiran asked for was initially mislabeled — the correct file is `data/tables/DTU_genes.csv` (14 unique genes), not `Table5-rhythmic_AS_events.xlsx` (228 genes, a different/larger category — rhythmic AS events, not DTU). Confirmed by Khushi.

**Still open, never started:** a supplementary polar-plot figure Kiran mentioned as "not critical if you can't find the code," and category/pathway annotations on the Fig 1G Venn (only exists in an xlsx tab, never added to the SVG itself).

---

## Constanza IF colocalization (DONE 2026-07-07)

`Constanzas_files/constanza_imaging/` — turns Constanza's Fiji/Coloc2 PER2-SRSF3 colocalization output into publication figures supporting the manuscript's IP result (Fig 4d). Two scripts: `tidy_data.py` (parses the 5-sheet source xlsx, drops confirmed QC-excluded cells) → `make_figures.py` (4 final figures: 3 HepG2 endogenous-vs-overexpression pairs + 1 mouse liver nuclei strip plot). A phospho-mutant screen panel (U2OS-HA-SRSF3-Mutants) was explicitly not needed per Constanza — data is tidied but no figure was built for it.

---

## Folder structure

```
orthogonal_validation/
├── data/                        reference files, Salmon quant, Khushi's Nanopore tables
│   ├── salmon_2024-09-25/       raw Salmon quant for the rhythmicity/rebuttal work (16 samples)
│   ├── nanopore_suppa/          Khushi's diffSplice output (WT-CT20 vs PerKO-CT20)
│   ├── tables/                  Khushi's Table2/3/5, DTU_genes.csv
│   ├── transcriptome_ext_index/ Salmon index for the custom transcriptome
│   └── files_from_khushi/
├── raw_data/                    downloaded SRA fastq (GSE130613 + GSE133398 accessions)
├── scripts/
│   ├── GSE130613/               full Salmon→STAR pipeline + rMATS + figures/tables
│   ├── GSE133398/               Salmon + SUPPA2 pipeline (Option A/B) + figures/tables
│   ├── GSE171975/               STAR+htseq pipeline (reuses GSE130613's index/GTF)
│   ├── rhythmicity_analysis/    rebuttal-specific scripts (Fig 1G, Fig 1B, CT8vsCT20)
│   └── reference_scripts/       Bharath's original bash scripts + Khushi's GitHub repo
├── results/                     per-dataset outputs, one subfolder per GSE + rhythmicity_analysis/
│   └── archive/                 superseded runs, kept for reference only
├── Constanzas_files/             IF colocalization figures (separate wet-lab-data subproject)
├── notes/                        workflow diagrams, deliverables snapshots, cluster connection notes
├── papers/                       source manuscripts/papers referenced throughout
└── meta/                         accession lists
```

Each GSE folder's script naming is dataset-prefixed where relevant (e.g. `GSE133398_*.R`) to avoid collisions when scripts get copied to the cluster.

## Known gotchas

- **iCloud sync**: this whole `Internship/` folder syncs via iCloud with Optimize Storage on — causes false git deletions, occasional file-read errors, and can lock `.git/index.lock`. Fix is Finder "Download Now" on the folder, or delete the lock manually.
- **PSMN gateway clutter**: any transfer that skips `allo-psmn`/`cl5218comp1` and hits `ssh.psmn.ens-lyon.fr` directly silently lands in a ~5MB stub, invisible from any real node. Worth a periodic `find ~ -type f` check there.
- **`$HOME` quota is shared across the whole igfl lab group** on PSMN — keep only scripts there, put bulk data (fastq, indices, BAMs) on `/scratch/Lake/$USER/` instead (no quota, but also no backup).
- **NaN gene-name bug**: `transcriptome_productivity.gtf` is missing `gene_name` for ~20% of genes (custom/non-canonical entries); a plain dict `.get(key, default)` won't catch a key that exists with a stored NaN value — always check `pd.isna()` explicitly when labeling figures by gene name from this GTF.
- **ggplot squishing**: to compress a heatmap panel's aspect ratio, use `theme(aspect.ratio=...)` in R, never post-hoc non-uniform PNG resizing (destroys text legibility).

## What's genuinely unfinished

1. **GSE133398**: 498-vs-84 splicing event gap, paused pending Kiran's direction on how to present it.
2. **GSE171975**: sent to the team, no reply from Bharath yet as of 2026-08-17.
3. **`GSE_datasets_evaluated.docx`** (a tracking table covering all 4 GEO datasets ever evaluated) has a stale row for GSE171975 still saying "not used" — needs updating, blocked on Gricey sharing the actual file (only a screenshot exists).
4. Reviewer 3 pt.1's optional polar-plot figure and Fig 1G category/pathway Venn annotations — never started, explicitly low priority.
