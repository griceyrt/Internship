# Fmr1 ONT Splicing and PolyA (Project 3)

**Dataset:** GSE188840, WT vs Fmr1 KO mouse brain cortex, ONT PromethION direct RNA-seq (SQK-RNA002), 3 replicates per genotype (`WT_Rep1-3`, `KO_Rep1-3`). This is Kiran/Lies's own dataset, published in Shin, Paek, Chikhaoui et al. 2022, *RNA* 28:756-765 ("Oppositional poly(A) tail length regulation by FMRP and CPEB1", PMID 35217597). See `papers/`.

**Goal (Kiran's 3-part ask):** (1) gene-level differential expression, (2) SUPPA splicing landscape/dPSI across all 7 event types, (3) overlap between DE and dPSI genes. A harder follow-up question (poly(A) tail length vs splicing correlation, especially alternative-last-exon events) was added once fast5 raw signal became available.

**Handed off:** 2026-08-21, my last internship day. See "Where this stands" below for exactly what's done vs. still open. The `06_de_methodology` investigation (section 7 below) was genuinely closed out on this same last day, via a reply to an unanswered email from Lies. See that section for the final numbers.

---

## Directory structure

```
fmr1_polya_splicing/
├── scripts/            cluster pipeline (sbatch scripts)
│   ├── 00_reference_download/       genome/GTF/transcriptome + fastq download
│   ├── 01_alignment_quantification/ Minimap2 + oarfish
│   ├── 02_deseq2_suppa_core/        DESeq2 (gene-level) + SUPPA (7 event types)
│   ├── 03_isoform_entropy/          Shannon entropy of isoform usage, WT vs KO
│   ├── 04_go_category_analysis/     GO/category tests on entropy shift
│   ├── 05_overlap1/                 DESeq2 ∩ SUPPA overlap
│   ├── 06_de_methodology/           DE methodology investigation (see below)
│   ├── 07_replicate_qc/             PCA/distance/Cook's distance outlier check
│   ├── 08_fast5_polya_prep/         fast5 download + migration to scratch
│   ├── 09_igv_validation/           Pik3cd IGV sanity check
│   └── 10_polya_length/             Nanopolish poly(A) length (IN PROGRESS)
├── plot_scripts/        local plotting scripts, same phase numbers as
│                         scripts/ where applicable
├── data/reference/       small local reference files only (big genome/GTF/
│                         transcriptome files stay cluster-only, no reason to
│                         duplicate hundreds of MB locally)
├── results/              pulled-down TSVs, stats tables, some figures
├── figures/              hand-picked FINAL figures for the pptx deck (not a
│                         comprehensive output archive, see results/ for that)
├── presentations/        workflow diagrams (Workflow_fmr1_updated.pdf/.key),
│                         entropy pptx deck
├── notes/                cluster connection how-to, IGV screenshots, cluster
│                         cleanup checklist, reference images from Kiran
├── papers/                reference PDFs (the source paper + comparison papers)
└── meta/                 GSE188840_Acc_List.txt (SRA accession → sample mapping)
```

**IMPORTANT: the `scripts/` numbering is BUILD ORDER, not the scientific workflow phase numbering** used in `presentations/Workflow_fmr1_updated.pdf` (my own planning diagram, Phase 0-7). They were deliberately kept separate rather than renumbered to match. Rough mapping if you need to cross-reference:

| Workflow.pdf phase | scripts/ folder |
|---|---|
| Phase 0 (reference) + Phase 1 (raw data) | `00_reference_download/` |
| Phase 2 (core quantification) | `01_alignment_quantification/` + `02_deseq2_suppa_core/` |
| Phase 3 (splicing/SUPPA) | `02_deseq2_suppa_core/` |
| Phase 3.1 (isoform landscape) | `03_isoform_entropy/` |
| Phase 4 (Overlap #1) | `05_overlap1/` |
| Phase 5 (custom transcriptome) | **not built yet** |
| Phase 6 (poly(A) length) | `10_polya_length/` (in progress) |
| Phase 7 (Overlap #2) | **not built yet** |

`04_go_category_analysis/`, `06_de_methodology/`, `07_replicate_qc/`, `09_igv_validation/` have no counterpart in the original diagram. They're follow-up investigative work that came out of Kiran's requests mid-project, not part of the original plan.

---

## Cluster access

Username `grangel`. Real home only starts at `allo-psmn`, NOT the external gateway (`ssh.psmn.ens-lyon.fr` has ~5MB max). Working node: `cl5218comp1`.

```
scp -o "ProxyJump=grangel@ssh.psmn.ens-lyon.fr,grangel@allo-psmn.psmn.ens-lyon.fr" \
    <local_path> grangel@cl5218comp1:<remote_path>
```

Shared conda env `biotools` has almost everything; a separate env **`go_entropy`** was created specifically for `org.Mm.eg.db`/`GO.db` (installing them into `biotools` fails the solver, an old R version pin conflict from deseq2/tximport). Nanopolish was newly added to `biotools` (see `10_polya_length/`, install command is in that script's comments if missing).

**Storage tiers** (check current usage with `df -h <path>` on `cl5218comp1`):
- `$HOME`: 41 TiB, SHARED across the whole `igfl` lab group, not per-user. Scripts/code only.
- `/scratch/Lake`: 520 TiB, no quota, but NO backup/snapshots. Large working data (BAMs, fast5, etc.) lives here.
- `/Xnfs/igfldb`: a lab-group NFS mount that IT staff said gets periodically wiped. Use only as brief staging, never final storage.

---

## Pipeline walkthrough (build order)

1. **`00_reference_download/`**: mm39 (GRCm39) genome fasta, GTF, stock cDNA transcriptome from Ensembl (uses the REST API to find the current release number, not a hardcoded one). Also pulls the 6 fastq samples via prefetch+fasterq-dump. **Genome build is mm39**, not the original 2022 paper's mm10. Kiran's explicit call to use the current build.

2. **`01_alignment_quantification/`**: Minimap2 (`-ax map-ont -N 100`) against the **transcriptome** (not genome), sorted+indexed with samtools. Then **oarfish** (not Salmon) for quantification. Salmon's own `--ont` flag redirects ONT users to oarfish, confirmed directly in its help text.

3. **`02_deseq2_suppa_core/`**: tximport (transcript-level → SUPPA TPM export; gene-level → DESeq2). **DESeq2 runs at gene-level** (Lies's call, simplifies Overlap #1; doesn't affect SUPPA since SUPPA's TPM comes from the quantifier directly). SUPPA `generateEvents` uses `-b V -t 1` ("variable boundary, 1nt threshold", this is what Kiran meant by "the +1 variable position"), all 7 event types (SE, A3, A5, MX, RI, AF, AL) tested **separately**, never pooled.

4. **`03_isoform_entropy/`**: per-gene Shannon entropy of isoform usage proportions (H = -Σp·log2p), WT vs KO. Ranked by entropy shift, tested at 3 TPM noise-filter thresholds (10/20/50, min50TPM primary per Lies's statistical rationale: low read depth makes proportion estimates unstable, not the entropy math itself). **Pik3cd and Rfx3** are the only genes in the top-15 at all 3 thresholds.

5. **`04_go_category_analysis/`**: two separate analyses, don't confuse them: (a) `entropy_by_gene_category.R` replicates PMC10099729's method (Wilcoxon test on entropy_diff, housekeeping-vs-other and high-vs-low-expression splits); (b) `entropy_go_enrichment_dotplot.R` replicates **Figure 4E of Kiran/Lies's own unpublished manuscript** (Chikhaoui/Mamgain/Seki et al., bioRxiv 2024.12.23.630108, this is what Kiran meant by "the Nature Comms manuscript", confirmed by matching the figure legend text exactly), a per-GO-term paired Wilcoxon test, NOT hypergeometric enrichment (the first version of this script got that wrong and was rewritten).

6. **`05_overlap1/`**: genes with both DESeq2 significance (padj<0.05) AND SUPPA significance (raw p<0.05) per event type. **142 DESeq2-significant genes** total; overlap ranges 0 (AL) to 8 (SE) genes per event type.

7. **`06_de_methodology/`**: a deep investigation triggered by Kiran asking why the DE volcano plot looked asymmetric (126 up/16 down) compared to the source paper's Figure 3F (which looks roughly symmetric). Built a 12-test grid in two rounds:
   - **Tests 1-8** (`deseq2_gene_level_2plus3_test.R`, `deseq2_isoform_level_test.R`): gene-level vs isoform-level DESeq2 × padj vs raw-p × with/without a post-hoc 2.5-fold filter. Test 8 (isoform, raw p, 2.5-fold) looked like the closest methodological match at first but gave **1170 vs the paper's stated ~100 (~12x)**. Lies then flagged (email 2026-08-10) that the paper's plot, despite being axis-labeled "p-value," may actually be padj. Recomputing with padj instead (test 7) dropped it to **61**, closing most of the gap in one move.
   - **Tests 9-12** (`deseq2_lfcthreshold_test.R`): tested whether DESeq2's proper `lfcThreshold`-based Wald test (against the fold-change boundary directly, via `altHypothesis="greaterAbs"`), as opposed to tests 1-8's post-hoc filtering (test against zero, then discard anything under the fold-change cutoff), closes the remaining gap. **Result: it made every combination smaller, never closer to 100** (isoform+padj: 61 → 18). This rules out the testing-methodology question as the explanation.
   - **Conclusion: closed out 2026-08-21, this is as far as this investigation gets.** Test 7 (isoform-level, padj<0.05, post-hoc 2.5-fold filter, 61 total) is the best match found across all 12 tests, real, about 40% short of the paper's ~100, not exactly reproduced. A separate check (`check_expressed_transcripts_10reads.R`) found the total expressed-transcript universe is close to theirs (33,914 vs their stated 30,975), ruling out a gross annotation-size mismatch too. **Remaining, untested candidates: genome build (mm10 vs mm39) and quantifier (their Salmon v1.5.2 alignment-mode vs oarfish)**. Both could shift which individual transcripts cross the significance line without changing the overall universe size, but testing either means rebuilding part of the pipeline, not a quick parameter swap.

   Full 12-test table saved at `results/de_lfcthreshold_test/all_12_tests_summary.txt`. I drafted a reply to Lies covering all of this (padj fix, expressed-transcript check, lfcThreshold negative result, final conclusion), but hadn't confirmed sending it. Check whether it went out before assuming this thread is closed on his end too.

8. **`07_replicate_qc/`**: Kiran said the paper excluded a noisy WT replicate (2 WT + 3 KO vs 3+3 here). Checked with PCA, sample-distance, and Cook's distance: **WT_Rep1 is confirmed an outlier on all 3 metrics** (large PCA separation on PC1, highest distance to its own replicate group, highest median Cook's distance). Excluding it improved the padj-based asymmetry (7.9:1 → 4.5:1) but not the raw-p-based one, so this explains *part* of the paper-comparison gap, not all of it.

9. **`08_fast5_polya_prep/`**: the 3.4TB fast5 raw signal download from Kiran's Japan collaborator (`kero.hgc.jp`). **This took many days and several real wget bugs to resolve** (all documented in the script's own comments): the server doesn't redirect no-slash directory URLs, doesn't handle HTTP Range requests, and never sends Last-Modified headers, so `-N`/`-c` (wget's normal resume flags) don't work here at all. `-nc` plus manual truncated-file cleanup was the real fix. **Actual final total: 3,701,803,782,411 bytes = 3.70TB, verified 3 independent ways, slightly EXCEEDS the ~3.4TB the email stated.** Don't be alarmed if `du -sh` looks short of that, it was a TiB-vs-TB display rounding issue, not missing data. Migrated to `/scratch/Lake/grangel/fmr1_polya_splicing/fast5/`.

10. **`09_igv_validation/`**: extracts Pik3cd reads from the transcriptome BAM, re-aligns to the genome (splice-aware), for IGV/Sashimi visualization. Confirmed 36 annotated isoforms are real (not a pipeline artifact), sent to Kiran, he confirmed.

11. **`10_polya_length/`, IN PROGRESS, not finished, actively running when the internship ended.** Nanopolish index + `nanopolish polya` per sample, using **fast5_pass reads only** and the EXISTING transcriptome-aligned BAMs from step 1 (no realignment needed). Sample-name mapping between the fast5 folders (`210126_DrLies_WT1_CTX` etc.) and this project's own naming (`WT_Rep1` etc.) is assumed 1:1 by replicate number, **not independently verified**, though `nanopolish index` would reveal a mismatch as very few successfully-indexed reads (watch for the script's own warning if any sample shows under 1000 reads). **This job needs multiple resubmissions** (same pattern as the fast5 download) since Nanopolish's signal-level analysis is slow. Both the indexing and polya output steps are skip-if-already-done, safe to resubmit.

    **Exact state at handoff (2026-08-21):** job `15614055` submitted, running on `c6420node091`. Deliberately left running rather than cancelled. 6 samples × 2 steps each is a genuinely multi-hour job that could not finish in the time remaining, but every completed step is saved and skipped on resubmission, so nothing is lost by leaving it to run as far as it gets (or die at the 24h SLURM limit) unattended. 30 minutes in, it was still on step 1/2 (`nanopolish index`) for sample 1/6 (`KO_Rep1`), expect it to not have gotten far. **To resume:** `cd ~/fmr1_polya_splicing/scripts/10_polya_length && sbatch nanopolish_polya.sh`. It will skip any sample/step already completed and continue with the rest. Check `squeue -u grangel` and the `*_nanopolish_polya.log` file first to see how far it actually got.

---

## NOT YET BUILT, the real remaining work

- **Custom transcriptome (Workflow.pdf Phase 5)**: IsoQuant (Kiran's tool choice over the lab's existing FLAIR pipeline) for novel isoform discovery, to build an expanded GTF. Needed only for a *refined* SUPPA re-run feeding Overlap #2, not needed for anything else already done.
- **Poly(A) aggregation**: once `10_polya_length/` produces real per-read tables, a binned aggregation script (per Lies's endorsement of "binned comparison" as the aggregation method) needs designing against actual data. Don't guess bin boundaries ahead of seeing the real distribution.
- **Overlap #2 (Workflow.pdf Phase 7, "Kiran's harder question")**: poly(A) length compared between genes with vs. without a significant splicing event, per event type (all 7, matching Overlap #1's convention), testing whether missplicing → differential polyadenylation. Needs both the custom-transcriptome refined SUPPA output AND the poly(A) aggregation above.

---

## Known gotchas worth knowing before touching this again

- **`ggsave()`'s SVG writer needs `svglite`**, not installed by default. Either install it or use base R's `grDevices::svg()` device instead (several scripts in this project already do the latter to avoid touching shared envs).

## Key contacts

- **Kiran Padmanabhan**, PI, INSERM/IGFL Lyon. Sets priorities, corresponds with the Japan collaborator and dataset originators.
- **Lies**, co-author on the source paper, made the original PerKO isoform-usage figure that inspired the entropy work here, reviews methodology choices.
- **Bharath**, pipeline architecture, teaching-heavy schedule so email is best.
- **Khushi**, PhD student, R/normalisation/SUPPA-format questions.
