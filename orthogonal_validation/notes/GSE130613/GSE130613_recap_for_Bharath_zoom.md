# GSE130613 (Gad Asher) — Recap for Zoom with Bharath

Prep notes for the MARS-seq call. All numbers below re-verified directly against current result files on 2026-07-14.

## Dataset
- GSE130613, Gad Asher lab, Weizmann. WT vs PerDKO liver, CT16-20, constant darkness (DD), n=4 each.
- WT: SRR9002567–570. PerDKO: SRR9002583–586. SRP194523.
- Confirmed by Asher (email, 2026-07-09): MARS-seq, "the common protocol in our sequencing facility."
- Reference: mm39/GRCm39. Custom transcriptome (Khushi's transcriptome_ext.fa) + custom GTF (transcriptome_productivity.gtf).

## Pipeline as run (all 4 stages complete)
1. **Salmon** — quasi-mapping, single-end (`-r`), library type `-l A` auto-detected as `-l U` (unstranded, compatible_fragment_ratio 1.0). Full trimming (Illumina adapter + polyA/T/G{12} tails) raised mapping rates: WT 82–88%, PerDKO 73–76% (from 44–70% untrimmed).
2. **tximport/DESeq2** (transcript-level, txOut=TRUE) + stageR two-stage testing.
3. **SUPPA2** — boundary `-b S` (strict), diffSplice empirical method, `-gc`.
4. **Figures/tables** — all generated, sent to Kiran/Bharath, in `DELIVERABLES/GSE130613_vs_Nanopore_v5_24:06/`.

## Current results
**DE (DESeq2 + stageR):** threshold padj<0.05 & |log2FC|>0.5.
- 447 significant genes (240 up in PerDKO, 207 down) — one representative transcript per gene (max |log2FC|).
- 769 significant transcript-level rows.

**Splicing (SUPPA, Illumina):** threshold |dPSI|>0.1 & p<0.3 (relaxed from 0.05 — low power at n=4, empirical method).
- 513 significant events across all 7 event types (A3, A5, AF, AL, MX, RI, SE).
- Overlap with Khushi's Nanopore diffSplice calls (exact event ID match): **32 events**.

**DE concordance vs Nanopore:**
- Gene level: 27 genes significant in both (background 9,850 shared genes).
- Transcript level: 16 transcripts significant in both, 11/16 consistent direction (background 12,966 shared transcripts).

## The open MARS-seq gap (why this Zoom exists)
- True MARS-seq puts the UMI/barcode on the forward read (R1); only R2 is the biological read to align.
- SRA metadata for these samples reports `layout=SINGLE` — only one read type was ever deposited publicly. Read headers checked directly (e.g. `@SRR9002584.1 1 length=75`) show no embedded barcode/UMI sequence.
- Net effect: the single-end approach isn't using the "wrong" read, but **no UMI-based deduplication has been done**, because there may be no UMI read available in the public data at all (not 100% confirmed — haven't checked GEO's raw supplementary files for a hidden multi-read structure).
- Bharath found **UTAP2** (Weizmann's own MARS-seq pipeline, Snakemake, gene counts + DE) as a candidate proper rework. Not yet implemented.
- Kiran's guidance so far: current DE/DESeq2 numbers above are fine to use as-is (possibly with tweaks Bharath mentioned verbally "this morning" — specifics not yet known, main reason for this call). Splicing/SUPPA numbers should not be over-interpreted given MARS-seq's 3' bias — junctions/internal exon structure aren't reliably captured.

## Questions to bring to the call
1. What are the specific modifications to the DESeq2/quantification step Bharath mentioned?
2. Does this replace the existing DESeq2 results (above), or refine them — do the 447/769 numbers need to be rerun?
3. The "ArrayExpress" comment in Kiran's email — GSE130613 came from SRA/GEO (SRP194523), not ArrayExpress. Is there additional/different raw data to source, or was that a generic offer?
4. Is UTAP2 worth starting now, or genuinely "after, if necessary" as Kiran said?
5. Confirm whether UMI/barcode reads exist anywhere in the public deposit (GEO supplementary files) before assuming UMI dedup is simply impossible.
