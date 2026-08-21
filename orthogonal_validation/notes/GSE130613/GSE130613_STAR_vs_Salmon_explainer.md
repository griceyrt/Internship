# GSE130613 MARS-seq: why we're testing STAR, and what it actually tests

Reading material for while things run on the cluster. Covers: the problem we're chasing, how Salmon and STAR fundamentally differ, why that difference matters for *this* dataset specifically, and the full decision tree for the investigation.

---

## 1. The problem, recapped

Our Illumina/MARS-seq pipeline (Salmon + custom transcriptome + DESeq2 + stageR) produces a WT-vs-PerDKO differential expression comparison with a suspiciously huge spread: log2 fold-changes ranging from about **-28 to +28**.

Two independent reference points suggest this is too wide to be real biology:

- **GEO's own official gene counts** for this exact comparison (Asher's team's own processing) give log2FC **-4.13 to +6.03** on the same biological samples.
- **Nanopore long-read data** (Khushi's parallel platform, same biological question) gives log2FC **-7.66 to +6.78** for the same set of genes/transcripts.

Both independent sources land in a tight, believable range. Only our own Illumina pipeline blows out to ±28. Something about *how we're processing the data* — not the underlying biology — is almost certainly responsible.

## 2. What we've ruled out so far

**Illumina vs. Nanopore as platforms:** ruled out. The GEO-official-counts check reused the *same* Illumina reads as our pipeline, just processed differently (their pipeline, not ours), and got a tight range. So it's not "Illumina is a noisier platform than Nanopore," it's something about our specific processing.

**PCR duplication exists, is a plausible contributor:** confirmed but not fully quantified. FastQC on our raw reads shows 48-54% sequence duplication, with meaningful fractions of reads duplicated 100-1,000+ times. MARS-seq needs heavy PCR amplification for such low RNA input, so this tracks. But we can't yet say *how much* this alone explains, since true UMI-based deduplication isn't possible for this dataset (see box below).

**Annotation boundaries (3'UTR extension):** ruled out as *the* fix. We extended every transcript's terminal exon (1000bp toward 5', 100bp toward 3', matching Asher's own published method exactly) so reads landing just past the current annotated end could still map. Result: -27.59 to +23.79, barely different from the original -28.08 to +28.01. This was the "cheapest" fix to test and it didn't work.

> **Why can't we just do UMI deduplication and be done with it?**
> True MARS-seq puts a random UMI (unique molecular identifier) on one read and the biological sequence on the other, so PCR duplicates — many reads that are all copies of the same original RNA molecule — can be collapsed back down to one count. But the raw data deposited on SRA for this dataset is **single-end only** (`layout=SINGLE` in the metadata, and the read headers have no barcode/UMI sequence in them). The UMI read appears to have never been publicly deposited, only the demultiplexed biological read. Both candidate MARS-seq-aware pipelines we researched (UTAP2, nf-core/marsseq) *require* that missing UMI read as input — there's no fallback mode for single-end, UMI-less data in either tool. So "just deduplicate properly" may not be achievable at all with the public data, regardless of which pipeline processes it.

## 3. How Salmon actually works (and why it might matter here)

Salmon doesn't align reads in the traditional sense. It builds an index off all reference sequences (in our case: 155,263 transcripts + genome decoy) using a **de Bruijn graph** — the same graph-theory structure used in genome assembly, going back to Euler's 1736 Königsberg bridges problem. Every possible 31-base chunk (k-mer) of every reference sequence becomes a node; nodes get connected if they overlap by 30 bases in some real sequence. This lets Salmon do **quasi-mapping**: instead of full alignment, it just checks which paths through this graph a read's k-mers are consistent with.

Here's the part that matters for us: a read can be consistent with *multiple* transcripts at once — most obviously when two isoforms of the same gene share an exon, but also when two barely-related transcripts happen to share a stretch of sequence (a real risk with our transcriptome, since it merges standard annotation with ~5,840 novel FLAIR-discovered transcripts that may be fragmentary or partially overlapping with known ones). Salmon resolves this ambiguity with an **EM (Expectation-Maximization) algorithm**: it estimates each transcript's abundance, uses those estimates to decide how to probabilistically split an ambiguous read's "weight" across the competing transcripts, re-estimates abundance, and repeats until it converges.

This is elegant and usually works well — it's *why* Salmon can quantify isoforms at all, rather than just genes. But it means every ambiguous read's contribution to a transcript's count is a probability-weighted guess, not a hard assignment. If a cluster of PCR-duplicated reads happens to map ambiguously, small differences in how those duplicates distribute between samples (pure PCR noise, not biology) could get amplified by the EM step into apparent expression differences. This is a real, plausible mechanism for exactly the kind of noisy, wide-spread log2FC we're seeing, on top of whatever raw duplication is doing.

## 4. How STAR + htseq-count works instead

STAR is a real aligner: it finds the actual best-scoring position(s) in the **genome** (not just the transcriptome) for each read, using a suffix-array-based seed-and-extend search, and is splice-aware (it knows to look for reads that span exon-exon junctions). Output is a standard BAM file with real genomic coordinates.

`htseq-count` then takes that BAM plus a GTF and does simple, rule-based counting: for each read, which gene's annotated region does it overlap? If it overlaps exactly one gene, count it there. If it overlaps zero or more-than-one gene ambiguously (default "union" mode), the read gets **discarded**, not probabilistically split.

This is what I mean by "deterministic": each read either gets a hard count toward one gene, or gets thrown out. No EM step, no probability redistribution, no dependence on relative abundance estimates from other samples. It's a blunter instrument than Salmon in some ways (you lose isoform-level resolution, and legitimately ambiguous reads are just dropped rather than intelligently allocated), but that bluntness means it can't have the specific failure mode described above.

This also happens to be **exactly the method Asher's own team used** to generate their official GEO counts (confirmed from their SI Methods: cutadapt → STAR against mm10 → htseq-count union mode, RefSeq annotation extended the same way we just tested). So running STAR + htseq-count with our extended GTF is, in effect, a full replication of their original pipeline — the cleanest possible test of whether the quantification *method* itself (not annotation, not literal UMI dedup) was the missing piece.

## 5. What the result will tell us

| Outcome | What it means |
|---|---|
| STAR + htseq-count gives a tight range (close to GEO's -4 to +6) | Salmon's EM-based transcript-level quantification was the real driver of the inflated spread. Fix: switch to gene-level STAR+htseq-count for this dataset (accepting the loss of isoform-level resolution). |
| STAR + htseq-count is *still* wide | Rules out both annotation boundaries and quantification method. Points squarely at PCR duplication itself as the remaining explanation — and since true UMI dedup isn't possible with this public data, this may be a hard limitation of the dataset rather than something fixable in our pipeline. |

## 6. Full decision tree so far

```
GSE130613 MARS-seq — inflated log2FC spread investigation
(Original Salmon pipeline: log2FC -28.08 to +28.01, 447 sig genes)
(Reference — GEO official counts: -4.13 to +6.03 | Nanopore: -7.66 to +6.78)
│
├─ Q: Is this an Illumina-vs-Nanopore platform difference?
│    └─ NO — GEO's own Illumina-based counts are tight too. Ruled out.
│
├─ STEP 1 (least effort): 3'UTR boundary extension
│    │  (extend GTF 1000bp 5' / 100bp 3', matching Asher's method,
│    │   rebuild Salmon index, requantify, same EM-based quantification)
│    │
│    └─ RESULT: log2FC -27.59 to +23.79, 395 sig genes
│         └─ FAILED to meaningfully tighten the spread.
│              └─ Annotation boundaries ruled out as the main driver.
│
├─ STEP 2 (current): STAR + htseq-count
│    │  (genome alignment + deterministic gene-level counting,
│    │   full replication of Asher's original method, same extended GTF,
│    │   still no UMI dedup — not possible with this public data)
│    │
│    ├─ IF tight range (~GEO's -4/+6) ──► Salmon's EM quantification was
│    │                                     the driver. Adopt STAR+htseq-count
│    │                                     for this dataset. DONE.
│    │
│    └─ IF still wide ──► Quantification method also ruled out.
│         │
│         ├─ STEP 3: UTAP2 (Weizmann's own MARS-seq pipeline)
│         │    └─ LIKELY BLOCKED: requires the UMI/barcode read as
│         │       separate input, which was never deposited publicly
│         │       for this dataset (layout=SINGLE, no barcode in
│         │       read headers). May not be runnable at all.
│         │
│         └─ STEP 4: zUMIs (fallback if UTAP2 doesn't work)
│              └─ SAME LIKELY BLOCKER — also needs the UMI read.
│                   │
│                   └─ IF both blocked: PCR duplication may be an
│                      inherent, uncorrectable limitation of this
│                      specific public data deposit. Team would need to
│                      decide how to caveat/handle this for the
│                      manuscript rather than "fix" it further.
```

## 7. Where things stand right now

- Step 1 (3'UTR extension): done, negative result, reported to Bharath.
- Step 2 (STAR + htseq-count): about to start. Need to confirm `STAR` and `htseq-count` are available in the `biotools` cluster environment before building the pipeline (index build, alignment × 8 samples, counting, then a small DESeq2 script adaptation since `htseq-count` output isn't the same format as Salmon's `quant.sf`).
- Steps 3-4 (UTAP2/zUMIs): not started, likely to hit the missing-UMI-read blocker described above if we get there.
