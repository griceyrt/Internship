# Orthogonal Validation — Project 2

> **One-line summary:** We used public Illumina short-read RNA-seq data to independently confirm splicing changes found in the lab's Nanopore long-read data, in mice without a functional circadian clock.

---

## What is "orthogonal validation"?

Orthogonal means "from a completely different angle." In science, when you find a result using one technology, you want to check if the same result appears using a completely different technology. If it does, the finding is much more trustworthy — it's unlikely that two different methods produce the same artefact.

In our case:
- **Kiran's lab** found differential splicing events in PerDKO vs WT mouse liver using **Nanopore long-read sequencing** (the new, expensive, reads-full-transcripts technology)
- **We** tried to confirm the same events using **Illumina short-read sequencing** (the older, cheaper, standard technology)

If the same splicing events appear in both datasets, independently, that is orthogonal validation.

---

## The Biology

### The circadian clock

Almost every cell in our body has a molecular clock that runs on a ~24-hour cycle. This clock is made of proteins that form a feedback loop — some activate genes, others repress them, creating a daily rhythm. The key repressor proteins are **PER1** and **PER2** (Period 1 and 2).

This clock controls a huge range of biological processes: metabolism, immune function, hormone release, gene expression — essentially the body's daily schedule.

### PerDKO mice

**PerDKO** = Per1/Per2 Double KnockOut. These are mice genetically engineered to lack both PER1 and PER2. Without these repressors, the circadian clock is broken — the feedback loop doesn't work and the mice have no functional circadian rhythm.

By comparing PerDKO mice to wild-type (WT, normal) mice, we can identify genes and processes that are under circadian clock control.

### Why liver?

The liver is one of the most clock-controlled organs in the body. A large fraction of its genes are rhythmically expressed — they go up and down over the 24-hour cycle. It's a particularly good tissue to study circadian biology.

### CT20 — what is circadian time?

CT = Circadian Time. CT0 is the start of the subjective day (even in constant darkness). CT20 is 20 hours into the cycle — roughly the subjective late night. We used mice harvested at CT16-20.

**Constant darkness (DD)** means the mice were kept in the dark without any light cues, so any rhythms we observe are truly endogenous (coming from inside the clock), not responses to light.

### What is alternative splicing?

When a gene is transcribed into RNA, not all parts of the RNA end up in the final messenger RNA (mRNA). Some segments (exons) are included, others (introns) are removed. **Alternative splicing** is when the same gene can produce different mRNA versions by including or skipping different exons. This is how one gene can produce many different protein variants.

**PSI (Percent Spliced In)** is the metric SUPPA uses: for a given exon, what fraction of transcripts include it? PSI = 1 means always included, PSI = 0 means always skipped.

**dPSI** is the difference in PSI between two conditions. A large dPSI means the splicing of that event changes significantly between PerDKO and WT.

---

## The Dataset

We used **GSE130613** from the NCBI Gene Expression Omnibus (GEO), a public repository of sequencing data.

| Property | Value |
|---|---|
| Source | Asher lab, Weizmann Institute |
| Tissue | Mouse liver |
| Conditions | WT vs PerDKO |
| Timepoint | CT16-20, constant darkness |
| Replicates | n=4 per condition |
| Technology | Illumina short-read (MARS-Seq) |
| Accession | GSE130613 / SRP194523 |

**Important caveat:** this dataset used MARS-Seq, a 3' end sequencing protocol. It was originally designed for gene expression quantification, not splicing analysis. We proceeded at Bharath's instruction, and the lower-than-ideal mapping rates (~60% WT, ~46% KO) likely reflect this limitation.

---

## The Pipeline

### Stage 1 — Setup (Bash)

**What we did:** Downloaded the raw sequencing data and quantified transcript expression.

1. **Downloaded 8 FASTQ files** from SRA using `fasterq-dump` (via an adapted version of Bharath's script). FASTQ files contain the raw sequencing reads — each read is a short DNA sequence with a quality score per base.

2. **Built a Salmon index** from the lab's custom mouse transcriptome (`transcriptome_ext.fa`, from Khushi) combined with the mm39 genome FASTA. The index is like a lookup table that lets Salmon map reads ultra-fast.
   - `gentrome.fa` = transcriptome + genome concatenated
   - `decoys.txt` = chromosome names so Salmon knows which sequences are "background"

3. **Ran Salmon quantification** on all 8 samples. For each sample, Salmon maps reads to the transcriptome and estimates how much of each transcript was present. Output: `quant.sf` per sample (columns: transcript ID, TPM, NumReads).

**Key files produced:**
- `raw_data/` — `.sra` files (raw binary format from NCBI)
- `temp/` — `.fastq.gz` files (readable sequencing reads)
- `results/SRP194523_salmon_*/SRR*/quant.sf` — quantification per sample

---

### Stage 2 — Normalisation (R)

**What we did:** Imported Salmon output into R, filtered noise, normalized, and prepared two parallel outputs.

**Script:** `scripts/normalisation_deseq2.R`

#### Track 1 — SUPPA input (Steps 5a–5d)

- **5a. tximport:** reads all 8 `quant.sf` files into R as a transcript-level counts matrix (transcripts × samples)
- **5b. Filter:** kept only transcripts with >5 counts in at least 2 samples — removes noise and very low-expressed transcripts (went from 155,263 → 20,754 transcripts)
- **5c. Normalize:** divided each transcript's counts by the total read count of that sample — makes samples comparable regardless of sequencing depth
- **5d. limma removeBatchEffect:** optional batch correction step — skipped because no batch information was available for this dataset

Output: `results/normalisation/combined_norm.tab` — a normalized expression table (transcripts × 8 samples) used as input to SUPPA

#### Track 2 — DESeq2 differential expression (Step 6)

In parallel, we used the raw counts from tximport (before normalization) to run **DESeq2** — a standard R package for finding differentially expressed genes between two conditions.

DESeq2 models count data statistically and outputs for each gene:
- **log2FoldChange** — how much the gene changes (log2 scale). +2 means 4× more in KO.
- **padj** — adjusted p-value (corrected for multiple testing across thousands of genes)

Output: `results/normalisation/deseq2_KO_vs_WT.tsv`

---

### Stage 3 — Splicing Analysis (SUPPA, Bash)

**What we did:** Ran SUPPA to detect differential splicing events.

**Script:** `scripts/run_suppa.sh`

- **Step 7 — generateEvents:** SUPPA reads the lab's custom GTF (gene annotation file) and generates a catalogue of all possible alternative splicing events. Output: `.ioe` files — one per event type (SE, A3, A5, MX, RI, AF, AL)
- **Step 9 — psiPerEvent:** for each event in the catalogue, SUPPA calculates the PSI value per sample using the normalized expression table. Output: `.psi` files
- **Step 10a — Split by condition:** divided PSI files and expression table into WT and KO groups
- **Step 10b — diffSplice:** compared WT vs KO PSI values to find significantly differentially spliced events. Method: empirical (permutation-based), `-gc` flag corrects p-values per gene

**Significance threshold:** |dPSI| > 0.1 AND p < 0.05

---

### Stage 4 — Results (Python/matplotlib)

**Script:** `scripts/plot_results.py` and `scripts/crossreference.py`

#### Step 11 — dPSI Volcano Plot

Shows all splicing events: x-axis = dPSI (PerDKO minus WT), y-axis = -log10(p-value). Events outside the dashed lines (|dPSI| > 0.1, p < 0.05) are significant. Coloured by event type.

**Result: 35 significant differential splicing events**

#### Step 12 — Cross-reference with Nanopore (Khushi's data)

We compared our 35 Illumina significant events against the significant events from Kiran's Nanopore analysis (provided by Khushi as `.dpsi` files).

Match criterion: exact event ID (both datasets used the same GTF, so IDs are directly comparable).

**Result: 3 events validated by both technologies**
- **Cpeb4** — RNA-binding protein involved in translational regulation, relevant to circadian biology
- **Mup2** — Major Urinary Protein 2, known to be rhythmically expressed in mouse liver
- **0610005C13Rik** — poorly characterized gene

#### Step 13 — DE Volcano Plot

Shows all genes from DESeq2: x-axis = log2FoldChange, y-axis = -log10(padj). Significant genes (|log2FC| > 1.5, padj < 0.05) coloured orange (up in PerDKO) or blue (down in PerDKO).

**Result: 121 significant DE genes (74 up, 47 down in PerDKO)**

#### Bonus — DE genes vs Spliced genes

Cross-reference between the 121 DE genes and the 30 uniquely spliced genes (from the 35 events).

**Result: 0 overlap** — the genes changing expression and the genes changing splicing are completely different sets. This suggests the circadian clock regulates transcription and splicing through independent mechanisms — at least in this dataset.

---

## Summary of Results

| Analysis | Result |
|---|---|
| Salmon mapping rate (WT) | ~60% |
| Salmon mapping rate (KO) | ~46% |
| Significant splicing events (Illumina) | 35 |
| Significant DE genes | 121 |
| Splicing events validated by Nanopore | 3 (Cpeb4, Mup2, 0610005C13Rik) |
| Overlap DE genes vs spliced genes | 0 |

---

## Limitations

- **MARS-Seq is 3' end sequencing** — not ideal for splicing analysis. Reads only cover the last ~75bp of transcripts. Splice junctions are in the middle of transcripts. This explains the lower mapping rates and lower number of splicing events compared to Nanopore.
- **Small sample size (n=4)** — the empirical method in SUPPA has limited statistical power with only 4 samples per condition. P-values can't go very low.
- **No batch correction** — we don't have batch information for this dataset, so limma's removeBatchEffect was not applied.
- **Dataset not designed for splicing** — GSE130613 was a hypoxia study. We used the normoxia controls (21% O2) as a WT vs PerDKO comparison.

---

## Tools Used

| Tool | Purpose |
|---|---|
| `fasterq-dump` (sra-tools) | Download FASTQ from SRA |
| `salmon 2.0.0` | Transcript quantification |
| `gffread` | GTF → transcriptome FASTA |
| `tximport` (R) | Import Salmon output into R |
| `DESeq2` (R) | Differential expression |
| `limma` (R) | Batch correction (optional) |
| `SUPPA 2.3` | Splicing event analysis |
| `matplotlib / pandas` (Python) | Figures |

---

## File Structure

```
orthogonal_validation/
  data/
    transcriptome_ext.fa          ← Khushi's transcriptome FASTA
    transcriptome_ext_index/      ← Salmon index
    gentrome_ext.fa               ← transcriptome + genome
    decoys.txt                    ← chromosome names
    nanopore_suppa/               ← Khushi's Nanopore diffSplice results
  meta/
    SRP194523_Acc_List.txt        ← 8 SRR accession IDs
  results/
    SRP194523_salmon_2026-06-16/  ← quant.sf per sample
    normalisation/
      combined_norm.tab           ← SUPPA expression input
      deseq2_KO_vs_WT.tsv         ← DESeq2 results
      overlap_events.tsv          ← 3 validated events
    suppa_2026-06-17/             ← PSI files, diffSplice output
    figures/                      ← all plots
  scripts/
    get_sra_dump_fastq_run_salmon.sh
    normalisation_deseq2.R
    run_suppa.sh
    plot_results.py
    crossreference.py
```
