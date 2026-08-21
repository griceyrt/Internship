# Orthogonal Validation — Full Pipeline Walkthrough

> **Goal:** Validate splicing events from the lab's Nanopore manuscript using public Illumina RNA-seq data (GSE130613). If the same splicing events appear in both technologies → orthogonal validation confirmed.

---

## 🧬 The Biological Question

| What we're asking | Why it matters |
|---|---|
| Do WT and Per1/2 double-knockout (PerDKO) mouse livers splice RNA differently? | Per1/2 are core circadian clock genes. Losing them disrupts rhythmicity — we want to know if this affects not just *how much* genes are expressed, but *how* their pre-mRNA is processed (spliced). |
| Can we reproduce Nanopore findings with Illumina data? | Nanopore is long-read (reads full transcripts). Illumina is short-read (reads fragments). If both detect the same splicing changes, the result is technology-independent = **orthogonal validation**. |
| What is splicing? | Genes contain exons (coding) and introns (non-coding). After transcription, introns are removed and exons are joined — this is splicing. Alternative splicing means different combinations of exons can be included/excluded, producing different protein isoforms from the same gene. |

---

## 📦 Libraries & Tools Overview

| Tool | Language | What it does in this project |
|---|---|---|
| `sra-tools` (prefetch, fasterq-dump) | Bash | Downloads raw sequencing data from NCBI SRA |
| `cutadapt` | Bash/Python | Removes adapter sequences from reads |
| `salmon` | Bash | Quantifies transcript expression (maps reads to transcriptome) |
| `tximport` | R | Imports Salmon output into R, converts to gene-level counts |
| `DESeq2` | R | Differential expression analysis (which genes change between WT and KO?) |
| `SUPPA2` | Python | Calculates PSI (splicing) per event and detects differential splicing |
| `pandas` | Python | Data wrangling — filtering, merging, crossreferencing tables |
| `matplotlib` | Python | Plotting — volcano plots, Venn diagrams |

---

---

# STAGE 1 — Download, Trim & Quantify

## Step 1 — Download raw sequencing data from NCBI SRA

> **Biology:** The raw data is stored in NCBI's Sequence Read Archive (SRA) as `.sra` files — a compressed format. We need to convert them to FASTQ, the standard format for raw sequencing reads.

| Command | What it does |
|---|---|
| `prefetch -O raw_data/ SRR9002567` | Downloads the `.sra` file for one sample into `raw_data/` |
| `fasterq-dump --threads 6 --outdir temp/ raw_data/SRR9002567/SRR9002567.sra` | Converts the `.sra` to `.fastq` (one file per sample — single-end) |
| `gzip temp/SRR9002567.fastq` | Compresses it to `.fastq.gz` to save disk space |

**What a FASTQ file looks like:**
```
@SRR9002567.1
ACGTACGTACGTACGT...     ← the actual DNA sequence read by the machine
+
IIIIIIIIIIIIIIII...     ← quality scores (how confident the sequencer was per base)
```

> **Why single-end?** MARS-Seq is a 3' end sequencing protocol — it only sequences one end of the fragment (the 3' end of the mRNA). This is cheaper and faster but gives less information than paired-end, which sequences both ends. This is also why our mapping rates are lower than typical RNA-seq.

---

## Step 1b — Adapter trimming (cutadapt)

> **Biology:** During library preparation, synthetic adapter sequences are ligated to RNA fragments so the sequencer can recognize them. If the RNA fragment is short, the sequencer reads past it into the adapter. These adapter bases are not biological — they need to be removed before alignment.

| Command | What it does |
|---|---|
| `cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA` | Scans each read and removes the Illumina TruSeq adapter sequence if found |
| `--minimum-length 20` | Discards reads that become shorter than 20 bases after trimming (too short to align reliably) |
| `-j 7` | Use 7 CPU threads (parallel processing) |
| `-o temp/SRR9002567_trimmed.fastq.gz temp/SRR9002567.fastq.gz` | Input → trimmed output |

> **Result:** SRR9002567 went from **56.71% → 60.55%** mapping rate after trimming. Modest improvement — the low rates are primarily due to MARS-Seq 3' end nature, not adapter contamination.

---

## Step 2 — Build Salmon index (done once)

> **Biology:** Salmon doesn't align reads to the full genome. Instead it maps them to a **transcriptome** — a catalogue of all known transcript sequences. To do this fast, it builds an index (like a book's index) of k-mers (short subsequences of length k=31). We use Khushi's custom transcriptome (`transcriptome_ext.fa`) which includes novel transcripts from the lab's manuscript.

```bash
# Concatenate transcriptome + genome into gentrome (genome acts as decoy)
cat data/transcriptome_ext.fa data/Mus_musculus.GRCm39.dna.primary_assembly.fa \
    > data/gentrome_ext.fa

# Extract chromosome names for decoy list
grep "^>" data/Mus_musculus.GRCm39.dna.primary_assembly.fa \
    | cut -d " " -f 1 | sed 's/>//' > data/decoys.txt

# Build index
salmon index \
    -t data/gentrome_ext.fa \
    -i data/transcriptome_ext_index \
    --decoys data/decoys.txt \
    -k 31 \
    --keepDuplicates
```

| Flag | Meaning |
|---|---|
| `-t gentrome_ext.fa` | Input: transcriptome + genome (genome = decoy to avoid false mappings) |
| `-i transcriptome_ext_index` | Output: the index folder |
| `--decoys decoys.txt` | Tells Salmon: reads mapping only to the genome (not transcriptome) are discarded |
| `-k 31` | K-mer size — 31 is the standard for Illumina short reads |
| `--keepDuplicates` | Keep duplicate transcript sequences (important for custom transcriptomes) |

---

## Step 3 — Salmon quantification

> **Biology:** Salmon reads each FASTQ file and for every read asks: "which transcript does this read come from?" It uses probabilistic mapping (quasi-mapping) — much faster than traditional alignment. The output is a table of estimated read counts and TPM (transcripts per million) per transcript.

```bash
salmon quant \
    -i data/transcriptome_ext_index \
    -l A \
    -p 7 \
    --seqBias --gcBias \
    -r temp/SRR9002567_trimmed.fastq.gz \
    -o results/SRP194523_salmon_2026-06-16/SRR9002567/
```

| Flag | Meaning |
|---|---|
| `-i` | Path to the index we just built |
| `-l A` | Auto-detect library type (Salmon detected **-l U** = unstranded) |
| `-p 7` | 7 threads |
| `--seqBias` | Corrects for sequence-specific bias in read starts |
| `--gcBias` | Corrects for GC content bias |
| `-r` | Single-end reads input (paired-end would use `-1` and `-2`) |
| `-o` | Output folder per sample |

**Output: `quant.sf`** — one per sample, looks like this:

```
Name          Length  EffectiveLength  TPM       NumReads
ENSMUST00001  1200    980.3            12.45     145.2
ENSMUST00002  800     600.1            0.00      0.0
...
```

> **TPM vs NumReads:** TPM (transcripts per million) is normalized for transcript length and sequencing depth — good for comparing between samples. NumReads (estimated counts) is what DESeq2 needs. tximport will handle the conversion.

---

---

# STAGE 2 — Normalization & Differential Expression (R)

## Step 5a — Import with tximport

> **Biology:** We have 8 separate `quant.sf` files. `tximport` reads them all into R and collapses transcript-level estimates to **gene level** using a tx2gene table (transcript ID → gene ID mapping). DESeq2 works at the gene level.

```r
library(tximport)

# Build tx2gene from the GTF
# transcript_id → gene_id mapping
tx2gene <- read.csv("tx2gene.csv")  # built by parsing transcriptome_productivity.gtf

# List all quant.sf files
files <- c(
  "results/.../SRR9002567/quant.sf",
  "results/.../SRR9002568/quant.sf",
  ...  # all 8 samples
)

# Import
txi <- tximport(files,
                type = "salmon",
                tx2gene = tx2gene,
                ignoreTxVersion = TRUE)
```

| Object | What it contains |
|---|---|
| `txi$counts` | Raw estimated counts per gene per sample (for DESeq2) |
| `txi$abundance` | TPM per gene per sample |
| `txi$length` | Effective gene lengths |

---

## Step 5b — Filter

> **Biology:** Many transcripts are detected at near-zero levels — noise, not real expression. We remove them to reduce multiple testing burden and improve statistical power.

```r
# Keep genes with > 5 counts in at least 2 samples
keep <- rowSums(txi$counts > 5) >= 2
counts_filtered <- txi$counts[keep, ]
# 155,263 → 20,754 transcripts retained
```

---

## Step 5c — Normalize (for SUPPA input)

> **Biology:** SUPPA needs normalized expression values to calculate PSI (splicing index). We do a simple library-size normalization: divide each sample's counts by its total read count. This makes samples comparable regardless of sequencing depth.

```r
# Normalize: counts / total counts per sample
norm_counts <- sweep(counts_filtered, 2, colSums(counts_filtered), "/")

# Write combined table for SUPPA
write.table(norm_counts, "results/normalisation/combined_norm.tab",
            sep="\t", quote=FALSE)
```

> **Note:** We do NOT use DESeq2-normalized values for SUPPA. DESeq2 normalization (size factors) is designed for DE, not for PSI calculation. The manuscript uses this simple normalization approach.

---

## Step 6 — DESeq2 differential expression

> **Biology:** DESeq2 models raw counts using a negative binomial distribution (accounts for overdispersion typical in RNA-seq). It estimates size factors (sequencing depth correction) and dispersion (gene-wise variability), then tests for differential expression between WT and PerDKO.

```r
library(DESeq2)

# Create DESeq2 object from tximport
coldata <- data.frame(
  condition = factor(c("WT","WT","WT","WT","KO","KO","KO","KO")),
  row.names = colnames(txi$counts)
)

dds <- DESeqDataSetFromTximport(txi, colData=coldata, design=~condition)

# Run DESeq2
dds <- DESeq(dds)

# Extract results (KO vs WT)
res <- results(dds, contrast=c("condition","KO","WT"))

# Filter significant: |log2FC| > 0.5 AND padj < 0.05
sig <- subset(res, abs(log2FoldChange) > 0.5 & padj < 0.05)
# → 343 significant DE genes
```

| Output column | Meaning |
|---|---|
| `log2FoldChange` | log2(KO / WT) — positive = higher in KO, negative = lower in KO |
| `padj` | p-value adjusted for multiple testing (Benjamini-Hochberg) |
| `baseMean` | Average expression across all samples |

> **Why log2?** A log2FC of 1 means the gene is 2× more expressed in KO. log2FC of -1 means 2× less. Using log2 makes up and down changes symmetric.

---

---

# STAGE 3 — Splicing Analysis (SUPPA2)

## Step 7 — generateEvents

> **Biology:** Before calculating splicing, SUPPA needs to know what alternative splicing events are *possible* given the transcriptome annotation. It reads the GTF file and maps out every possible exon skipping, alternative splice site, intron retention, etc.

```bash
suppa.py generateEvents \
    -i data/transcriptome_productivity.gtf \
    -o results/suppa_2026-06-17/transcriptome_ext \
    -f ioe \
    -e SE SS MX RI FL \
    -b S
```

| Flag | Meaning |
|---|---|
| `-f ioe` | Output format: IOE (isoform-to-event) |
| `-e SE SS MX RI FL` | Event types: **SE**=exon skipping, **SS**=alt splice site, **MX**=mutually exclusive exons, **RI**=intron retention, **FL**=first/last exon |
| `-b S` | Strict boundary matching — confirmed by Khushi and manuscript methods |

**Output:** One `.ioe` file per event type, e.g. `transcriptome_ext_SE.ioe`
Each line = one possible splicing event with the transcripts that include vs exclude it.

---

## Step 9 — psiPerEvent

> **Biology:** PSI = **Percent Spliced In**. For a given splicing event, PSI measures the fraction of transcripts that *include* the alternative sequence. PSI = 1.0 means 100% of transcripts include it; PSI = 0.0 means none do. SUPPA calculates PSI per sample using the normalized expression table.

```bash
suppa.py psiPerEvent \
    -i results/suppa_2026-06-17/transcriptome_ext_SE.ioe \
    -e results/normalisation/combined_norm.tab \
    -o results/suppa_2026-06-17/
```

**Output:** `.psi` files — a matrix of PSI values (events × samples)

```
event_id                    SRR9002567  SRR9002568  ...
Cpeb4;SE:...                0.82        0.79        ...
Mup2;SE:...                 0.45        0.41        ...
```

---

## Step 10a — Split by condition

```bash
# Split PSI files: WT columns vs KO columns
cut -f1,2,3,4,5 SE.psi > SE_WT.psi      # samples 1-4 = WT
cut -f1,6,7,8,9 SE.psi > SE_KO.psi      # samples 5-8 = KO

# Same for the expression table
cut -f1,2,3,4,5 combined_norm.tab > WT_norm.tab
cut -f1,6,7,8,9 combined_norm.tab > KO_norm.tab
```

---

## Step 10b — diffSplice

> **Biology:** diffSplice compares PSI values between conditions. For each event, it computes dPSI = PSI(WT) - PSI(KO) and a p-value. A significant event means the splicing pattern changed between WT and PerDKO mice — the clock loss affects how that gene is spliced.

```bash
suppa.py diffSplice \
    --method empirical \
    --input results/suppa_2026-06-17/transcriptome_ext_SE.ioe \
    --psi results/suppa_2026-06-17/SE_WT.psi \
           results/suppa_2026-06-17/SE_KO.psi \
    --tpm WT_norm.tab KO_norm.tab \
    -gc \
    -o results/suppa_2026-06-17/diff/diff_SE
```

| Flag | Meaning |
|---|---|
| `--method empirical` | Uses bootstrapping to estimate p-values (better for small n) |
| `-gc` | Generates cluster files for the empirical method |

**Output:** `diff_SE.dpsi.temp.0`

```
event_id           dPSI     p-value
Cpeb4;SE:...       -0.23    0.001
Mup2;SE:...         0.18    0.031
```

> **Sign convention:** dPSI in the file = WT − KO. We negate it in the plots to show KO − WT (positive = more inclusion in KO).

**Significance thresholds:**
- Illumina (this project): `|dPSI| > 0.1` AND `p < 0.3` ← relaxed because n=4 with empirical method has low power
- Nanopore (Khushi's data): `|dPSI| > 0.1` AND `p < 0.05` ← stricter, more samples

---

---

# STAGE 4 — Results (Python)

## Step 11 — dPSI Volcano plot

> **Biology:** A volcano plot visualizes all splicing events at once. Events far right or left = big splicing change. Events high up = statistically significant. Events in both = the ones we care about.

```python
import pandas as pd
import matplotlib.pyplot as plt

# Load all event types and concatenate
dfs = []
for event_type in ["SE", "A3", "A5", "MX", "RI", "AF", "AL"]:
    df = pd.read_csv(f"diff_{event_type}.dpsi.temp.0", sep="\t")
    df["event_type"] = event_type
    dfs.append(df)
all_events = pd.concat(dfs)

# Negate dPSI: file has WT-KO, we want KO-WT
all_events["dPSI"] = -all_events["dPSI"]

# Significance filter
sig = all_events[
    (all_events["dPSI"].abs() > 0.1) &
    (all_events["p-value"] < 0.3)
]
# → 532 significant events, 457 unique genes

# Plot
plt.scatter(all_events["dPSI"], -np.log10(all_events["p-value"]), ...)
plt.axvline(x=0.1); plt.axvline(x=-0.1)   # dPSI threshold
plt.axhline(y=-np.log10(0.3))              # p threshold
```

---

## Step 12 — Crossreference with Nanopore

> **Biology:** Khushi's Nanopore data went through the same SUPPA pipeline but on long reads. We match events by their exact SUPPA event ID (e.g. `Cpeb4;SE:chr15:52100000-52101000:+`). An exact match means both technologies detected the same splicing change in the same location.

```python
# Load Illumina significant events
illumina_sig = all_events[
    (all_events["dPSI"].abs() > 0.1) &
    (all_events["p-value"] < 0.3)
]

# Load Nanopore events (Khushi's files)
nano_dfs = []
for f in glob("data/nanopore_suppa/res_*.dpsi.temp.0"):
    nano_dfs.append(pd.read_csv(f, sep="\t"))
nanopore_all = pd.concat(nano_dfs)

# Filter Nanopore at stricter threshold
nanopore_sig = nanopore_all[
    (nanopore_all["dPSI"].abs() > 0.1) &
    (nanopore_all["p-value"] < 0.05)
]
# → 971 Nanopore events

# Find overlap by event ID
overlap = illumina_sig[illumina_sig.index.isin(nanopore_sig.index)]
# → 32 overlapping events across 31 unique genes
```

**Validated at strict p<0.05 in both:** `Cpeb4`, `Mup2`, `0610005C13Rik`

---

## Step 13 — DE Volcano

```python
de = pd.read_csv("results/normalisation/deseq2_KO_vs_WT.tsv", sep="\t")

sig_de = de[
    (de["log2FoldChange"].abs() > 0.5) &
    (de["padj"] < 0.05)
]
# → 343 DE genes

plt.scatter(de["log2FoldChange"], -np.log10(de["padj"]), ...)
plt.axvline(x=0.5); plt.axvline(x=-0.5)
plt.axhline(y=-np.log10(0.05))

# Label 27 genes selected by Kiran (circadian, Cyp450, immune)
for gene in kiran_genes:
    plt.annotate(gene, ...)
```

---

## Step 14 — DE vs Splicing Venn

> **Biology:** A gene can be affected in two independent ways: its total expression level can change (DE) and/or how it's spliced can change. These are not the same thing — a gene can be spliced differently without changing total expression. The Venn shows how much overlap there is between the two phenomena.

```python
de_genes    = set(sig_de["gene_name"])           # 343 genes
splice_genes = set(illumina_sig["gene_name"])     # 457 genes

only_de      = de_genes - splice_genes            # 325
both         = de_genes & splice_genes            # 18
only_splice  = splice_genes - de_genes            # 439
```

---

---

# 🔑 Key Concepts Glossary

| Term | Definition |
|---|---|
| **PSI** | Percent Spliced In — fraction of transcripts including a given exon/sequence. Range 0–1. |
| **dPSI** | Delta PSI — difference in PSI between two conditions. We show KO − WT. |
| **TPM** | Transcripts Per Million — normalized expression unit accounting for transcript length and sequencing depth. |
| **log2FC** | log2 fold change — log2(KO/WT). Positive = higher in KO. |
| **padj** | Adjusted p-value — corrected for multiple testing across thousands of genes. |
| **quasi-mapping** | Salmon's fast mapping method — doesn't do full alignment, uses k-mer matching. |
| **tx2gene** | A table linking transcript IDs to gene IDs — needed to collapse transcript counts to gene level. |
| **IOE file** | Isoform-to-event — SUPPA's format listing which transcripts include/exclude each splicing event. |
| **empirical method** | SUPPA's bootstrapping approach for p-values — better suited for small sample sizes than the asymptotic method. |
| **MARS-Seq** | Massively parallel RNA single-cell sequencing — 3' end protocol, single-end, lower mapping rates than standard RNA-seq. |
| **Orthogonal validation** | Confirming a finding using a completely different technology or dataset — strengthens confidence in the result. |
