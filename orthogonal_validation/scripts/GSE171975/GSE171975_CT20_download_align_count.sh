#!/bin/bash
#SBATCH --job-name=gse171975_ct20
#SBATCH --partition=Lake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=%j_gse171975_ct20.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biotools

# GSE171975 (Aviram/Manella/Asher/Golik, PLOS Biology 2021), CT20 window only,
# second independent Asher dataset (WT n=4, PerDKO n=3). Full chain: SRA
# download -> cutadapt trim -> STAR align (extended-GTF index) -> htseq-count
# (against the original non-extended GTF, matching the adopted GSE130613
# approach). Author: Gricey
#
# Assumes single-end reads (not independently confirmed from the GEO record --
# if fasterq-dump produces _1/_2 files instead, this is paired-end and the
# script needs adjusting).
#
# MultiQC in this env needs conda-forge kaleido (not pip) or plots render
# blank -- see reference_biotools_env_gotchas.

BASE_DIR="$HOME/orthogonal_validation"
DATA_DIR="$BASE_DIR/data"
INDEX_DIR="$DATA_DIR/star_index_extended"
GTF_ORIGINAL="$DATA_DIR/transcriptome_productivity.gtf"
OUTDIR="$BASE_DIR/results/GSE171975/GSE171975_CT20_$(date +"%Y-%m-%d")"
TEMP_DIR="$BASE_DIR/temp/GSE171975"
FASTQC_RAW_DIR="$BASE_DIR/results/GSE171975/GSE171975_CT20_fastqc_raw"
FASTQC_TRIMMED_DIR="$BASE_DIR/results/GSE171975/GSE171975_CT20_fastqc_trimmed"
CORES="$SLURM_CPUS_PER_TASK"

mkdir -p "$OUTDIR" "$TEMP_DIR" "$FASTQC_RAW_DIR" "$FASTQC_TRIMMED_DIR"

# sample_label:SRR_accession -- SRR (run), not SRX (experiment): prefetch
# resolves and downloads under the SRR name
SAMPLES=(
    "WT_A:SRR14227295"
    "WT_B:SRR14227296"
    "WT_C:SRR14227297"
    "WT_D:SRR14227298"
    "PerDKO_A:SRR14227343"
    "PerDKO_B:SRR14227344"
    "PerDKO_C:SRR14227345"
)

for entry in "${SAMPLES[@]}"; do
    LABEL="${entry%%:*}"
    ACC="${entry##*:}"

    echo ""
    echo "======================================"
    echo "Processing $LABEL ($ACC)"
    echo "======================================"

    echo "[1/6] Downloading ($ACC)..."
    prefetch "$ACC" -O "$TEMP_DIR"
    fasterq-dump "$TEMP_DIR/$ACC" -O "$TEMP_DIR" --threads "$CORES"

    RAW_FASTQ="$TEMP_DIR/${ACC}.fastq"
    if [ ! -f "$RAW_FASTQ" ]; then
        echo "ERROR: expected single-end file $RAW_FASTQ not found."
        echo "       If ${ACC}_1.fastq / ${ACC}_2.fastq exist instead, this"
        echo "       sample is paired-end -- script needs adjusting, skipping."
        continue
    fi
    gzip -f "$RAW_FASTQ"

    echo "[2/6] FastQC (raw)..."
    fastqc "${RAW_FASTQ}.gz" -o "$FASTQC_RAW_DIR" --quiet

    echo "[3/6] Trimming (cutadapt)..."
    TRIMMED="$TEMP_DIR/${LABEL}_trimmed_full.fastq.gz"
    cutadapt \
        -a AGATCGGAAGAG -a "A{12}" -a "T{12}" -a "G{12}" \
        --minimum-length 20 \
        -o "$TRIMMED" \
        "${RAW_FASTQ}.gz" \
        > "$TEMP_DIR/${LABEL}_cutadapt.log"

    echo "[4/6] FastQC (trimmed)..."
    fastqc "$TRIMMED" -o "$FASTQC_TRIMMED_DIR" --quiet

    SAMPLE_OUT="$OUTDIR/$LABEL"
    mkdir -p "$SAMPLE_OUT"

    echo "[5/6] STAR alignment..."
    STAR --runMode alignReads \
        --genomeDir "$INDEX_DIR" \
        --readFilesIn "$TRIMMED" \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix "${SAMPLE_OUT}/${LABEL}_" \
        --runThreadN "$CORES"

    BAM="${SAMPLE_OUT}/${LABEL}_Aligned.sortedByCoord.out.bam"

    echo "[6/6] htseq-count (original GTF)..."
    htseq-count \
        -f bam \
        -r pos \
        -s no \
        -m union \
        --additional-attr=gene_name \
        "$BAM" \
        "$GTF_ORIGINAL" \
        > "${SAMPLE_OUT}/${LABEL}_counts.txt"

    echo "Done with $LABEL."
done

echo ""
echo "Running MultiQC (raw)..."
multiqc "$FASTQC_RAW_DIR" -o "$FASTQC_RAW_DIR" --quiet

echo "Running MultiQC (trimmed)..."
multiqc "$FASTQC_TRIMMED_DIR" -o "$FASTQC_TRIMMED_DIR" --quiet

echo ""
echo "======================================"
echo "All samples processed (or check log above for per-sample errors)."
echo "Results in: $OUTDIR"
echo "Raw QC:     $FASTQC_RAW_DIR/multiqc_report.html"
echo "Trimmed QC: $FASTQC_TRIMMED_DIR/multiqc_report.html"
echo "======================================"
