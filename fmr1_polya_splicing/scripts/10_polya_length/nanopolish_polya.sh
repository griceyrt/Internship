#!/bin/bash
#SBATCH --job-name=nanopolish_polya
#SBATCH --partition=Lake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=%j_nanopolish_polya.log

# Phase 6 (poly(A) length): nanopolish index + nanopolish polya per sample.
# Uses fast5_pass reads only (higher-confidence basecalls) and the EXISTING
# Phase 2 transcriptome-aligned, sorted+indexed BAMs -- no realignment
# needed, Nanopolish polya works directly off what quantification already
# produced.
#
# Sample name mapping (fast5 folders from Japan use a different naming
# convention than this project's own WT_Rep1/KO_Rep1 style -- verified
# by nanopolish index itself: a genuine mismatch would show as very few
# or zero reads successfully indexed for that pair, watch the per-sample
# summary below for that signal):
#   WT_Rep1 -> 210126_DrLies_WT1_CTX   KO_Rep1 -> 210126_DrLies_KO1_CTX
#   WT_Rep2 -> 210126_DrLies_WT2_CTX   KO_Rep2 -> 210126_DrLies_KO2_CTX
#   WT_Rep3 -> 210126_DrLies_WT3_CTX   KO_Rep3 -> 210126_DrLies_KO3_CTX

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biotools

# nanopolish has not been used in this project before -- if the command
# below fails with "command not found", install it into biotools first:
#   conda install -c bioconda -y nanopolish
if ! command -v nanopolish &> /dev/null; then
    echo "ERROR: nanopolish not found in the biotools env. Install with:"
    echo "  conda install -c bioconda -y nanopolish"
    exit 1
fi

BASE_DIR="$HOME/fmr1_polya_splicing"
FAST5_BASE="/scratch/Lake/grangel/fmr1_polya_splicing/fast5"
ALIGN_BASE="$BASE_DIR/results/gse188840_minimap2_2026-08-03"
TRANSCRIPTOME="$BASE_DIR/data/reference/Mus_musculus.GRCm39.cdna.all.fa"
TEMP_DIR="$BASE_DIR/temp"
OUTPUT_DIR="$BASE_DIR/results/polya_length"
CORES="$SLURM_CPUS_PER_TASK"

mkdir -p "$OUTPUT_DIR"

if [ ! -f "${TRANSCRIPTOME}.fai" ]; then
    echo "ERROR: ${TRANSCRIPTOME}.fai not found -- run prep_transcriptome_faidx.sh first."
    exit 1
fi

declare -A SAMPLE_TO_FAST5=(
    ["WT_Rep1"]="210126_DrLies_WT1_CTX"
    ["WT_Rep2"]="210126_DrLies_WT2_CTX"
    ["WT_Rep3"]="210126_DrLies_WT3_CTX"
    ["KO_Rep1"]="210126_DrLies_KO1_CTX"
    ["KO_Rep2"]="210126_DrLies_KO2_CTX"
    ["KO_Rep3"]="210126_DrLies_KO3_CTX"
)

for SAMPLE in "${!SAMPLE_TO_FAST5[@]}"; do
    FAST5_SAMPLE="${SAMPLE_TO_FAST5[$SAMPLE]}"
    echo ""
    echo "======================================"
    echo "Processing $SAMPLE (fast5 dir: $FAST5_SAMPLE)"
    echo "======================================"

    FASTQ="$TEMP_DIR/${SAMPLE}.fastq.gz"
    FAST5_DIR="$FAST5_BASE/$FAST5_SAMPLE/fast5_pass"
    BAM="$ALIGN_BASE/$SAMPLE/${SAMPLE}.sorted.bam"
    OUT_TSV="$OUTPUT_DIR/${SAMPLE}_polya.tsv"

    if [ ! -f "$FASTQ" ]; then
        echo "ERROR: fastq not found at $FASTQ -- skipping $SAMPLE."
        continue
    fi
    if [ ! -d "$FAST5_DIR" ]; then
        echo "ERROR: fast5_pass dir not found at $FAST5_DIR -- skipping $SAMPLE."
        continue
    fi
    if [ ! -f "$BAM" ]; then
        echo "ERROR: sorted BAM not found at $BAM -- skipping $SAMPLE."
        continue
    fi

    if [ ! -f "${FASTQ}.index" ]; then
        echo "[1/2] nanopolish index..."
        nanopolish index -d "$FAST5_DIR" "$FASTQ"
    else
        echo "[1/2] Index already exists, skipping."
    fi

    if [ -f "$OUT_TSV" ]; then
        echo "[2/2] $OUT_TSV already exists, skipping."
    else
        echo "[2/2] nanopolish polya..."
        nanopolish polya --threads="$CORES" \
            --reads="$FASTQ" \
            --bam="$BAM" \
            --genome="$TRANSCRIPTOME" \
            > "$OUT_TSV"

        N_READS=$(tail -n +2 "$OUT_TSV" | wc -l)
        N_PASS=$(tail -n +2 "$OUT_TSV" | awk -F'\t' '$NF=="PASS"' | wc -l)
        echo "$SAMPLE: $N_READS reads processed, $N_PASS with qc_tag=PASS."
        if [ "$N_READS" -lt 1000 ]; then
            echo "WARNING: very few reads for $SAMPLE -- check the fast5/fastq sample mapping above,"
            echo "this could mean $FAST5_SAMPLE does not actually correspond to $SAMPLE."
        fi
    fi
done

echo ""
echo "======================================"
echo "All samples processed. Per-read poly(A) tables in: $OUTPUT_DIR"
echo "======================================"
