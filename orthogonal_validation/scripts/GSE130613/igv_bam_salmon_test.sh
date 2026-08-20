#!/bin/bash
# Generates a Salmon-based BAM for ONE sample for IGV, using Salmon's
# --writeMappings flag. Author: Gricey
#
# NOT a STAR genome alignment -- Salmon maps to transcriptome_ext.fa
# (transcripts + genome decoys), so the BAM has TRANSCRIPT coordinates, not
# chromosome coordinates. Load transcriptome_ext.fa itself as the IGV
# reference genome, NOT mm39.
#
# USAGE: bash scripts/igv_bam_salmon_test.sh <SRR_id>
# Run from: orthogonal_validation/   (needs salmon + samtools on PATH,
#           e.g. `conda activate suppa_env`)
#
# Example: bash scripts/igv_bam_salmon_test.sh SRR9002585

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${SCRIPT_DIR}/.."

SRR="$1"
if [ -z "$SRR" ]; then
    echo "Usage: bash scripts/igv_bam_salmon_test.sh <SRR_id>"
    echo "Example: bash scripts/igv_bam_salmon_test.sh SRR9002585"
    exit 1
fi

SALMON_INDEX="${BASE_DIR}/data/transcriptome_ext_index"
TRIMMED="${BASE_DIR}/temp/${SRR}_trimmed_full.fastq.gz"
OUTDIR="${BASE_DIR}/results/GSE130613/igv_test_${SRR}_$(date +"%Y-%m-%d")"
SAM="${OUTDIR}/${SRR}.sam"
BAM="${OUTDIR}/${SRR}.bam"
SORTED_BAM="${OUTDIR}/${SRR}.sorted.bam"

if [ ! -f "$TRIMMED" ]; then
    echo "Trimmed FASTQ not found: $TRIMMED"
    echo "This sample hasn't been downloaded/trimmed yet. Run"
    echo "  bash scripts/get_sra_dump_fastq_run_salmon.sh SRP194523"
    echo "first, or pick a different SRR that already has a *_trimmed_full.fastq.gz in temp/."
    exit 1
fi

mkdir -p "$OUTDIR"

echo "[1/4] Running Salmon quant with --writeMappings for $SRR..."
salmon quant -i "$SALMON_INDEX" \
    -l A \
    -p 4 \
    --seqBias --gcBias \
    -r "$TRIMMED" \
    --writeMappings="$SAM" \
    -o "$OUTDIR" \
    2>&1 | tee "${OUTDIR}/${SRR}_salmon_igv.log"

echo "[2/4] Converting SAM -> BAM..."
samtools view -bS "$SAM" > "$BAM"

echo "[3/4] Sorting BAM..."
samtools sort "$BAM" -o "$SORTED_BAM"

echo "[4/4] Indexing sorted BAM..."
samtools index "$SORTED_BAM"

# Drop the intermediate SAM/unsorted BAM, keep only the sorted+indexed one
rm -f "$SAM" "$BAM"

echo ""
echo "======================================"
echo "Done. To view in IGV:"
echo "  1. Genomes > Load Genome from File...  ->  ${BASE_DIR}/data/transcriptome_ext.fa"
echo "  2. File > Load from File...             ->  ${SORTED_BAM}"
echo "  3. In the location box, jump to transcript: c86ce1ed-b815-4b2a-8559-9a2751272523 (novel Saa1 isoform)"
echo "======================================"
