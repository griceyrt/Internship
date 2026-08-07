#!/bin/bash
#SBATCH --job-name=igv_check_pik3cd
#SBATCH --partition=Lake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=%j_igv_check_pik3cd.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biotools

# =============================================================================
# IGV sanity check -- Pik3cd, requested by Kiran 2026-08-06
# Author: Gricey
#
# Kiran wants a visual (IGV) confirmation that Pik3cd genuinely has the many
# isoforms shown in the stacked-bar figure (n=14 in WT, n=16 in KO at
# min50TPM). The existing pipeline's BAM (minimap2_align_sort_index.sh) is
# aligned to the TRANSCRIPTOME for oarfish/Salmon quantification -- each
# transcript is its own flat linear reference there, so it can't show
# intron-spanning splice structure the way a genome browser view can.
# Classic IGV isoform validation needs reads aligned to the GENOME with a
# splice-aware preset, viewed against the GTF gene model (ideally as a
# Sashimi plot, which shows per-junction read support directly).
#
# Rather than re-align a full sample's ~4-6M reads to the genome (expensive
# and unnecessary for a single-gene check), this script:
#   1) pulls Pik3cd's transcript IDs from the GTF
#   2) extracts only the reads that mapped to those transcripts from the
#      EXISTING transcriptome BAM (fast -- tiny subset of reads)
#   3) re-aligns just that subset to the genome with minimap2's direct-RNA
#      spliced preset (-ax splice -uf -k14, the standard ONT dRNA-to-genome
#      settings; -uf assumes reads are all sense-strand, true for direct RNA)
#   4) sorts + indexes the result -- small enough to download and view in
#      IGV directly on your Mac
#
# Runs for ONE representative sample per condition by default (WT_Rep1,
# KO_Rep1) -- enough for a qualitative "do we see many isoforms" check.
# Add more names to SAMPLES below if you want more replicates.
# =============================================================================

BASE_DIR="$HOME/fmr1_polya_splicing"
GTF="$BASE_DIR/data/reference/Mus_musculus.GRCm39.116.gtf.gz"
GENOME="$BASE_DIR/data/reference/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
MINIMAP2_DIR=$(find "$BASE_DIR/results" -maxdepth 1 -name "gse188840_minimap2_*" | sort | tail -1)
OUT_DIR="$BASE_DIR/results/igv_check_pik3cd"
GENE_ID="ENSMUSG00000039936"   # Pik3cd
CORES="$SLURM_CPUS_PER_TASK"
SAMPLES=("WT_Rep1" "KO_Rep1")

mkdir -p "$OUT_DIR"

echo "Using Minimap2 BAM dir: $MINIMAP2_DIR"
if [ -z "$MINIMAP2_DIR" ]; then
    echo "ERROR: could not find results/gse188840_minimap2_* -- check the path."
    exit 1
fi

echo "[1/3] Getting Pik3cd transcript IDs from the GTF..."
TX_IDS=$(zcat "$GTF" | awk -F'\t' -v gid="$GENE_ID" '$3=="transcript" && $9 ~ gid {print $9}' \
    | grep -oP 'transcript_id "\K[^"]+' | sort -u)
N_TX=$(echo "$TX_IDS" | grep -c .)
echo "Transcript IDs found ($N_TX total):"
echo "$TX_IDS"

for SAMPLE in "${SAMPLES[@]}"; do
    echo ""
    echo "=== $SAMPLE ==="
    BAM="$MINIMAP2_DIR/$SAMPLE/$SAMPLE.sorted.bam"
    if [ ! -f "$BAM" ]; then
        echo "ERROR: $BAM not found -- skipping $SAMPLE."
        continue
    fi

    # BAM reference names carry Ensembl version suffixes (from the cDNA
    # fasta headers, e.g. ENSMUST00000037596.9) but the GTF's transcript_id
    # doesn't -- same version mismatch as the SUPPA/tximport issue, resolved
    # the same way: look up the exact versioned name from the BAM header.
    REGIONS=()
    for TX in $TX_IDS; do
        FULL=$(samtools view -H "$BAM" | grep -oP "SN:${TX}\.[0-9]+" | sed 's/SN://' | head -1)
        [ -n "$FULL" ] && REGIONS+=("$FULL")
    done
    echo "Matched ${#REGIONS[@]} of $N_TX transcript IDs to references in $SAMPLE's BAM."

    if [ ${#REGIONS[@]} -eq 0 ]; then
        echo "ERROR: no matching transcript references found for $SAMPLE -- skipping."
        continue
    fi

    echo "[2/3] Extracting Pik3cd reads + re-aligning to the genome (splice-aware)..."
    samtools view -b "$BAM" "${REGIONS[@]}" \
        | samtools fastq - 2>/dev/null \
        | minimap2 -ax splice -uf -k14 -t "$CORES" "$GENOME" - 2>"$OUT_DIR/${SAMPLE}_minimap2.log" \
        | samtools sort -@ "$CORES" -o "$OUT_DIR/${SAMPLE}_pik3cd_genome.sorted.bam" -

    samtools index "$OUT_DIR/${SAMPLE}_pik3cd_genome.sorted.bam"
    N_READS=$(samtools view -c "$OUT_DIR/${SAMPLE}_pik3cd_genome.sorted.bam")
    echo "$SAMPLE: $N_READS reads aligned to the Pik3cd genomic locus."
done

echo ""
echo "[3/3] Done. These are small -- download to your Mac for IGV:"
ls -lh "$OUT_DIR"/*.bam "$OUT_DIR"/*.bai 2>/dev/null
echo ""
echo "Pik3cd locus (GRCm39/mm39): chr4:149,732,411-149,787,028 (- strand)"
echo "(Note: this genome's contigs are named Ensembl-style, e.g. '4', not 'chr4' --"
echo " use whichever your IGV genome load ends up showing in the chromosome list.)"
