#!/bin/bash
#SBATCH --job-name=gse188840_suppa
#SBATCH --partition=Lake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=%j_gse188840_suppa.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biotools

# =============================================================================
# SUPPA pipeline for GSE188840
# Author: Gricey
# Adapted from: boundary_analysis/scripts/run_suppa.sh (project 1)
#
# Differences from the project 1 template:
#   - 2 conditions (WT, KO), 3 reps each -- not 3 conditions/9 samples
#   - TPM matrix comes directly from normalisation_deseq2_oarfish.R's
#     tximport export (TPM_all_samples.tpm) -- no need to rebuild it from
#     quant.sf, R already produced the combined, non-versioned-ID matrix
#   - GTF is the stock mm39 Ensembl GTF (gzipped -- decompressed below),
#     not project 1's custom transcriptome_productivity.gtf
#   - Boundary mode fixed at V/1nt (Kiran's "+1 variable position" ask,
#     already settled), not parameterized like the original script
#   - No clustering step -- not needed for Overlap #1
#   - SUPPA itself is cloned fresh here (self-contained), same
#     clone-if-missing pattern used in the project 4 (HFD) script, rather
#     than depending on wherever project 4 happened to put its own copy
# =============================================================================

BASE_DIR="$HOME/fmr1_polya_splicing"
CODE_DIR="$BASE_DIR/code"
GTF_GZ="$BASE_DIR/data/reference/Mus_musculus.GRCm39.116.gtf.gz"
GTF="$BASE_DIR/data/reference/Mus_musculus.GRCm39.116.gtf"
TPM_DIR="$BASE_DIR/results/normalisation"
OUTPUT_DIR="$BASE_DIR/results/gse188840_suppa_$(date +"%Y-%m-%d")"
EVENTS_DIR="${OUTPUT_DIR}/events"
PSI_DIR="${OUTPUT_DIR}/psi"
DIFF_DIR="${OUTPUT_DIR}/diff"

mkdir -p "$CODE_DIR" "$EVENTS_DIR" "$PSI_DIR" "$DIFF_DIR"

echo "=== Setting up SUPPA2 (clone once, reused by future scripts) ==="
if [ ! -f "$CODE_DIR/SUPPA/suppa.py" ]; then
    git clone https://github.com/comprna/SUPPA.git "$CODE_DIR/SUPPA"
else
    echo "SUPPA already present at $CODE_DIR/SUPPA, skipping clone."
fi
SUPPA="python3 $CODE_DIR/SUPPA/suppa.py"

echo "=== Decompressing GTF (SUPPA needs a plain .gtf, not .gz) ==="
if [ ! -f "$GTF" ]; then
    gunzip -k "$GTF_GZ"
else
    echo "$GTF already exists, skipping decompression."
fi

echo ""
echo "=== STEP 1: Splitting TPM matrix by condition ==="
awk 'BEGIN{OFS="\t"} NR==1{print "WT_Rep1","WT_Rep2","WT_Rep3"} NR>1{print $1,$2,$3,$4}' \
    "$TPM_DIR/TPM_all_samples.tpm" > "$TPM_DIR/TPM_WT.tpm"
awk 'BEGIN{OFS="\t"} NR==1{print "KO_Rep1","KO_Rep2","KO_Rep3"} NR>1{print $1,$5,$6,$7}' \
    "$TPM_DIR/TPM_all_samples.tpm" > "$TPM_DIR/TPM_KO.tpm"

echo ""
echo "=== STEP 2: Generating events (boundary=V, threshold=1nt) ==="
$SUPPA generateEvents \
    -i "$GTF" \
    -o "${EVENTS_DIR}/events" \
    -f ioe \
    -e SE SS MX RI FL \
    -b V \
    -t 1

echo ""
echo "=== STEP 3: PSI per event ==="
for EVENT in SE A3 A5 MX RI AF AL; do
    echo "  Processing ${EVENT}..."
    $SUPPA psiPerEvent \
        --ioe-file "${EVENTS_DIR}/events_${EVENT}_variable_1.ioe" \
        --expression-file "$TPM_DIR/TPM_all_samples.tpm" \
        -o "${PSI_DIR}/${EVENT}"
done

echo ""
echo "=== STEP 4: Splitting PSI files by condition ==="
for EVENT in SE A3 A5 MX RI AF AL; do
    awk 'BEGIN{OFS="\t"} NR==1{print "WT_Rep1","WT_Rep2","WT_Rep3"} NR>1{print $1,$2,$3,$4}' \
        "${PSI_DIR}/${EVENT}.psi" > "${PSI_DIR}/${EVENT}_WT.psi"
    awk 'BEGIN{OFS="\t"} NR==1{print "KO_Rep1","KO_Rep2","KO_Rep3"} NR>1{print $1,$5,$6,$7}' \
        "${PSI_DIR}/${EVENT}.psi" > "${PSI_DIR}/${EVENT}_KO.psi"
done

echo ""
echo "=== STEP 5: diffSplice (WT vs KO, method=empirical) ==="
cd "${OUTPUT_DIR}"
for EVENT in SE A3 A5 MX RI AF AL; do
    echo "  Processing ${EVENT}..."
    $SUPPA diffSplice \
        --method empirical \
        --input "events/events_${EVENT}_variable_1.ioe" \
        --psi "psi/${EVENT}_WT.psi" "psi/${EVENT}_KO.psi" \
        --tpm "$TPM_DIR/TPM_WT.tpm" "$TPM_DIR/TPM_KO.tpm" \
        -gc \
        -o "diff/diff_${EVENT}"
done

echo ""
echo "=== SUMMARY (p < 0.05) ==="
echo "NOTE: with only 2 conditions (vs project 1's 3), the diffSplice output"
echo "filename may be diff_EVENT.dpsi rather than .dpsi.temp.0 -- check"
echo "'ls diff/' if this summary loop comes back empty/wrong and adjust."
for EVENT in SE A3 A5 MX RI AF AL; do
    COUNT=$(awk 'NR>1 && $3!="nan" && $3+0 < 0.05' "diff/diff_${EVENT}.dpsi.temp.0" 2>/dev/null | wc -l)
    echo "  ${EVENT}: ${COUNT} events"
done

echo ""
echo "======================================"
echo "SUPPA pipeline complete."
echo "Results in: $OUTPUT_DIR"
echo "======================================"
