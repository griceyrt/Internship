#!/bin/bash
# =============================================================================
# SUPPA Pipeline Script
# Author: Gricey
# Description: Runs the full SUPPA pipeline from generateEvents to clusterEvents
#              All parameters are defined at the top for easy modification
# =============================================================================

# =============================================================================
# PARAMETERS — modify these to run different analyses
# =============================================================================

# Boundary mode: S = Strict, V = Variable
BOUNDARY="V"

# Threshold (only used when BOUNDARY=V, ignored when BOUNDARY=S)
THRESHOLD=1

# Minimum total expression filter for psiPerEvent (0 = no filter)
TPM_FILTER=0

# diffSplice method: empirical or classical
DIFF_METHOD="empirical"

# clusterEvents parameters
CLUSTER_EPS=0.5
CLUSTER_MINPTS=5
CLUSTER_SIG=0.05
CLUSTER_DPSI=0.1
CLUSTER_METHOD="OPTICS"

# =============================================================================
# PATHS — modify these if your folder structure changes
# =============================================================================

# Base directory (where salmon/ and GTF file are)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${SCRIPT_DIR}/../data"

# GTF file
GTF="${BASE_DIR}/transcriptome_productivity.gtf"

# Salmon TPM files (individual per sample)
SALMON_DIR="${BASE_DIR}/salmon"

# SUPPA script
SUPPA="python3 ${SCRIPT_DIR}/../SUPPA/suppa.py"

# =============================================================================
# OUTPUT FOLDER — named automatically based on parameters
# =============================================================================

if [ "$BOUNDARY" == "S" ]; then
    RUN_NAME="strict"
else
    RUN_NAME="variable_${THRESHOLD}nt"
fi

OUTPUT_DIR="${BASE_DIR}/output_${RUN_NAME}"
EVENTS_DIR="${OUTPUT_DIR}/events"
TPM_DIR="${OUTPUT_DIR}/tpm"
PSI_DIR="${OUTPUT_DIR}/psi"
DIFF_DIR="${OUTPUT_DIR}/diff"
CLUSTER_DIR="${OUTPUT_DIR}/cluster"

# =============================================================================
# CREATE OUTPUT FOLDERS
# =============================================================================

echo "Creating output folders..."
mkdir -p "${EVENTS_DIR}" "${TPM_DIR}" "${PSI_DIR}" "${DIFF_DIR}" "${CLUSTER_DIR}"

# =============================================================================
# STEP 1: PREPARE TPM FILES
# =============================================================================

echo ""
echo "=== STEP 1: Preparing TPM files ==="

# Extract transcript ID and TPM column from each sample
for i in 1 2 3 4 5 6 7 8 9; do
    echo "  Processing KM${i}..."
    echo -e "KM${i}" > "${TPM_DIR}/KM${i}.tpm"
    awk 'NR>1{print $1"\t"$4}' "${SALMON_DIR}/KM${i}_quant/quant.sf" >> "${TPM_DIR}/KM${i}.tpm"
done

# Join all samples into one TPM file
echo "  Joining all samples..."
$SUPPA joinFiles -f tpm \
    -i "${TPM_DIR}/KM1.tpm" "${TPM_DIR}/KM2.tpm" "${TPM_DIR}/KM3.tpm" \
       "${TPM_DIR}/KM4.tpm" "${TPM_DIR}/KM5.tpm" "${TPM_DIR}/KM6.tpm" \
       "${TPM_DIR}/KM7.tpm" "${TPM_DIR}/KM8.tpm" "${TPM_DIR}/KM9.tpm" \
    -o "${TPM_DIR}/TPM_all_samples"

# Fix header (remove duplicate sample name prefix)
sed -i '' '1s/.*/KM1\tKM2\tKM3\tKM4\tKM5\tKM6\tKM7\tKM8\tKM9/' "${TPM_DIR}/TPM_all_samples.tpm"

# Split TPM by condition
echo "  Splitting TPM by condition..."
awk 'BEGIN{OFS="\t"} NR==1{print "KM1","KM2","KM3"} NR>1{print $1,$2,$3,$4}' "${TPM_DIR}/TPM_all_samples.tpm" > "${TPM_DIR}/TPM_CT8.tpm"
awk 'BEGIN{OFS="\t"} NR==1{print "KM4","KM5","KM6"} NR>1{print $1,$5,$6,$7}' "${TPM_DIR}/TPM_all_samples.tpm" > "${TPM_DIR}/TPM_CT20.tpm"
awk 'BEGIN{OFS="\t"} NR==1{print "KM7","KM8","KM9"} NR>1{print $1,$8,$9,$10}' "${TPM_DIR}/TPM_all_samples.tpm" > "${TPM_DIR}/TPM_PerKO.tpm"

echo "  TPM files ready!"

# =============================================================================
# STEP 2: GENERATE EVENTS
# =============================================================================

echo ""
echo "=== STEP 2: Generating events (boundary=${BOUNDARY}, threshold=${THRESHOLD}nt) ==="

$SUPPA generateEvents \
    -i "${GTF}" \
    -o "${EVENTS_DIR}/events" \
    -f ioe \
    -e SE SS MX RI FL \
    -b "${BOUNDARY}" \
    -t "${THRESHOLD}"

echo "  Events generated!"

# =============================================================================
# STEP 3: PSI PER EVENT
# =============================================================================

echo ""
echo "=== STEP 3: Calculating PSI per event ==="

# Determine IOE file suffix based on boundary mode
if [ "$BOUNDARY" == "S" ]; then
    IOE_SUFFIX="strict"
else
    IOE_SUFFIX="variable_${THRESHOLD}"
fi

for EVENT in SE A3 A5 MX RI AF AL; do
    echo "  Processing ${EVENT}..."
    $SUPPA psiPerEvent \
        --ioe-file "${EVENTS_DIR}/events_${EVENT}_${IOE_SUFFIX}.ioe" \
        --expression-file "${TPM_DIR}/TPM_all_samples.tpm" \
        -f "${TPM_FILTER}" \
        -o "${PSI_DIR}/${EVENT}"
done

echo "  PSI files ready!"

# =============================================================================
# STEP 4: SPLIT PSI BY CONDITION
# =============================================================================

echo ""
echo "=== STEP 4: Splitting PSI files by condition ==="

for EVENT in SE A3 A5 MX RI AF AL; do
    awk 'BEGIN{OFS="\t"} NR==1{print "KM1","KM2","KM3"} NR>1{print $1,$2,$3,$4}' "${PSI_DIR}/${EVENT}.psi" > "${PSI_DIR}/${EVENT}_CT8.psi"
    awk 'BEGIN{OFS="\t"} NR==1{print "KM4","KM5","KM6"} NR>1{print $1,$5,$6,$7}' "${PSI_DIR}/${EVENT}.psi" > "${PSI_DIR}/${EVENT}_CT20.psi"
    awk 'BEGIN{OFS="\t"} NR==1{print "KM7","KM8","KM9"} NR>1{print $1,$8,$9,$10}' "${PSI_DIR}/${EVENT}.psi" > "${PSI_DIR}/${EVENT}_PerKO.psi"
done

echo "  PSI files split by condition!"

# =============================================================================
# STEP 5: DIFF SPLICE
# =============================================================================

echo ""
echo "=== STEP 5: Running diffSplice (method=${DIFF_METHOD}) ==="

cd "${OUTPUT_DIR}"

for EVENT in SE A3 A5 MX RI AF AL; do
    echo "  Processing ${EVENT}..."
    $SUPPA diffSplice \
        --method "${DIFF_METHOD}" \
        --input "events/events_${EVENT}_${IOE_SUFFIX}.ioe" \
        --psi "psi/${EVENT}_CT8.psi" "psi/${EVENT}_CT20.psi" "psi/${EVENT}_PerKO.psi" \
        --tpm "tpm/TPM_CT8.tpm" "tpm/TPM_CT20.tpm" "tpm/TPM_PerKO.tpm" \
        -gc \
        -o "diff/diff_${EVENT}"
done

echo "  diffSplice done!"

# =============================================================================
# STEP 6: CLUSTER EVENTS (SE only for now)
# =============================================================================

echo ""
echo "=== STEP 6: Clustering events ==="

# Clean psivec (remove nan rows)
awk 'NR==1 || !/nan/' "diff/diff_SE.dpsi.temp.0" > "diff/diff_SE_clean.dpsi"
awk 'NR==1 || !/nan/' "diff/diff_SE.psivec" > "diff/diff_SE_clean.psivec"

$SUPPA clusterEvents \
    --dpsi "diff/diff_SE_clean.dpsi" \
    --psivec "diff/diff_SE_clean.psivec" \
    --sig-threshold "${CLUSTER_SIG}" \
    --eps "${CLUSTER_EPS}" \
    --min-pts "${CLUSTER_MINPTS}" \
    -dt "${CLUSTER_DPSI}" \
    -g 1-3,4-6,7-9 \
    -c "${CLUSTER_METHOD}" \
    -o "cluster/cluster_SE"

echo "  Clustering done!"

# =============================================================================
# SUMMARY
# =============================================================================

echo ""
echo "=== PIPELINE COMPLETE ==="
echo "Results saved in: ${OUTPUT_DIR}"
echo ""
echo "Significant events (p < 0.05):"
echo "--- CT8 vs CT20 ---"
for EVENT in SE A3 A5 MX RI AF AL; do
    COUNT=$(awk 'NR>1 && $3!="nan" && $3+0 < 0.05' "diff/diff_${EVENT}.dpsi.temp.0" | wc -l)
    echo "  ${EVENT}: ${COUNT} events"
done
echo "--- CT20 vs PerKO ---"
for EVENT in SE A3 A5 MX RI AF AL; do
    COUNT=$(awk 'NR>1 && $3!="nan" && $3+0 < 0.05' "diff/diff_${EVENT}.dpsi.temp.1" | wc -l)
    echo "  ${EVENT}: ${COUNT} events"
done
