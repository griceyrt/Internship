#!/bin/bash
# =============================================================================
# SUPPA Pipeline — Orthogonal Validation Project 2
# Author: Gricey
# Adapted from: BHARATH_run_suppa_sra.sh + boundary_analysis/scripts/run_suppa.sh
#
# Steps covered:
#   7.  generateEvents from GTF (strict boundary -b S)
#   9.  psiPerEvent (uses normalized expression table from R Step 5)
#   10a. Split PSI and expression table by condition (WT / KO)
#   10b. diffSplice (WT CT16-20 vs KO CT16-20, empirical method)
#
# USAGE: bash scripts/run_suppa.sh
# Run from: orthogonal_validation/
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${SCRIPT_DIR}/.."

# SUPPA
SUPPA="python3 ${SCRIPT_DIR}/../../SUPPA/suppa.py"

# Input files
GTF="${BASE_DIR}/../boundary_analysis/data/transcriptome_productivity.gtf"
NORM_TABLE="${BASE_DIR}/results/normalisation/combined_norm.tab"

# Output directory
OUTDIR="${BASE_DIR}/results/suppa_$(date +%Y-%m-%d)"
mkdir -p "${OUTDIR}"

LOGFILE="${OUTDIR}/suppa.log"
echo "SUPPA run on $(date)" > $LOGFILE

# Sample column names (must match headers in combined_norm.tab)
WT_SAMPLES="WT_a\tWT_b\tWT_c\tWT_d"
KO_SAMPLES="KO_a\tKO_b\tKO_c\tKO_d"

# =============================================================================
# STEP 7 — generateEvents from GTF
# =============================================================================
echo "=== STEP 7: generateEvents ===" | tee -a $LOGFILE

$SUPPA generateEvents \
    -i "${GTF}" \
    -o "${OUTDIR}/transcriptome_ext" \
    -f ioe \
    -e SE SS MX RI FL \
    -b S \
    && echo "generateEvents complete." | tee -a $LOGFILE

# =============================================================================
# STEP 8 — SUPPA expression input is already combined_norm.tab from R Step 5
# No action needed here, just confirming the file exists
# =============================================================================
if [ ! -f "${NORM_TABLE}" ]; then
    echo "ERROR: combined_norm.tab not found at ${NORM_TABLE}" | tee -a $LOGFILE
    exit 1
fi
echo "=== STEP 8: Using expression table: ${NORM_TABLE} ===" | tee -a $LOGFILE

# =============================================================================
# STEP 9 — psiPerEvent
# =============================================================================
echo "=== STEP 9: psiPerEvent ===" | tee -a $LOGFILE

for ioefile in "${OUTDIR}/transcriptome_ext_"*.ioe; do
    EVENT=$(basename -s .ioe $ioefile | sed 's/transcriptome_ext_//')
    echo "  Processing ${EVENT}..." | tee -a $LOGFILE
    $SUPPA psiPerEvent \
        -i "${ioefile}" \
        -e "${NORM_TABLE}" \
        -o "${OUTDIR}/${EVENT}" \
        && echo "  psiPerEvent ${EVENT} complete." | tee -a $LOGFILE
done

# =============================================================================
# STEP 10a — Split PSI and expression table by condition (WT / KO)
# =============================================================================
echo "=== STEP 10a: Splitting by condition ===" | tee -a $LOGFILE

# Split expression table
echo -e "transcript_id\t${WT_SAMPLES}" > "${OUTDIR}/WT_norm.tab"
awk 'NR>1 {print $1"\t"$2"\t"$3"\t"$4"\t"$5}' "${NORM_TABLE}" >> "${OUTDIR}/WT_norm.tab"

echo -e "transcript_id\t${KO_SAMPLES}" > "${OUTDIR}/KO_norm.tab"
awk 'NR>1 {print $1"\t"$6"\t"$7"\t"$8"\t"$9}' "${NORM_TABLE}" >> "${OUTDIR}/KO_norm.tab"

echo "Expression tables split." | tee -a $LOGFILE

# Split PSI files by condition
for psifile in "${OUTDIR}/"*_strict.psi; do
    [ -f "$psifile" ] || continue
    EVENT=$(basename -s .psi $psifile)

    # WT: columns 1 (event_id), 2-5 (WT_a,b,c,d)
    awk 'NR==1{print $1"\tWT_a\tWT_b\tWT_c\tWT_d"} NR>1{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' \
        "${psifile}" > "${OUTDIR}/${EVENT}_WT.psi"

    # KO: columns 1 (event_id), 6-9 (KO_a,b,c,d)
    awk 'NR==1{print $1"\tKO_a\tKO_b\tKO_c\tKO_d"} NR>1{print $1"\t"$6"\t"$7"\t"$8"\t"$9}' \
        "${psifile}" > "${OUTDIR}/${EVENT}_KO.psi"
done

echo "PSI files split by condition." | tee -a $LOGFILE

# =============================================================================
# STEP 10b — diffSplice (WT vs KO, empirical, gene-corrected)
# =============================================================================
echo "=== STEP 10b: diffSplice ===" | tee -a $LOGFILE

DIFF_DIR="${OUTDIR}/diff"
mkdir -p "${DIFF_DIR}"

cd "${OUTDIR}"

for ioefile in transcriptome_ext_*.ioe; do
    EVENT=$(basename -s .ioe $ioefile | sed 's/transcriptome_ext_//')
    WT_PSI="${EVENT}_WT.psi"
    KO_PSI="${EVENT}_KO.psi"

    if [ ! -f "${WT_PSI}" ] || [ ! -f "${KO_PSI}" ]; then
        echo "  Skipping ${EVENT} — PSI files not found." | tee -a $LOGFILE
        continue
    fi

    echo "  diffSplice ${EVENT}..." | tee -a $LOGFILE
    $SUPPA diffSplice \
        --method empirical \
        --input "${ioefile}" \
        --psi "${WT_PSI}" "${KO_PSI}" \
        --tpm "WT_norm.tab" "KO_norm.tab" \
        -gc \
        -o "diff/diff_${EVENT}" \
        && echo "  diffSplice ${EVENT} complete." | tee -a $LOGFILE
done

echo "" | tee -a $LOGFILE
echo "=== SUPPA PIPELINE COMPLETE ===" | tee -a $LOGFILE
echo "Results in: ${OUTDIR}" | tee -a $LOGFILE
