#!/bin/bash
# GSE133398 -- SUPPA pipeline (generateEvents, psiPerEvent, split by
# condition, diffSplice). WT = GFP control, KO = CRE / SRSF3 KO.
# Author: Gricey
# USAGE: bash scripts/GSE133398_run_suppa.sh
# Run from: orthogonal_validation/

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${SCRIPT_DIR}/.."

SUPPA="python3 ${SCRIPT_DIR}/../../SUPPA/suppa.py"
GTF="${BASE_DIR}/../boundary_analysis/data/transcriptome_productivity.gtf"
NORM_TABLE="${BASE_DIR}/results/GSE133398/normalisation_GSE133398/combined_norm.tab"
OUTDIR="${BASE_DIR}/results/GSE133398/suppa_GSE133398_$(date +%Y-%m-%d)"
mkdir -p "${OUTDIR}"

LOGFILE="${OUTDIR}/suppa.log"
echo "SUPPA run on $(date)" > "$LOGFILE"

WT_SAMPLES="WT_a\tWT_b\tWT_c\tWT_d"
KO_SAMPLES="KO_a\tKO_b\tKO_c\tKO_d"

echo "=== STEP 7: generateEvents ===" | tee -a "$LOGFILE"
$SUPPA generateEvents \
    -i "${GTF}" \
    -o "${OUTDIR}/transcriptome_ext" \
    -f ioe \
    -e SE SS MX RI FL \
    -b S \
    && echo "generateEvents complete." | tee -a "$LOGFILE"

if [ ! -f "${NORM_TABLE}" ]; then
    echo "ERROR: combined_norm.tab not found at ${NORM_TABLE}" | tee -a "$LOGFILE"
    echo "  Run GSE133398_normalisation_deseq2.R first." | tee -a "$LOGFILE"
    exit 1
fi
echo "=== STEP 8: Using expression table: ${NORM_TABLE} ===" | tee -a "$LOGFILE"

echo "=== STEP 9: psiPerEvent ===" | tee -a "$LOGFILE"
for ioefile in "${OUTDIR}/transcriptome_ext_"*.ioe; do
    EVENT=$(basename -s .ioe "$ioefile" | sed 's/transcriptome_ext_//')
    echo "  Processing ${EVENT}..." | tee -a "$LOGFILE"
    $SUPPA psiPerEvent \
        -i "${ioefile}" \
        -e "${NORM_TABLE}" \
        -o "${OUTDIR}/${EVENT}" \
        && echo "  psiPerEvent ${EVENT} complete." | tee -a "$LOGFILE"
done

echo "=== STEP 10a: Splitting by condition ===" | tee -a "$LOGFILE"
echo -e "transcript_id\t${WT_SAMPLES}" > "${OUTDIR}/WT_norm.tab"
awk 'NR>1 {print $1"\t"$2"\t"$3"\t"$4"\t"$5}' "${NORM_TABLE}" >> "${OUTDIR}/WT_norm.tab"

echo -e "transcript_id\t${KO_SAMPLES}" > "${OUTDIR}/KO_norm.tab"
awk 'NR>1 {print $1"\t"$6"\t"$7"\t"$8"\t"$9}' "${NORM_TABLE}" >> "${OUTDIR}/KO_norm.tab"

echo "Expression tables split." | tee -a "$LOGFILE"

for psifile in "${OUTDIR}/"*_strict.psi; do
    [ -f "$psifile" ] || continue
    EVENT=$(basename -s .psi "$psifile")

    awk 'NR==1{print $1"\tWT_a\tWT_b\tWT_c\tWT_d"} NR>1{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' \
        "${psifile}" > "${OUTDIR}/${EVENT}_WT.psi"

    awk 'NR==1{print $1"\tKO_a\tKO_b\tKO_c\tKO_d"} NR>1{print $1"\t"$6"\t"$7"\t"$8"\t"$9}' \
        "${psifile}" > "${OUTDIR}/${EVENT}_KO.psi"
done

echo "PSI files split by condition." | tee -a "$LOGFILE"

echo "=== STEP 10b: diffSplice ===" | tee -a "$LOGFILE"
DIFF_DIR="${OUTDIR}/diff"
mkdir -p "${DIFF_DIR}"

cd "${OUTDIR}"

for ioefile in transcriptome_ext_*.ioe; do
    EVENT=$(basename -s .ioe "$ioefile" | sed 's/transcriptome_ext_//')
    WT_PSI="${EVENT}_WT.psi"
    KO_PSI="${EVENT}_KO.psi"

    if [ ! -f "${WT_PSI}" ] || [ ! -f "${KO_PSI}" ]; then
        echo "  Skipping ${EVENT} — PSI files not found." | tee -a "$LOGFILE"
        continue
    fi

    echo "  diffSplice ${EVENT}..." | tee -a "$LOGFILE"
    $SUPPA diffSplice \
        --method empirical \
        --input "${ioefile}" \
        --psi "${WT_PSI}" "${KO_PSI}" \
        --tpm "WT_norm.tab" "KO_norm.tab" \
        -gc \
        -o "diff/diff_${EVENT}" \
        && echo "  diffSplice ${EVENT} complete." | tee -a "$LOGFILE"
done

echo "" | tee -a "$LOGFILE"
echo "=== SUPPA PIPELINE COMPLETE ===" | tee -a "$LOGFILE"
echo "Results in: ${OUTDIR}" | tee -a "$LOGFILE"
