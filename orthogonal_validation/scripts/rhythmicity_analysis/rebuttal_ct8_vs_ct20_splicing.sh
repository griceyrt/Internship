#!/bin/bash
# Rebuttal point (d) -- WT CT8 vs CT20 splicing analysis (SUPPA). Same pattern
# used to build data/nanopore_suppa/ (abundance_forsuppa.R + run_suppa.sh),
# adapted to compare WT-CT8 vs WT-CT20 instead of WT-CT20 vs KO-CT20.
# Author: Gricey
# GTF reused from boundary_analysis/, reference annotation only.
# USAGE: bash scripts/rebuttal_ct8_vs_ct20_splicing.sh
# Run from: orthogonal_validation/

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${SCRIPT_DIR}/.."
SALMON_DIR="${PROJECT_ROOT}/data/salmon_2024-09-25"
GTF="${PROJECT_ROOT}/../boundary_analysis/data/transcriptome_productivity.gtf"
OUTDIR="${PROJECT_ROOT}/results/rhythmicity_analysis/rebuttal_ct8_vs_ct20"
SUPPA="python3 ${PROJECT_ROOT}/../SUPPA/suppa.py"

DPSI_CUTOFF=0.1
PVAL_CUTOFF=0.05  # Nanopore convention (matches workflow_decisions.md), same
                  # threshold used for the original WT-CT20-vs-KO-CT20 comparison

mkdir -p "${OUTDIR}/events" "${OUTDIR}/psi" "${OUTDIR}/diff"

echo "=== Step 1: generateEvents (strict boundary, matching Khushi's own run_suppa.sh) ==="
$SUPPA generateEvents -i "$GTF" -o "${OUTDIR}/events/events" -f ioe -e SE SS MX RI FL -b S

echo "=== Step 2: build TPM table for CT8 + CT20 ==="
python3 - "$SALMON_DIR" "$OUTDIR" << 'EOF'
import sys
import pandas as pd

salmon_dir, outdir = sys.argv[1], sys.argv[2]
samples = {
    "CT8_1": f"{salmon_dir}/quant_CT8_1/quant_CT8_1.sf",
    "CT8_2": f"{salmon_dir}/quant_CT8_2/quant_CT8_2.sf",
    "CT20_1": f"{salmon_dir}/quant_CT20_1/quant_CT20_1.sf",
    "CT20_2": f"{salmon_dir}/quant_CT20_2/quant_CT20_2.sf",
}
tpm = pd.DataFrame()
for s, path in samples.items():
    df = pd.read_csv(path, sep="\t", index_col="Name")
    tpm[s] = df["TPM"]
tpm.index.name = None
tpm.to_csv(f"{outdir}/tpm_ct8_ct20.tab", sep="\t")
print(tpm.shape)
EOF
# strip pandas' leading-tab header (SUPPA misparses it as a field-count mismatch)
sed -i '1s/.*/CT8_1\tCT8_2\tCT20_1\tCT20_2/' "${OUTDIR}/tpm_ct8_ct20.tab"

echo "=== Step 3: psiPerEvent per event type ==="
for EVENT in SE A3 A5 MX RI AF AL; do
  $SUPPA psiPerEvent --ioe-file "${OUTDIR}/events/events_${EVENT}_strict.ioe" \
    --expression-file "${OUTDIR}/tpm_ct8_ct20.tab" \
    -o "${OUTDIR}/psi/${EVENT}"
done

echo "=== Step 4: split PSI + TPM by condition ==="
for EVENT in SE A3 A5 MX RI AF AL; do
  awk 'BEGIN{OFS="\t"} NR==1{print "CT8_1","CT8_2"} NR>1{print $1,$2,$3}' "${OUTDIR}/psi/${EVENT}.psi" > "${OUTDIR}/psi/${EVENT}_CT8.psi"
  awk 'BEGIN{OFS="\t"} NR==1{print "CT20_1","CT20_2"} NR>1{print $1,$4,$5}' "${OUTDIR}/psi/${EVENT}.psi" > "${OUTDIR}/psi/${EVENT}_CT20.psi"
done
awk 'BEGIN{OFS="\t"} NR==1{print "CT8_1","CT8_2"} NR>1{print $1,$2,$3}' "${OUTDIR}/tpm_ct8_ct20.tab" > "${OUTDIR}/tpm_CT8.tab"
awk 'BEGIN{OFS="\t"} NR==1{print "CT20_1","CT20_2"} NR>1{print $1,$4,$5}' "${OUTDIR}/tpm_ct8_ct20.tab" > "${OUTDIR}/tpm_CT20.tab"

echo "=== Step 5: diffSplice (empirical method, -gc, matching established project convention) ==="
for EVENT in SE A3 A5 MX RI AF AL; do
  $SUPPA diffSplice --method empirical \
    --input "${OUTDIR}/events/events_${EVENT}_strict.ioe" \
    --psi "${OUTDIR}/psi/${EVENT}_CT8.psi" "${OUTDIR}/psi/${EVENT}_CT20.psi" \
    --tpm "${OUTDIR}/tpm_CT8.tab" "${OUTDIR}/tpm_CT20.tab" \
    -gc \
    -o "${OUTDIR}/diff/diff_${EVENT}"
done

echo ""
echo "=== Summary: significant events (|dPSI|>${DPSI_CUTOFF} & p<${PVAL_CUTOFF}) ==="
for EVENT in SE A3 A5 MX RI AF AL; do
  COUNT=$(awk -v dc="$DPSI_CUTOFF" -v pc="$PVAL_CUTOFF" \
    'NR>1 && $2!="nan" && $3!="nan" && ($2>dc || $2<-dc) && $3<pc' "${OUTDIR}/diff/diff_${EVENT}.dpsi.temp.0" | wc -l)
  echo "  ${EVENT}: ${COUNT} events"
done

echo ""
echo "Done. Results in ${OUTDIR}/diff/"
