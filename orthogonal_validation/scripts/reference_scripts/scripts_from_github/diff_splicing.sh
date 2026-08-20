#!/usr/bin/env bash
# This script takes the transcriptome structure (GTF) and expression levels (TPM) as inputs and produces the PSI values per sample.

# USAGE: sh diff_splicing.sh 

TXDIR="../results/suppa_2024-11-06"
EXPDIR="../data"
OUTDIR="../results/Differential_suppa_$(date +%Y-%m-%d)"
LOGFILE="${OUTDIR}/suppa.log"

mkdir -p $OUTDIR

#Making expression files for CT 20 WT vs KO only no need now
#awk -v FS='\t' -v OFS='\t' '{print $1,$8,$9}' "${EXPDIR}/" > "${EXPDIR}/WT20.txt"
#awk -v FS='\t' -v OFS='\t' '{print $1,$14,$15}' "${EXPDIR}/" > "${EXPDIR}/KO20.txt"

#take abundance column x y to WT_20

for ioefile in "${TXDIR}/transcriptome_ext_"*.ioe; do
tmp=$(echo "$ioefile" | awk -F '_' '{print $4}')
  echo ${tmp}_WT 
  echo ${tmp}_KO
suppa.py psiPerEvent -i $ioefile -e "${EXPDIR}/WT_abundance.txt" -o "${OUTDIR}/${tmp}_WT"
suppa.py psiPerEvent -i $ioefile -e "${EXPDIR}/KO_abundance.txt" -o "${OUTDIR}/${tmp}_KO"
suppa.py diffSplice --method empirical -gc -p "${OUTDIR}/${tmp}_WT.psi" "${OUTDIR}/${tmp}_KO.psi" -e "${EXPDIR}/WT_abundance.txt" "${EXPDIR}/KO_abundance.txt" --input $ioefile -o "${OUTDIR}/res_${tmp}"
done 
