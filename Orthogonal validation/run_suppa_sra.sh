#!/usr/bin/env bash
# This script takes the transcriptome structure (GTF), an SRA experiment accession number (SRPXXXXXX) of the quantified data and date of quantification as inputs and produces the PSI values per sample.

# USAGE: ./run_suppa_sra.sh <path/to/gtf/file> <SRP accession number> <date in yyyy-mm-dd>

OUTDIR="../results/${2}_suppa_$(date +%Y-%m-%d)"
LOGFILE="${OUTDIR}/suppa.log"

mkdir -p $OUTDIR

echo "using suppa v2.3" > $LOGFILE

suppa.py generateEvents -i $1  -o "${OUTDIR}/transcriptome_ext" -f ioe -e SE SS MX RI FL -b S && echo "Generating events complete" >> $LOGFILE

FIRST=1
while read line 
do
	if [ $FIRST -eq 1 ]; then
		awk -F '\t' 'NR>1 {print $1}' "../results/$2_salmon_${3}/${line}/quant.sf" > "${OUTDIR}/combined_tpm.tab"
		FIRST=0
		HEADER=$line
	else
		HEADER="${HEADER}\t${line}"

	fi
	awk -F '\t' 'NR>1 {print $4}' "../results/$2_salmon_${3}/${line}/quant.sf" > "${OUTDIR}/temp.tab"
	paste "${OUTDIR}/combined_tpm.tab" "${OUTDIR}/temp.tab" | column -s $'\t' -o $'\t'| sponge "${OUTDIR}/combined_tpm.tab"
done <<< $(cat "../meta/$2_Acc_List.txt")
sed -i -E 's/[_.][A-Za-z0-9:]+//' "${OUTDIR}/combined_tpm.tab"
sed -i "1i${HEADER}" "${OUTDIR}/combined_tpm.tab" && echo "Compiled TPM values for all samples." >> $LOGFILE
rm "${OUTDIR}/temp.tab"

for ioefile in "${OUTDIR}/transcriptome_ext_"*.ioe; do
	suppa.py psiPerEvent -i $ioefile \
			     -e "${OUTDIR}/combined_tpm.tab" \
			     -o "${OUTDIR}/$(basename -s .ioe $ioefile)" && echo "Computing PSI values for $(basename -s .ioe $ioefile) complete." >> $LOGFILE
done
