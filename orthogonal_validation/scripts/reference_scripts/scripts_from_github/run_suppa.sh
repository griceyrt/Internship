#!/usr/bin/env sh
# This script takes the transcriptome structure (GTF) and expression levels (TPM) as inputs and produces the PSI values per sample.
# USAGE: sh run_suppa.sh 
TXDIR="../results/FLAIR_2024-09-11"
TPMDIR="../results/salmon_2024-09-25"
OUTDIR="../results/suppa_$(date +%Y-%m-%d)"
LOGFILE="${OUTDIR}/suppa.log"

. /home/kmamgain/miniconda3/etc/profile.d/conda.sh
conda activate /home/kmamgain/miniconda3/envs/suppa

mkdir -p $OUTDIR
echo "using suppa v2.3" > $LOGFILE
suppa.py generateEvents -i "${TXDIR}/transcriptome_ext.gtf" -o "${OUTDIR}/transcriptome_ext" -f ioe -e SE SS MX RI FL -b S && echo "Generating events complete" >> $LOGFILE

FIRST=1
for sffile in "${TPMDIR}/quant_"*"/quant_"*.sf; do
	if [ $FIRST -eq 1 ]; then
		awk 'NR>1 {print $1}' $sffile > "${OUTDIR}/combined_tpm.tab"
		FIRST=0
		HEADER=$(basename -s .sf $sffile)
	else
		HEADER="${HEADER}\t$(basename -s .sf $sffile)"

	fi
	awk 'NR>1 {print $4}' $sffile > "${OUTDIR}/temp.tab"
	paste "${OUTDIR}/combined_tpm.tab" "${OUTDIR}/temp.tab" | column -s $'\t' -t | sponge "${OUTDIR}/combined_tpm.tab"
done
sed -i -E 's/[_.][A-Za-z0-9:]+//' "${OUTDIR}/combined_tpm.tab"
sed -i "1i$HEADER" "${OUTDIR}/combined_tpm.tab" && echo "Compiled TPM values for all samples." >> $LOGFILE
rm "${OUTDIR}/temp.tab"

for ioefile in "${OUTDIR}/transcriptome_ext_"*.ioe; do
	suppa.py psiPerEvent -i $ioefile \
			     -e "../data/abundance_no_batch_for_suppa.txt"\
			     -o "${OUTDIR}/$(basename -s .ioe $ioefile)" && echo "Computing PSI values for $(basename -s .ioe $ioefile) complete." >> $LOGFILE
done