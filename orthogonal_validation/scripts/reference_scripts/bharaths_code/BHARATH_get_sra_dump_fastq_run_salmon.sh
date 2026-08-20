#!/bin/bash
# This script takes an SRA experiment accession (SRPXXXXXX) of a paired RNA-seq experiment and 
# retrieves it using sratools and quantifies it with Salmon using the chosen transcriptome index 
#
# USAGE: ./get_sra_dump_fastq_run_salmon.sh <SRA experiment accession> <path/to/salmon/index>

GENTROME="/extra/Mm/gentrome.fa"
DECOY="/extra/Mm/decoys.txt"
SALMON_INDEX="/extra/Mm/transcriptome_ext_index"
OUTDIR="../results/${1}_salmon_$(date +"%Y-%m-%d")"
LOG="../results/${1}_salmon_$(date +"%Y-%m-%d")/${1}.log"
CORES=7

if [ -z $2 ]; then
	salmon index -t $GENTROME \
		-i $SALMON_INDEX \
		--decoys $DECOY \
		-k 31 \
		--keepDuplicates
else 
	SALMON_INDEX=$2
fi

mkdir -p $OUTDIR
echo "Run on $(date)" > $LOG
salmon --version >> $LOG
parallel-fastq-dump --version >> $LOG

cat "../meta/$1_Acc_List.txt" | while read line
do
	# Fetch the .sra files if not available in <raw_data> directory
	prefetch -O ../raw_data/ $line \
		>> $LOG

	# Extract the fastq files if they are not present in the <temp> directory
	if [ ! -f "../temp/${line}_1.fastq.gz" ] || [ ! -f "../temp/${line}_2.fastq.gz" ]; then
		parallel-fastq-dump --threads 6 \
			--outdir "../temp/" \
			--split-files \
			--gzip \
			-s "${line}" \
			>> $LOG
	fi

	# Run salmon quantification (assumes indexing of txome is already done)
	salmon quant -i $SALMON_INDEX \
		-l A --validateMappings \
		-p $CORES \
		--seqBias --gcBias -q \
		-1 "../temp/${line}_1.fastq.gz" \
		-2 "../temp/${line}_2.fastq.gz" \
		-o "${OUTDIR}/${line}/" \
		>> $LOG
done



