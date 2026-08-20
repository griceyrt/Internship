#! usr/bin/sh
outdir="../results/FLAIR_2024-09-11" 

# make an extended transcriptome fasta file
#grep -F ">" "${outdir}/collapse_2024-09-11.isoforms.fa" \
#       	| sed -E -e 's/_.*$//' -e 's/(\-[0-9])*$//' \
#	| cut -c2- \
#	| sort \
#	| uniq -u \
#	| grep -v -E "^ENSMUST" > "${outdir}/novel_fasta_tx_list.txt"
#status=$?

echo ${outdir}
#grep -F -f "../results/FLAIR_2024-09-11/novel_fasta_tx_list.txt" "../results/FLAIR_2024-09-11/collapse_2024-09-11.isoforms.gtf" | cat "../data/Mus_musculus.GRCm39.108.chr.gtf" - > "../results/FLAIR_2024-09-11/transcriptome_productivity.gtf"
#conda activate productivity
#gtfToGenePred "../results/FLAIR_2024-09-11/transcriptome_productivity.gtf" "../data/transcriptome_productivity.genepred"
genePredToBed "../data/transcriptome_productivity.genepred" "../data/for_productivity.bed"
Python /mnt/c/Users/kmamgain/Desktop/T/flair/flair-master/bin/predictProductivity.py -i "../data/for_productivity.bed"  -f "../data/Mus_musculus.GRCm39.dna.primary_assembly.fa"  -g "../results/FLAIR_2024-09-11/transcriptome_productivity.gtf" --longestORF --append_column > "../results/FLAIR_2024-09-11/productivity.txt"
