# RUN RNA SEQ COMPLETE PIPELINE.

threads=2
num_cores=$(($(nproc) - 2))
max_parallel=$((num_cores/threads))
echo $max_parallel
echo "--------------------"

############## QUALITY CONTROL ####################
# Fastqc and MultiQC
working_dir="/mnt/d/Kiran/organoide/rna_seq"
fastq_dir="${working_dir}/data/Fastq_files"
qc_dir="${working_dir}/results/qc"
fastqc -t $nproc -o $qc_dir ${fastq_dir}/*.fastq.gz
multiqc $qc_dir -o $qc_dir -n raw_mqc_report

trim poor quality with fastp
fastp_qc(){
    fastq_file="$1"
    basename="${fastq_file%_R1_001.fastq.gz}"
    sample_name="$(basename $basename)"
    outfile="${working_dir}/results/fastp/${sample_name}_trimmed.fastq.gz"
    echo "fastp processing ${sample_name}"
    echo "output dirctory is ${outfile}"
    fastp \
        -i $fastq_file \
        -o $outfile \
        --trim_poly_g \
        --qualified_quality_phred 28 \
        --length_required 35 \
        --cut_front --cut_tail \
        --cut_window_size 4 \
        --cut_mean_quality 28 \
        --thread $threads \
        -h ${working_dir}/results/logs/${sample_name}.html
}

export fastq_dir working_dir
export -f fastp_qc
export threads

find "$fastq_dir" -name "*.fastq.gz" | \
    parallel -j $max_parallel --bar fastp_qc

echo "fastp qc complete"
echo "-----------------------"


echo "Rerunning fastqc and multiqc"
fastp_dir="${working_dir}/results/fastp"
fastqc -t $num_cores -o ${fastp_dir} ${fastp_dir}/*.fastq.gz
multiqc $fastp_dir -o $fastp_dir -n trimmed_mqc_report

################# INFER STRANDEDNESS ##############
echo "Find read strand"
data_dir="${working_dir}/data"
gtfToGenePred $data_dir/Mus_musculus.GRCm39.115.gtf /dev/stdout | \
 genePredToBed /dev/stdin $data_dir/genes.bed

# check from bam file aligned by psi
psi_bam_dir="/mnt/d/Kiran/organoide/rna_seq/Projet_70_25/Analyses/Lib1_70_25_S1/STAR/Lib1_70_25_S1.bam"
infer_experiment.py -i $psi_bam_dir -r $data_dir/genes.bed -s 200000
######### Output results from reqseq ################
# This is SingleEnd Data
# Fraction of reads failed to determine: 0.1097
# Fraction of reads explained by "++,--": 0.8858
# Fraction of reads explained by "+-,-+": 0.0044
############## SO DATA IS SE FORWARD STRANDED ##########



###################### ALIGNMENT ####################
#   First build index and align one bam file to infer from it.
echo "Generating STAR indexes..."
STAR --runMode genomeGenerate \
  --genomeDir ${data_dir}/index/ \
  --genomeFastaFiles ${data_dir}/Mus_musculus.GRCm39.dna.primary_assembly.fa \
  --sjdbGTFfile ${data_dir}/Mus_musculus.GRCm39.115.gtf \
  --sjdbOverhang 75 \
  --runThreadN $num_cores

### Run aligner parrallely
threads=6
max_parallel=4
alignment(){
    fastq_file="$1"
    basename="${fastq_file%.fastq.gz}"
    sample_name="$(basename $basename)"
    outfile_prefix="${working_dir}/results/star/${sample_name}_"
    echo "Mapping file is ${sample_name}"
    echo "output file is ${outfile_prefix}"
    STAR --runMode alignReads \
        --genomeDir ${data_dir}/index/ \
        --readFilesIn $fastq_file \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --twopassMode Basic \
        --outFileNamePrefix $outfile_prefix \
        --outTmpDir /tmp/STAR_${sample_name} \
        --runThreadN $threads
    
    rm -rf /tmp/STAR_${sample_name}
    samtools index ${outfile_prefix}Aligned.sortedByCoord.out.bam
}
export fastp_dir working_dir data_dir
export -f alignment
export threads

find "$fastp_dir" -name "*.fastq.gz" | \
    parallel -j $max_parallel --bar alignment

# Quick summary of all STAR logs
star_dir="${working_dir}/results/star"
multiqc $star_dir -o $star_dir -n alignment_report

###################### Salmon ####################
### I choosed ot use allignment free sofware to work directly with the 
### transcriptome and also because it's less computation generating 
### larger index files and bam files which I don't need downstream.
### so you can comment the START aligner above.

# Build Salmon index (do once)
# For selective alignment (recommended), create a decoy-aware index
echo "Building Salmon index"
grep '^>' ${data_dir}/Mus_musculus.GRCm39.dna.primary_assembly.fa | cut -d ' ' -f1 | \
    sed 's/>//' > ${data_dir}/decoys.txt
cat ${data_dir}/Mus_musculus.GRCm39.cdna.all.fa \
    ${data_dir}/Mus_musculus.GRCm39.dna.primary_assembly.fa > ${data_dir}/gentrome.fa

salmon index -t ${data_dir}/gentrome.fa -d ${data_dir}/decoys.txt \
    -i ${data_dir}/salmon_index -p $num_cores

# Quantify each sample (SE reads)
echo "Running Salmon..."
echo "Parallelly running ${max_parallel} jobs with ${threads} cores each"
salmon_reads_counts(){
    fastq_file="$1"
    basename="${fastq_file%_trimmed.fastq.gz}"
    sample_name="$(basename $basename)"
    outfile_dir="${working_dir}/results/salmon/${sample_name}"
    echo "Salmon is processing file : ${sample_name}"
    echo "output dir is ${outfile_dir}"

    salmon quant \
        -i ${data_dir}/salmon_index \
        -l SF \
        -r $fastq_file \
        --seqBias --gcBias \
        --validateMappings \
        -p $threads \
        -o $outfile_dir/
}
export fastp_dir working_dir data_dir
export -f salmon_reads_counts
export threads

find "$fastp_dir" -name "*.fastq.gz" | \
    parallel -j $max_parallel --bar salmon_reads_counts


echo "Running comprehensive multiqc"
multiqc \
  $fastq_dir $fastp_dir ${working_dir}/results/logs/ \
  ${working_dir}/results/salmon -o ${working_dir}/results/multiqc/ \
  -n complete_qc_report