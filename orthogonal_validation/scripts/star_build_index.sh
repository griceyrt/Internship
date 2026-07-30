#!/bin/bash
#SBATCH --job-name=star_index
#SBATCH --partition=Lake-short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=%j_star_index.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biotools

# Builds a STAR genome index using the genome FASTA already reconstructed
# during the 3'UTR extension step (genome_only.fa, pulled from gentrome_ext.fa's
# decoy sequences) and the same 3'UTR-extended GTF, so this run is directly
# comparable to that earlier test and replicates Asher's own STAR+htseq-count
# method as closely as possible.

DATA_DIR="$HOME/orthogonal_validation/data"
GENOME="$DATA_DIR/genome_only.fa"
GTF="$DATA_DIR/transcriptome_productivity_extended.gtf"
INDEX_DIR="$DATA_DIR/star_index_extended"

mkdir -p "$INDEX_DIR"

STAR --runMode genomeGenerate \
    --genomeDir "$INDEX_DIR" \
    --genomeFastaFiles "$GENOME" \
    --sjdbGTFfile "$GTF" \
    --sjdbOverhang 74 \
    --runThreadN "$SLURM_CPUS_PER_TASK"

echo "Done. Index at: $INDEX_DIR"
