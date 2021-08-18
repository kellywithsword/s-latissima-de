#!/bin/bash
#SBATCH -p cegs
#SBATCH --mem=8gb
#SBATCH --time=1-0
#SBATCH -c 12

sample_id=$1

source activate trimgalore
trim_galore --fastqc --paired --retain_unpaired --phred33 --output_dir trimmed_reads --length 36 --stringency 1 -e 0.1 ${sample_id}_R1.fastq.gz ${sample_id}_R2.fastq.gz

