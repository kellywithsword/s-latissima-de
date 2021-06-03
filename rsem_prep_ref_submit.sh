#!/bin/bash
#SBATCH -p cegs
#SBATCH --mem=30gb
#SBATCH --time=30:00:00
#SBATCH -c 12

source ~/bin/anaconda3/etc/profile.d/conda.sh
conda activate rsem

reference=$1
reference_no_ext=`echo "$reference" | sed 's/\..*//g'`
sample_id=$2
echo "Aligning sample $sample_id to $reference_no_ext"

rsem-calculate-expression --bowtie2 --sort-bam-by-coordinate -p 12 --paired-end ${sample_id}_R1.fastq.gz ${sample_id}_R2.fastq.gz $reference_no_ext $sample_id

