#!/bin/bash
#SBATCH -p cegs
#SBATCH --mem=20gb
#SBATCH --time=03:00:00

#source ~/bin/anaconda3/etc/profile.d/conda.sh
source activate rsem

reference=$1
reference_no_ext=`echo "$reference" | sed 's/\..*//g'`
echo "Preparing reference $reference_no_ext"

rsem-prepare-reference --bowtie2 $reference $reference_no_ext


