#!/bin/bash
#SBATCH -p cegs
#SBATCH --mem=8gb
#SBATCH --time=10:00:00
#SBATCH -c 12

reference=$1
bam=$2
bam_base=`sed 's/\..*//g' $bam`

picard CollectMultipleMetrics -I $bam -R $reference -O ${bam_basename}.picard_metrics

