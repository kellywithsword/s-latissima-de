#!/bin/bash
#SBATCH -p cegs
#SBATCH --mem=15gb
#SBATCH --time=10:00:00
#SBATCH -c 12

bam=$1
bam_basename=`echo $bam; sed 's/.*\.//g'`

qualimap bamqc -nt 12 -sd -sdmode 1 -bam $bam -outdir bamqc_${bam_basename} --java-mem-size=14G

