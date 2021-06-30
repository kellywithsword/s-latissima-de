#!/bin/bash
#SBATCH -p cegs
#SBATCH --mem=20gb
#SBATCH --time=03:00:00

# Help message
if [[ $1 = "-h" ]] || [[ $1 = "--help" ]] \
|| [[ -z $1 ]]; then
  printf "Usage: sbatch rsem_prep_ref_submit.sh \
reference.fa \
[--gff3/--gtf/--transcript-to-gene-map/--... \
reference.(gff3/gtf/txt/...)]\n\n\
Note: FASTA & annotation files must be unzipped.\n\n\
Requires RSEM.\n\
(https://github.com/deweylab/RSEM)\n"
  exit 0
fi

source activate rsem

reference=$1
reference_no_ext=`echo "$reference" | sed 's/\..*//g'`
printf "Preparing reference $reference_no_ext.\n"

if [[ -z $2 ]]; then
  printf "Treating reference FASTA $reference \
as transcripts.\n"
  rsem-prepare-reference --bowtie2 $reference \
$reference_no_ext
elif [[ -z $3 ]]; then
  printf "Error: no annotation file provided.\n"
  exit 1  
else
  transcript_extract_flag=$2
  transcript_extract_file=$3
  printf "Treating reference $reference as genome, \
extracting transcripts with $transcript_extract_file.\n"
  if [[ $transcript_extract_flag = "--gff3" ]]; then
    rsem-prepare-reference --bowtie2 \
$transcript_extract_flag $transcript_extract_file \
--gff3-RNA-patterns mRNA,rRNA \
$reference $reference_no_ext
  else
    rsem-prepare-reference --bowtie2 \
$transcript_extract_flag $transcript_extract_file \
$reference $reference_no_ext
  fi
fi

