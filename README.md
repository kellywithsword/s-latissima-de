# s-latissima-de
## Differential expression analysis of _S. latissima_ RNAseq data


### What you'll need
- A reference genome or transcriptome (FASTA)
- RNAseq read files (FASTQ)

### 1. Create "samples file" containing list of sample IDs

For example, with paired-end samples whose FASTQ files are named 
*sampleid\_R1.fastq* and *sampleid\_R2.fastq*, you can run this
command in the directory containing the FASTQ files to create 
a compatible samples file:

`ls *fastq* | sed 's/_R[12].*//g' | sort -u > samples_file.txt`

### 2. Prepare reference for RSEM

`sbatch rsem_prep_ref_submit.sh reference.fasta \
[--gff3/--gtf/--transcript-to-gene-map/--... \
reference.(gff3/gtf/txt/...)`

