#!/bin/bash

# prepare directories!!!! #TODO
# !assuming in working directory!
# mkdir alignment BAM CNS fastq fastq/raw fastq/interleaved QC QC/fastqc refs SAM test Trees
# mkdir QC/fastqc
# (1) merge / sort reads ?

# (2) QC (fastqc & multiqc) + trim (trimmomatic)
# fastqc
# export PATH=$PATH:/data/software/FastQC/fastqc ## not realy working
# new_fastqc=/data/software/FastQC/fastqc <- gives error. 'fastqc' command works.
chmod 755 $new_fastqc
$new_fastqc fastq/raw/*.fastq.gz --outdir=QC/fastqc # NOTICE: EXPECTING fastq.gz POSTFIX
# multiqc
# export export PATH=$PATH:/data/software/multiqc/MultiQC/multiqc
multiqc QC/fastqc