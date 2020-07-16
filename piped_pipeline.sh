#!/bin/bash
# TODO manage output during run. errors for error log and standart output to output log
# Echo to screen the pipeline's stage currently working
# !assuming in working directory!
# mkdir -p alignment BAM CNS fastq/{raw, trimmed} QC/fastqc refs test Trees
# (1) merge / sort reads ?

############ so far as general pipeline
# start pipe

# index reference
# export?
bwa index refs/REF_NC_045512.2.fasta
new_samtools=/data/software/samtools/samtools-1.10_new/samtools-1.10/samtools
new_bcftools=/data/software/bcftools/bcftools-1.10.2_new/bcftools-1.10.2/bcftools

for r1 in fastq/trimmed/*R1*.paired.fastq.gz; do
  r2=${r1/R1/R2} # ${var/find/replace}
  output=${r1/_R1/}
  bwa mem -v1 -t4 refs/REF_NC_045512.2.fasta $r1 $r2 | samtools view -@ 8 -Sb - | $new_samtools view -@ 8 -b -F 260 - |\
  samtools sort -@ 8 - BAM/`basename $output .paired.fastq.gz`.mapped.sorted
done


# index
for file in BAM/*.mapped.sorted.bam; do
  samtools index $file
done


# (7) create concensus sequence
# https://github.com/samtools/bcftools/wiki/HOWTOs#consensus-calling
# http://samtools.github.io/bcftools/howtos/consensus-sequence.html
for file in BAM/*.mapped.sorted.bam; do
  $new_samtools mpileup -uf refs/REF_NC_045512.2.fasta $file | $new_bcftools call -mv -Oz --threads 8 \
  -o CNS/calls.vcf.gz # change to bcftools mpileup??
  $new_bcftools index --threads 8 CNS/calls.vcf.gz
  cat refs/REF_NC_045512.2.fasta | $new_bcftools consensus CNS/calls.vcf.gz > \
  CNS/`basename $file .mapped.sorted.bam`.fasta
done

