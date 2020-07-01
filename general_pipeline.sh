#!/bin/bash

# prepare directories!!!! #TODO
# assuming in working directory
mkdir alignment BAM CNS fastq fastq/raw fastq/interleaved QC QC/fastqc refs SAM test Trees
cd fastq/raw # location of raw fastq files

# (1) merge / sort reads ?

# (2) QC (fastqc & multiqc) + trim (trimmomatic)
# fastqc
# export PATH=$PATH:/data/software/FastQC/fastqc ## not realy working
fastqc_new=/data/software/FastQC/fastqc
chmod 755 /data/software/FastQC/fastqc
$fastqc_new ./*.fastq.gz --outdir=QC/fastqc # NOTICE: EXPECTING fastq.gz POSTFIX
# multiqc
# export export PATH=$PATH:/data/software/multiqc/MultiQC/multiqc
multiqc QC/fastqc

# trimmomatic
# export PATH=$PATH:/data/software/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar
trimmomatic_path=/data/software/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar
for file in fastq/interleaved/*.fastq.gz; do
  java -jar $trimmomatic_path PE -phred33 -threads 32 $file fastq/interleaved_trimmed/`basename $file` TRAILING:28
done

# TODO: Maybe instead of running all samples together, for each sample create a pipe for better runtime?
# (3) Map reads to corona virus (REF_NC_045512.2)
# index reference # TODO should be in the pipeline or in a different script for first use?
# export?
bwa index refs/REF_NC_045512.2.fasta
# map reads to reference
for file in fastq/interleaved/*fastq.gz; do
  bwa mem refs/REF_NC_045512.2.fasta $file > SAM/`basename $file fastq.gz`.sam
done

# (4) sam to bam
for file in SAM/*.sam; do
  samtools view -Sb $file > BAM/`basename $file .sam`.bam
done

new_samtools=/data/software/samtools/samtools-1.10_new/samtools-1.10/samtools
# (5) keep only mapped reads
for file in BAM*/.bam; do
   $new_samtools view -b -F 260 $file > BAM/`basename $file bam`.mapped.bam
done

# (6) sort and index bam files
# sort
for $file in BAM/*.mapped.bam; do
  samtools sort $file BAM/`basename $file mapped.bam`.mapped.sorted
done
# index
for file in BAM/*.mapped.sorted.bam; do
  samtools index $file
done

# (7) create concensus sequence
new_bcftools=/data/software/bcftools/bcftools-1.10.2_new/bcftools-1.10.2/bcftools
# https://github.com/samtools/bcftools/wiki/HOWTOs#consensus-calling
# http://samtools.github.io/bcftools/howtos/consensus-sequence.html
for file in BAM/*.mapped.sorted.bam; do
  $new_samtools mpileup -uf   refs/REF_NC_045512.2.fasta $file | $new_bcftools call -mv -Oz -o CNS/calls.vcf.gz
  $new_bcftools index CNS/calls.vcf.gz
  cat refs/REF_NC_045512.2.fasta | $new_bcftools consensus CNS/calls.vcf.gz > CNS/`basename $file .mapped.sorted.bam`.fastq
done

# (8) fastq => fasta
chmod 755 /data/software/seqtk/seqtk
# export PATH=$PATH:/data/software/seqtk/seqtk
for file in CNS/*.fastq; do
  seqtk seq -a $file > CNS/`basename $file .fastq`.fasta
done
