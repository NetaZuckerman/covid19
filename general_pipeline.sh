#!/bin/bash

# TODO prepare directories
# !assuming in working directory!
# mkdir alignment BAM CNS fastq fastq/raw fastq/trimmed QC QC/fastqc refs SAM test Trees
# mkdir QC/fastqc
# TODO: Maybe instead of running all samples together, for each sample create a pipe for better runtime?
# (1) merge / sort reads ?

# (2) QC (fastqc & multiqc) + trim (trimmomatic)
# fastqc
# export PATH=$PATH:/data/software/FastQC/fastqc ## not realy working
# new_fastqc=/data/software/FastQC/fastqc <- gives error. 'fastqc' command works.
chmod 755 $new_fastqc
$new_fastqc fastq/raw/*.fastq.gz --outdir=QC/fastqc # NOTICE: EXPECTING fastq.gz POSTFIX
# multiqc
# export export PATH=$PATH:/data/software/multiqc/MultiQC/multiqc
multiqc QC/fastqc -o QC/

# trimmomatic
# export PATH=$PATH:/data/software/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar
# trimmomatic expect existing output folder
trimmomatic_path=/data/software/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar
for r1 in fastq/raw/*R1*; do
	r2=${r1/R1/R2}
	r1_base=$(basename -s .fastq.gz $r1)
	r2_base=$(basename -s .fastq.gz $r2)
	singles1=fastq/trimmed/"$r1_base".unpaired.fastq.gz
	singles2=fastq/trimmed/"$r2_base".unpaired.fastq.gz
	paired1=fastq/trimmed/"$r1_base".paired.fastq.gz
	paired2=fastq/trimmed/"$r2_base".paired.fastq.gz
	java -jar $trimmomatic_path PE -phred33 -threads 32 $r1 $r2 $paired1 $singles1 $paired2 $singles2 TRAILING:28
done
########################################################################################################################
# TODO: copy reference sequence from /data/projects/Michal/CoV2019/data/Eritreas/refs # already indexed, may skip index

# (3) Map reads to corona virus (REF_NC_045512.2)
# index reference # TODO should be in the pipeline or in a different script for first use?
# export?
bwa index refs/REF_NC_045512.2.fasta
# map reads to reference -> PE
for r1 in fastq/trimmed/*R1*.paired.fastq.gz; do
  r2=${r1/R1/R2} # ${var/find/replace}
  output=${r1/_R1/}
  output=${output/.paired/}
  bwa mem -v1 -t4 refs/REF_NC_045512.2.fasta $r1 $r2 > SAM/`basename $output fastq.gz`.sam
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
for file in BAM/*.mapped.bam; do
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
