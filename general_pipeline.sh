#!/bin/bash

# how to get input files??
# !assuming in working directory!
# mkdir -p alignment BAM/mapped CNS fastq/{raw, trimmed} QC/fastqc refs SAM test Trees
# TODO: for each sample create a pipe for better runtime
# (1) merge / sort reads ?

# (2) QC (fastqc & multiqc) + trim (trimmomatic)
# fastqc
# export PATH=$PATH:/data/software/FastQC/fastqc ## not realy working
# new_fastqc=/data/software/FastQC/fastqc <- gives error. 'fastqc' command works.
# chmod 755 $new_fastqc
# fastqc fastq/raw/*.fastq.gz --outdir=QC/fastqc # NOTICE: EXPECTING fastq.gz POSTFIX
# multiqc
# export export PATH=$PATH:/data/software/multiqc/MultiQC/multiqc
# multiqc QC/fastqc -o QC/

# trimmomatic
# export PATH=$PATH:/data/software/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar
# trimmomatic expect existing output folder
#trimmomatic_path=/data/software/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar
#for r1 in fastq/raw/*R1*; do
#	r2=${r1/R1/R2}
#	r1_base=$(basename -s .fastq.gz $r1)
#	r2_base=$(basename -s .fastq.gz $r2)
#	singles1=fastq/trimmed/"$r1_base".unpaired.fastq.gz
#	singles2=fastq/trimmed/"$r2_base".unpaired.fastq.gz
#	paired1=fastq/trimmed/"$r1_base".paired.fastq.gz
#	paired2=fastq/trimmed/"$r2_base".paired.fastq.gz
#	java -jar $trimmomatic_path PE -phred33 -threads 32 $r1 $r2 $paired1 $singles1 $paired2 $singles2 TRAILING:28
#done

# (3) Map reads to corona virus (REF_NC_045512.2)
# index reference
# export?
bwa index refs/REF_NC_045512.2.fasta
# map reads to reference -> PE
# TODO: May be ran on trimmed or on raw. get wanter option form user!!
for r1 in fastq/trimmed/*R1*.paired.fastq.gz; do
  r2=${r1/R1/R2} # ${var/find/replace}
  output=${r1/_R1/}
  bwa mem -v1 -t4 refs/REF_NC_045512.2.fasta $r1 $r2 > SAM/`basename $output .paired.fastq.gz`.sam
done

# samtools <command> -@: number of threads (default 1)
# (4) sam to bam
for file in SAM/*.sam; do
  samtools view -@ 8 -Sb $file > BAM/`basename $file .sam`.bam
done

new_samtools=/data/software/samtools/samtools-1.10_new/samtools-1.10/samtools
new_bcftools=/data/software/bcftools/bcftools-1.10.2_new/bcftools-1.10.2/bcftools
# (5) keep only mapped reads
for file in BAM*/.bam; do
   $new_samtools view -@ 8 -b -F 260 $file > BAM/`basename $file .bam`.mapped.bam
done

# (6) sort and index bam files
# sort
for file in BAM/*.mapped.bam; do
  samtools sort -@ 8 $file BAM/`basename $file .mapped.bam`.mapped.sorted
done

# index
for file in BAM/*.mapped.sorted.bam; do
  samtools index $file
done

# (7) create concensus sequence
# https://github.com/samtools/bcftools/wiki/HOWTOs#consensus-calling
# http://samtools.github.io/bcftools/howtos/consensus-sequence.html
# (the commands)
#/data/software/samtools/samtools-1.10_new/samtools-1.10/samtools mpileup -uf refs/REF_NC_045512.2.fasta BAM/last/s11.mapped.sorted.bam | /data/software/bcftools/bcftools-1.10.2_new/bcftools-1.10.2/bcftools call -mv -Oz -o CNS/last/calls.vcf.gz
#/data/software/bcftools/bcftools-1.10.2_new/bcftools-1.10.2/bcftools index CNS/last/calls.vcf.gz
#cat refs/REF_NC_045512.2.fasta | /data/software/bcftools/bcftools-1.10.2_new/bcftools-1.10.2/bcftools consensus CNS/last/calls.vcf.gz > CNS/last/s11.fastq
# runs now:
for file in BAM/*.mapped.sorted.bam; do
  $new_samtools mpileup -uf refs/REF_NC_045512.2.fasta $file | $new_bcftools call -mv -Oz --threads 8 -o CNS/calls.vcf.gz # change to bcftools mpileup??
  $new_bcftools index --threads 8 CNS/calls.vcf.gz
  # cat refs/REF_NC_045512.2.fasta | $new_bcftools consensus CNS/calls.vcf.gz > CNS/`basename $file .mapped.sorted.bam`.fastq # output in fasta format already
  cat refs/REF_NC_045512.2.fasta | $new_bcftools consensus CNS/calls.vcf.gz > CNS/`basename $file .mapped.sorted.bam`.fasta
done

rm CNS/calls.vcf.gz CNS/calls.vcf.gz.csi


# (8) FASTQ TO FASTA
#chmod 755 /data/software/seqtk

## (8) fastq => fasta ## already fasta, no need
#chmod 755 /data/software/seqtk/seqtk
## export PATH=$PATH:/data/software/seqtk/seqtk
#for file in CNS/*.fastq; do
#  seqtk seq -a $file > CNS/`basename $file .fastq`.fasta
#done

for file in CNS/*.fasta; do
  # -i edits file in place
  name=${file/CNS\//} # ${var/find/replace} => remove 'CNS/' prefix
  sed -i "s/>.*/>${name%%.*}/" "$file"
done

### EXTRA:

# align with MAFFT
# dana, before aligning all consensus sequences against the reference, you have to gather all .fasta CNS files into one file, along with the reference. that’s the input.
# preper input file:
cat CNS/*.fasta refs/REF_NC_045512.s.fasta > alignment/all_not_aligned.fasta
#https://towardsdatascience.com/how-to-perform-sequence-alignment-on-2019-ncov-with-mafft-96c1944da8c6
mafft --clustalout alignment/all_notAligned.fasta > alignment/all_aligned.clustalout
mafft alignment/all_notAligned.fasta > alignment/all_aligned.fasta

# TODO: create report!!!
# make sure report files are empty
''>TotalNumReads.txt
''>BAM/TotalMappedReads.txt
for file in BAM/*.bam; do
  if [[ $file == *.mapped*.bam ]]; then
    continue
  fi
  samtools view -c $file >> BAM/TotalNumReads.txt # total num of reads
  samtools view -c -F 260 $file >> BAM/TotalMappedReads.txt
done
