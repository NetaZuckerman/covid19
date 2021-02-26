######### PIPELINE TO PROCESS RAW FASTQ TO CONSENSUS FASTA   ####################
# notes:
# 1. should use bwa mem -p for inteleaved, or two files.
## (3) Map reads to corona virus (REF_NC_045512.2)

old_samtools=/usr/bin/samtools

# index reference - no need right now
#export PATH=$PATH:/Users/netazuck/Documents/software/BWA/bwa-0.7.12
#bwa index data/refs/REF_NC_045512.2.fasta
## interleaved -> if paired end, command should be 'bwa mem -p'! so it is a mistake. try both ways and see difference.
# map reads to reference
#(the command) bwa mem refs/REF_NC_045512.2.fasta fastq/interleaved_trimmed/s79.fastq.gz > SAM/s79.sam
for r1 in fastq/raw/*R1*.fastq.gz; do
  r2=${r1/R1/R2}
  file=${r1/_R1/}
  bwa mem refs/REF_NC_045512.2.fasta $r1 $r2 > SAM/`basename $file .fastq.gz`.sam
done


#  (4) SAM to BAM
#(the command) samtools view -Sb SAM/s11.sam  >  BAM/s11.bam
for file in SAM/*.sam; do
  echo `basename $file`
  $old_samtools view -Sb $file  >  BAM/`basename $file .sam`.bam
done


# (5) KEEP ONLY MAPPED READS
# (the command) samtools view -b -F 260 BAM/s11.bam > BAM/s11.mapped.bam
# neta used new samtools for this command
for file in BAM/*.bam; do
  echo `basename $file`
  samtools view -b -F 260 $file > BAM/`basename $file .bam`.mapped.bam
done

# (6) SORT AND INDEX BAM FILES

# SORT
# (the command) samtools sort BAM/s11.mapped.bam BAM/s11.mapped.sorted
for file in BAM/*.mapped.bam; do
  echo `basename $file`
  samtools sort $file > BAM/`basename $file .mapped.bam`.mapped.sorted
done

# INDEX
# (the command) samtools index BAM/s11.mapped.sorted.bam
for file in BAM/*.mapped.sorted.bam; do
  echo `basename $file`
  $old_samtools index $file
done


# (7) CREATE CONSENSUS SEQUENCE
# https://github.com/samtools/bcftools/wiki/HOWTOs#consensus-calling
# http://samtools.github.io/bcftools/howtos/consensus-sequence.html
# bcftools has to be new - commands not in the old version at all. the old here is the use of samtools mpileup instead of bcftools mpileup command.
# maybe Ns missing are in BCFtools and not samtools command. However, bcftools consensus gets a bed file with locations to add Ns, that is depth from samtools
# but only lines that has depth 0 (3rd col).

for file in BAM/*.mapped.sorted.bam; do
  echo `basen ame $file`
  samtools mpileup -uf refs/REF_NC_045512.2.fasta $file | bcftools call -mv -Oz -o CNS/calls.vcf.gz
  bcftools index CNS/calls.vcf.gz
  cat refs/REF_NC_045512.2.fasta | bcftools consensus CNS/calls.vcf.gz > CNS/`basename $file .mapped.sorted.bam`.fastq
done

# (the commands)
#/data/software/samtools/samtools-1.10_new/samtools-1.10/samtools mpileup -uf refs/REF_NC_045512.2.fasta BAM/last/s11.mapped.sorted.bam | /data/software/bcftools/bcftools-1.10.2_new/bcftools-1.10.2/bcftools call -mv -Oz -o CNS/last/calls.vcf.gz
#/data/software/bcftools/bcftools-1.10.2_new/bcftools-1.10.2/bcftools index CNS/last/calls.vcf.gz
#cat refs/REF_NC_045512.2.fasta | /data/software/bcftools/bcftools-1.10.2_new/bcftools-1.10.2/bcftools consensus CNS/last/calls.vcf.gz > CNS/last/s11.fastq

# (8) FASTQ TO FASTA
#chmod 755 /data/software/seqtk
#export PATH=$PATH:/Users/netazuck/Documents/software/seqtk    #export PATH=$PATH:/data/software/seqtk
# (the command) seqtk seq -a CNS/s11.fastq > CNS/s11.fasta

for file in CNS/*.fastq; do
  echo `basename $file`
  seqtk seq -a $file > CNS/`basename $file .fastq`.fasta
done


### EXTRA:

# align with MAFFT

# dana, before aligning all consensus sequences against the reference, you have to gather all .fasta CNS files into one file, along with the reference. that’s the input.
cat CNS/*.fasta refs/REF_NC_045512.2.fasta > alignment/notAligned.fasta
#https://towardsdatascience.com/how-to-perform-sequence-alignment-on-2019-ncov-with-mafft-96c1944da8c6
mafft --clustalout alignment/notAligned.fasta > alignment/aligned.clustalout
mafft alignment/notAligned.fasta > alignment/aligned.fasta


# total number of reads  (optional - for the report?)
# (the command) samtools view -c BAM/s79.bam
for file in BAM/last/*.bam; do
  echo `basename $file`
  $old_samtools view -c $file >> BAM/last/TotalNumReads.txt
done


# number of reads mapped to corona ref    (optional - for the report?)
#samtools view -c -F 260 BAM/s79.bam
for file in BAM/*.bam; do
  echo `basename $file`
  samtools view -c -F 260 $file >> BAM/TotalMappedReads.txt
done


#calculate breadth and coverage    (that’s the perl code I told you about that i don’t like..)
#chmod 755 /data/projects/Michal/nCoV2019/code/weeSAMv.pl
#/usr/local/bin/perl /data/projects/Michal/nCoV2019/code/weeSAMv.pl -b BAM/last/s11.mapped.sorted.bam -out QC/CoverageStats/last/s11_coverage_stats.txt
for file in BAM/*.mapped.sorted.bam; do
  echo `basename $file`
  perl /data/projects/Michal/nCoV2019/code/weeSAMv.pl -b $file -out QC/CoverageStats/`basename $file .mapped.sorted.bam`.CoverageStats.txt
done
#/usr/local/bin/perl /data/projects/Michal/nCoV2019/code/weeSAMv.pl -b $file -out QC/CoverageStats/`basename $file .mapped.sorted.bam`.CoverageStats.txt  #local computer

###########

# try parallel for loops
cores=3

for r1 in fastq/raw/*R1*.fastq.gz; do
    r2=${r1/R1/R2}
    output=${r1/_R1/}
    output=${output/_paired/}
    output=${output/.gz/}
    bwa mem -v1 refs/REF_NC_045512.fasta "$r1" "$r2" | samtools view -b - > BAM/`basename $output.fastq`.bam &

    background=( $(jobs -p) )
    if (( ${#background[@]} == cores )); then
        wait -n
    fi
done
