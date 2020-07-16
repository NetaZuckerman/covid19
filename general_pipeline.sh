#!/bin/bash
# TODO: for each sample create a pipe for better runtime

# (3) Map reads to corona virus (REF_NC_045512.2)
# index reference
# export?
bwa index refs/REF_NC_045512.2.fasta
# map reads to reference -> PE
# TODO: May be ran on trimmed or on raw. get wanter option form user!!
for r1 in fastq/trimmed/*R1*_paired.fastq.gz; do
  r2=${r1/R1/R2} # ${var/find/replace}
  output=${r1/_R1/}
  # old:   bwa mem -v1 -t4 refs/REF_NC_045512.2.fasta $r1 $r2 > SAM/`basename $output _paired.fastq.gz`.sam
  bwa mem -v1 -t4 refs/REF_NC_045512.2.fasta "$r1" "$r2" | samtools view -@ 8 -Sb - > BAM/`basename $output _paired.fastq.gz`.bam
done

## samtools <command> -@: number of threads (default 1)
## (4) sam to bam
#for file in SAM/*.sam; do
#  samtools view -@ 8 -Sb $file > BAM/`basename $file .sam`.bam
#done

new_samtools=/data/software/samtools/samtools-1.10_new/samtools-1.10/samtools
new_bcftools=/data/software/bcftools/bcftools-1.10.2_new/bcftools-1.10.2/bcftools
# (5) keep only mapped reads
for file in BAM/*.bam; do
   $new_samtools view -@ 8 -b -F 260 $file > BAM/`basename $file .bam`.mapped.bam
done

# (6) sort and index bam files
# sort
for file in BAM/*.mapped.bam; do
  sorted=${file/.mapped.bam/mapped.sorted.bam}
  samtools sort -@ 8 $file BAM/`basename $file .mapped.bam`.mapped.sorted
  samtools index "$sorted"
done

## index
#for file in BAM/*.mapped.sorted.bam; do
#  samtools index $file
#done

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
  # cat refs/REF_NC_045512.2.fasta | $new_bcftools consensus CNS/calls.vcf.gz > CNS/`basename $file .mapped.sorted.bam`.fasta
  $new_bcftools consensus -f refs/REF_NC_045512.2.fasta CNS/calls.vcf.gz > CNS/`basename $file .mapped.sorted.bam`.fasta
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
cat CNS/*.fasta refs/REF_NC_045512.2.fasta > alignment/all_not_aligned.fasta
#https://towardsdatascience.com/how-to-perform-sequence-alignment-on-2019-ncov-with-mafft-96c1944da8c6
mafft --clustalout alignment/all_not_aligned.fasta > alignment/all_aligned.clustalout
mafft alignment/all_not_aligned.fasta > alignment/all_aligned.fasta

# TODO: create report!!!
# make sure report files are empty
"">stats.txt

#for file in BAM/*.bam; do
#  if [[ $file == *.mapped*.bam ]]; then
#    continue
#    #TODO: coverage - move here?
#  fi
#  tot_reads=$(samtools view -c $file)
#  mapped_reads=$(samtools view -c -F 260 $file)
#  echo -e "$file\t$tot_reads\t$mapped_reads\t$((mapped_reads / tot_reads))" >> stats.txt
#done
# coverage headers: 1#rname  2startpos  3endpos    4numreads  5covbases  6coverage  7meandepth  8meanbaseq  9meanmapq
echo -e "sample\tmappedreads\tcovbases\tcoverage\%\tmeandepth\tmaxdepth\tmindepth\ttot_reads" > report.txt
for file in BAM/*.mapped.sorted.bam; do
  sample_name=${file/BAM\//}
  original_bam=${file/.mapped.sorted.bam/.bam} #{var/find/replae}
  tot_reads=$(samtools view -c "$original_bam")
  line=$($new_samtools coverage -H "$file" | cut -f4,5,6)
  $new_samtools depth -a "$file" > depth.txt
  max_depth=$(awk 'BEGIN{max=0} {if ($3>max) max=$3} END {print max}' depth.txt) # calculate max depth
  min_depth=$(awk 'BEGIN {min=0} {if ($3<min) min=$3} END {print min}' depth.txt) # calculate min depth
  mean_depth=$(awk '{c++;s+=$3}END{print s/c}' depth.txt)
  line="${sample_name}\t${line}\t${mean_depth}\t${max_depth}\t${min_depth}\t${tot_reads}"
  echo -e "$line" >> report.txt
done

