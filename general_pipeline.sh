#!/bin/bash

trim_flag=false
refseq=refs/REF_NC_045512.2.fasta
new_samtools=/data/software/samtools/samtools-1.10_new/samtools-1.10/samtools
new_bcftools=/data/software/bcftools/bcftools-1.10.2_new/bcftools-1.10.2/bcftools
# parse input wiht flags
while (( "$#" )); do
  case "$1" in
    -t|--trimmed_fq)
      trim_flag=true
      shift
      ;;
    -*|--*=) # unsupported flags
      echo "Error: Unsupported flag $1" >&2
      exit 1
      ;;
  esac
done

# (3) Map reads to corona virus (REF_NC_045512.2)
# index reference
if [ ! -f "$refseq" ]; then # TODO: make reference changeable!
  echo "$refseq" does not exist. Please provide the reference sequence and try again.
  exit 1
fi
bwa index "$refseq"

# map reads to reference. straight to bam format! assuming PE for now # TODO
if $trim_flag # data is trimmed
then
  for r1 in fastq/trimmed/*R1*_paired.fastq.gz; do
    r2=${r1/R1/R2} # ${var/find/replace}
    output=${r1/_R1/}
    bwa mem -v1 -t4 refs/REF_NC_045512.2.fasta "$r1" "$r2" | samtools view -@ 8 -Sb - > BAM/`basename $output _paired.fastq.gz`.bam
  done
else # data is raw
  for r1 in fastq/raw/*R1*.fastq.gz; do
    r2=${r1/R1/R2}
    output=${r1/_R1/}
    bwa mem -v1 -t4 refs/REF_NC_045512.2.fasta "$r1" "$r2" | samtools view -@ 8 -Sb - > BAM/`basename $output .fastq.gz`.bam
  done
fi

# (5) keep only mapped reads
for file in BAM/*.bam; do
   $new_samtools view -@ 8 -b -F 260 $file > BAM/`basename $file .bam`.mapped.bam
done

# (6) sort and index bam files
for file in BAM/*.mapped.bam; do
  sorted=${file/.mapped.bam/mapped.sorted.bam}
  samtools sort -@ 8 $file BAM/`basename $file .mapped.bam`.mapped.sorted
  samtools index "$sorted"
done

# (7) create concensus sequence
for file in BAM/*.mapped.sorted.bam; do
  $new_samtools mpileup -uf refs/REF_NC_045512.2.fasta $file | $new_bcftools call -mv -Oz --threads 8 -o CNS/calls.vcf.gz # change to bcftools mpileup??
  $new_bcftools index --threads 8 CNS/calls.vcf.gz
  # cat refs/REF_NC_045512.2.fasta | $new_bcftools consensus CNS/calls.vcf.gz > CNS/`basename $file .mapped.sorted.bam`.fasta
  $new_bcftools consensus -f refs/REF_NC_045512.2.fasta CNS/calls.vcf.gz > CNS/`basename $file .mapped.sorted.bam`.fasta
done

rm CNS/calls.vcf.gz CNS/calls.vcf.gz.csi

# change fasta header from ref to sample name
for file in CNS/*.fasta; do
  name=${file/CNS\//} # ${var/find/replace} => remove 'CNS/' prefix
  sed -i "s/>.*/>${name%%.*}/" "$file"
done

### EXTRA:
# align with MAFFT
# https://towardsdatascience.com/how-to-perform-sequence-alignment-on-2019-ncov-with-mafft-96c1944da8c6
cat CNS/*.fasta refs/REF_NC_045512.2.fasta > alignment/all_not_aligned.fasta
mafft --clustalout alignment/all_not_aligned.fasta > alignment/all_aligned.clustalout
mafft alignment/all_not_aligned.fasta > alignment/all_aligned.fasta

mkdir -p results/
report=results/report.txt
# samtools coverage headers: 1#rname  2startpos  3endpos    4numreads  5covbases  6coverage  7meandepth  8meanbaseq  9meanmapq
echo -e "sample\tmapped%\tmappedreads\ttotreads\tcovbases\tcoverage%\tmeandepth\tmaxdepth\tmindepth" > "$report"
for file in BAM/*.mapped.sorted.bam; do
  sample_name=${file/BAM\//}
  original_bam=${file/.mapped.sorted.bam/.bam} #{var/find/replace}
  tot_reads=$(samtools view -c "$original_bam") # do not use -@n, no output inside variable
  line=( $($new_samtools coverage -H "$file" | cut -f4,5,6) )
  mapped_num=${line[0]}
  percentage_mapped=$(awk -v m="$mapped_num" -v t="$tot_reads" 'BEGIN {print (m/t)*100}')
  $new_samtools depth "$file" > depth.txt
  max_depth=$(awk 'BEGIN{max=0} {if ($3>max) max=$3} END {print max}' depth.txt)
  min_depth=$(awk 'BEGIN {min=$max_depth} {if ($3<min) min=$3} END {print min}' depth.txt)
  mean_depth=$(awk '{c++;s+=$3}END{print s/c}' depth.txt)
  # echo to report
  line="${sample_name}\t${percentage_mapped}\t${mapped_num}\t${tot_reads}\t${line[1]}\t${line[2]}\t${mean_depth}\t${max_depth}\t${min_depth}"
  echo -e "$line" >> "$report"
done

rm depth.txt

# for parallel: TODO
