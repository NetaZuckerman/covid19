#!/bin/bash

# map
for r1 in fastq/raw/*R1*.fastq.gz; do
  mkdir -p onlyR1/BAM
  output=${r1/_R1/}
  bwa mem -v1 -t10 refs/REF_NC_045512.2.fasta "$r1" | samtools view -b - > onlyR1/BAM/`basename $output .fastq.gz`.bam
done

# keep mapped reads
for file in onlyR1/BAM/*.bam; do
  samtools view -F 260 -h "$file" > onlyR1/BAM/`basename "$file" .bam`.mapped.bam
done

# sort and index bam files
for file in onlyR1/BAM/*.mapped.bam; do
  samtools sort "$file" -o onlyR1/BAM/`basename "$file" .mapped.bam`.mapped.sorted.bam
done
for file in onlyR1/BAM/*.mapped.sorted.bam; do
  samtools index "$file"
done

#depth files
mkdir -p onlyR1/QC/depth
for file in onlyR1/BAM/*.mapped.sorted.bam; do
  samtools depth -a "$file" > onlyR1/QC/depth/`basename "$file" .mapped.sorted.bam`.txt
done

mkdir -p onlyR1/CNS onlyR1/CNS_5
# consenus
for file in onlyR1/BAM/*.mapped.sorted.bam; do
  file_name=`basename "$file" .mapped.sorted.bam`
  file_name=$( echo "$file_name" | cut -d'_' -f 1 )
  samtools mpileup -A "$file" | ivar consensus -m 1 -p onlyR1/CNS/"$file_name"
  samtools mpileup -A "$file" | ivar consensus -m 5 -p onlyR1/CNS_5/"$file_name"
done
rm onlyR1/CNS/*.qual.txt onlyR1/CNS_5/*.qual.txt

# change fasta headers
for file in onlyR1/CNS/*.fa*; do
  name=`basename "$file"`
  sed -i "s/>.*/>${name%%.*}/" "$file"
done
for file in onlyR1/CNS_5/*.fa*; do
  name=`basename "$file"`
  sed -i "s/>.*/>${name%%.*}/" "$file"
done

# align with augur (mafft)
mkdir -p onlyR1/alignment
cat CNS_5/*.fa* > onlyR1/alignment/all_not_aligned.fasta
augur align \
--sequences onlyR1/alignment/all_not_aligned.fasta \
--reference-sequence refs/REF_NC_045512.2.fasta \
--output onlyR1/alignment/all_aligned.fasta

mkdir -p onlyR1/results
python /data/projects/Dana/scripts/covid19/MutTable.py onlyR1/alignment/all_aligned.fasta onlyR1/results/muttable.csv
python /data/projects/Dana/scripts/covid19/variants.py onlyR1/alignment/all_aligned.fasta onlyR1/results/variants.csv

report=onlyR1/QC/report.txt
echo -e "sample\tmapped%\tmappedreads\ttotreads\tcovbases\tcoverage%\tcoverageCNS_5%\tmeandepth\tmaxdepth\tmindepth" > "$report"
for file in onlyR1/BAM/*.mapped.sorted.bam; do
  sample_name=`basename "$file" .mapped.sorted.bam`
  sample_name=$( echo "$sample_name" | cut -d'_' -f 1 )

  original_bam=${file/.mapped.sorted.bam/.bam}
  tot_reads=$(samtools view -c "$original_bam")
  coverage_stats=( $(samtools coverage -H "$file" | cut -f4,5,6) )
  breadth_cns5=$(cut -f3 onlyR1/QC/depth/`basename "$file" .mapped.sorted.bam`.txt | awk '$1>5{c++} END{print c+0}')
  genome_size=$(cat onlyR1/QC/depth/`basename $file .mapped.sorted.bam`.txt | wc -l)
  coverage_cns5=$(echo "$breadth_cns5/$genome_size"*100 | bc -l)
  mapped_num=${coverage_stats[0]}
  percentage_mapped=$(awk -v m="$mapped_num" -v t="$tot_reads" 'BEGIN {print (m/t)*100}')
  depth=$(awk '{if($3==0){next}; if(min==""){min=max=$3}; if($3>max) {max=$3}; if($3< min) {min=$3}; total+=$3; count+=1} END {print total/count"\t"max"\t"min}' onlyR1/QC/depth/`basename $file .mapped.sorted.bam`.txt)
  echo -e "${sample_name}\t${percentage_mapped}\t${mapped_num}\t${tot_reads}\t${coverage_stats[1]}\t${coverage_stats[2]}\t${coverage_cns5}\t${depths}" >> "$report"
done
#############################################################################################


