#!/bin/bash
###############################
# Dana Bar-Ilan
# 21.07.20
# The pipeline:
# mapping fq files to ref-seq <bwa> -> keep mapped reads, sort, index <samtools> -> create consensus sequence <samtools & bcftools>
# -> align consensuses to ref-seq <mafft>
# + produce report


trim_flag=false
dirs_flag=false
refseq=refs/REF_NC_045512.2.fasta # default refseq
new_samtools=/data/software/samtools/samtools-1.10_new/samtools-1.10/samtools
new_bcftools=/data/software/bcftools/bcftools-1.10.2_new/bcftools-1.10.2/bcftools
# parse input wiht flags
function usage() {
    cat <<EOF
Usage: $0 [options]

-h| --help                      Print this usage message and exit. Ignore the rest
-t|--trimmed_fq                 Run the pipeline with trimmed fsatq data (instead of raw).
-d|--create_dirs                Create all project's directories in current working directory.
-r|--refseq      refseq/path/   User defined reference. Required: refseq/path/ - path to reference fasta file.
                                default: refs/REF_NC_045512.2.fasta
EOF
exit 0
}

while (( "$#" )); do
  case "$1" in
    -t|--trimmed_fq)
      trim_flag=true
      shift
      ;;
    -d|--create_dirs) # just create all directories
    dirs_flag=true
    shift
    ;;
  -r|--refseq) # user provided refseq
    shift
    refseq="$1"
    shift
    ;;
  -h|--help)
    usage
    shift
    ;;
    -*|--*=) # unsupported flags
      echo "Error: Unsupported flag $1" >&2
      exit 1
      ;;
  esac
done

if $dirs_flag; then
  mkdir -p fastq/{raw,trimmed} QC/{fastqc} refs BAM CNS alignment Trees results
  echo "Created project directories. Please download your data to fastq/raw and/or fastq/trimmed, and your reference sequence to refs/. "
  exit 0
fi

# (3) Map reads to corona virus (REF_NC_045512.2)
# index reference
if [ ! -f "$refseq" ]; then
  mkdir -p refs
  echo "$refseq" does not exist. Please provide the reference sequence and try again.
  exit 1
fi

bwa index "$refseq"

mkdir -p BAM CNS alignment results Trees
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

### RESULTS REPORT
report=results/report.txt
# samtools coverage headers: 1#rname  2startpos  3endpos    4numreads  5covbases  6coverage  7meandepth  8meanbaseq  9meanmapq
echo -e "sample\tmapped%\tmappedreads\ttotreads\tcovbases\tcoverage%\tmeandepth\tmaxdepth\tmindepth" > "$report"
for file in BAM/*.mapped.sorted.bam; do
  sample_name=${file/BAM\//}
  original_bam=${file/.mapped.sorted.bam/.bam} #{var/find/replace}
  tot_reads=$(samtools view -c "$original_bam") # do not use -@n when capturing output in variable
  line=( $($new_samtools coverage -H "$file" | cut -f4,5,6) )
  mapped_num=${line[0]}
  percentage_mapped=$(awk -v m="$mapped_num" -v t="$tot_reads" 'BEGIN {print (m/t)*100}')
  $new_samtools depth "$file" > depth.txt
  depths=$(awk '{if(min==""){min=max=$3}; if($3>max) {max=$3}; if($3< min) {min=$3}; total+=$3; count+=1} END {print total/count"\t"max"\t"min}' depth.txt)
  echo -e "${sample_name}\t${percentage_mapped}\t${mapped_num}\t${tot_reads}\t${line[1]}\t${line[2]}\t${depths}" >> "$report"
done

rm depth.txt

# for parallel: TODO manage up to N parallel runs to save time
