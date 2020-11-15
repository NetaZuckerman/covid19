#!/bin/bash
###############################
# Dana Bar-Ilan
# 21.07.20
#
# mapping fq files to ref-seq <bwa> -> keep mapped reads, sort, index <samtools> -> create consensus sequence <samtools & bcftools>
# -> align consensuses to ref-seq <mafft>
# + produce report
# version: samtools-1.10
#
trap "kill 0" EXIT

function initialize_globals() {
  dirs_flag=false
# refseq=refs/REF_NC_045512.2.fasta # default refseq
  new_bcftools=/data/software/bcftools/bcftools-1.10.2_new/bcftools-1.10.2/bcftools
  threads=32
  input_path=""
}

# parse input wiht flags
function usage() {
    cat <<EOF
Usage: $0 [options]

required: [-h | -d | -i AND -r]
-h| --help                      print this usage message and exit. Ignore the rest
-d|--create_dirs                create all project's directories in current working directory.
-i              [fastq.gz/path] input path to fastq.gz files location.
-r|--refseq     [refseq/path/]  user defined reference. required: refseq/path/ - path to reference fasta file.

optional:
--threads       [int]           number of threads. default: 32
EOF
exit 0
}

function get_user_input() {
  while (( "$#" )); do
    case "$1" in
      -d|--create_dirs) # just create all directories
        dirs_flag=true
        shift
        ;;
      -r|--refseq) # user provided refseq
        shift
        refseq="$1"
        shift
        ;;
      --threads)
        shift
        threads="$1"
        shift
        ;;
      -i) # path to fastq.gz files location
        shift
        input_path="$1"
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
}

function check_flags() {
  if $dirs_flag; then
    mkdir -p fastq/{raw,trimmed} QC/fastqc refs BAM CNS CNS_5 alignment
    echo "Created project directories. Please download your data to fastq/raw and/or fastq/trimmed, and your reference sequence to refs/. "
    exit 0
  fi

  if [ -z "$input_path" ]; then
    echo "Please provide -i <input_path> to your fastq.gz files location."
    usage
    exit 1
  fi

  if [ -z "$refseq" ]; then
    mkdir -p refs
    echo please provide reference sequence --refseq!
    exit 1
  fi
}

# map reads to reference. straight to bam format! assuming PE for now # TODO
function map_to_ref() {
  # index reference
  bwa index "$refseq"

  if [ -d BAM/ ]; then # not first run, rm BAM files to avoid mixups
    rm BAM/* 2> /dev/null # if file not found it's ok. no need to see on screen
  fi
  if [ -d CNS/ ]; then
    rm CNS/* 2> /dev/null
  fi
  if [ -d CNS_5/ ]; then
    rm CNS_5/ 2> /dev/null
  fi
  mkdir -p BAM CNS alignment results Trees

  for r1 in "$input_path"*R1*.fastq*; do
    if [[ $r1 == *Undetermined*.fastq* || $r1 == *unpaired*.fastq* ]]; then # ignore undetermined
        continue
    fi
    r2=${r1/R1/R2} # ${var/find/replace}
    output=${r1/_R1/}
    output=${output/_paired/}
    output=${output/.gz/}
    bwa mem -v1 -t"$threads" "$refseq" "$r1" "$r2" | samtools view -@ "$threads" -b - > BAM/`basename $output .fastq`.bam
  done
}

function keep_mapped_reads() {
  for file in BAM/*.bam; do
     samtools view -@ "$threads" -b -F 260 $file > BAM/`basename $file .bam`.mapped.bam
  done
}

function sort_index_bam() {
  for file in BAM/*.mapped.bam; do
    samtools sort -@ "$threads" "$file" -o BAM/`basename "$file" .mapped.bam`.mapped.sorted.bam
    samtools index -@ "$threads" BAM/`basename "$file" .mapped.bam`.mapped.sorted.bam
  done
}


function depth() {
    mkdir -p QC/depth/
    # create depth files from each mapped sorted and indexed bam file. -a: include all sample's positions including 0 depth.
    for file in BAM/*.mapped.sorted.bam; do
      samtools depth -a "$file" > QC/depth/`basename "$file" .mapped.sorted.bam`.txt
    done
    # create bed file with 0 depth regions from each depth file
    # multi process TODO
    # in the future maybe pipe the depth into the rest of the command?
#    for file in QC/depth/*.txt; do
#      cut -f 3 "$file" | sort | uniq | while read X; do awk -v X="$X" '($3==X && $3==0) { printf("%s\t%d\t%d\n",$1,int($2)-1,int($2)+1);}' "$file" | sort -k1,1 -k2,2n | bedtools merge -i - | sed "s/\$/\t${X}/" ; done > QC/depth/`basename "$file" .txt`.bed
#    done
}

function consensus() {
  for file in BAM/*.mapped.sorted.bam; do
    sample_name=`basename $file .mapped.sorted.bam`
    #bed_file=QC/depth/`basename $file .mapped.sorted.bam`.bed
    bcftools mpileup -f "$refseq" "$file" | bcftools call -mv -Oz -o CNS/"$sample_name"_calls.vcf.gz
    bcftools index CNS/"$sample_name"_calls.vcf.gz
#    bcftools consensus -f "$refseq" -m "$bed_file" CNS/"$sample_name"_calls.vcf.gz > CNS/"$sample_name".fasta
    bcftools consensus -f "$refseq" CNS/"$sample_name"_calls.vcf.gz > CNS/"$sample_name".fasta
    # mask low depth:
    # mask 0 depth:
    python /home/dana/covid19/mask_fasta.py CNS/"$sample_name".fasta QC/depth/`basename $file .mapped.sorted.bam`.txt -n 1
    # mask 5 depth: every position under 5 depth is N.
    mkdir -p CNS_5/
    python /home/dana/covid19/mask_fasta.py CNS/"$sample_name".fasta QC/depth/`basename $file .mapped.sorted.bam`.txt -n 5 -o CNS_5/"$sample_name".fasta
    rm CNS/"$sample_name"_calls.vcf.gz CNS/"$sample_name"_calls.vcf.gz.csi
  done

#  for file in BAM/*.mapped.sorted.bam; do
#    samp_name=`basename $file .mapped.sorted.bam`
#    bed_file=QC/depth/"$sample_name".bed  # matching bed file
#    $new_bcftools index --threads "$threads" CNS/"$samp_name"_calls.vcf.gz
#    $new_bcftools consensus -f "$refseq" -m "$bed_file" CNS/"$samp_name"_calls.vcf.gz > CNS/`basename $file .mapped.sorted.bam`.fasta
#    rm CNS/"$samp_name"_calls.vcf.gz CNS/"$samp_name"_calls.vcf.gz.csi
#  done
}

# change fasta header from ref to sample name
function change_fasta_header() {
  for file in CNS/*.fasta; do
    # change header to sample name:
    name=${file/CNS\//} # ${var/find/replace} => remove 'CNS/' prefix
    sed -i "s/>.*/>${name%%.*}/" "$file"
  done
}

function mafft_alignment() {
  # align with MAFFT
  # https://towardsdatascience.com/how-to-perform-sequence-alignment-on-2019-ncov-with-mafft-96c1944da8c6
  cat CNS/*.fasta "$refseq" > alignment/all_not_aligned.fasta
  mafft --clustalout alignment/all_not_aligned.fasta > alignment/all_aligned.clustalout
  mafft alignment/all_not_aligned.fasta > alignment/all_aligned.fasta
}

function results_report() {
  report=QC/report.txt
  # samtools coverage headers: 1#rname  2startpos  3endpos    4numreads  5covbases  6coverage  7meandepth  8meanbaseq  9meanmapq
  echo -e "sample\tmapped%\tmappedreads\ttotreads\tcovbases\tcoverage%\tmeandepth\tmaxdepth\tmindepth" > "$report"
  for file in BAM/*.mapped.sorted.bam; do
    sample_name=`basename $file .mapped.sorted.bam`
    original_bam=${file/.mapped.sorted.bam/.bam} #{var/find/replace}
    tot_reads=$(samtools view -c "$original_bam") # do not use -@n when capturing output in variable
    coverage_stats=( $(samtools coverage -H "$file" | cut -f4,5,6) ) # number of mapped reads, covered bases, coverage
    mapped_num=${coverage_stats[0]}
    percentage_mapped=$(awk -v m="$mapped_num" -v t="$tot_reads" 'BEGIN {print (m/t)*100}')
    depths=$(awk '{if($3==0){next}; if(min==""){min=max=$3}; if($3>max) {max=$3}; if($3< min) {min=$3}; total+=$3; count+=1} END {print total/count"\t"max"\t"min}' QC/depth/"$sample_name".txt)
    echo -e "${sample_name}\t${percentage_mapped}\t${mapped_num}\t${tot_reads}\t${coverage_stats[1]}\t${coverage_stats[2]}\t${depths}" >> "$report"
  done
}

########################### MAIN ###############################
# trap ctrl-c to end in the same directory as started even if user ended the program

# call all functions
# user input:
initialize_globals
get_user_input "$@"
check_flags
# start workflow:
map_to_ref
keep_mapped_reads
sort_index_bam
depth
consensus
change_fasta_header
mafft_alignment
results_report
wait
echo "pipeline finished! (:"

