#!/bin/bash
###############################
# Dana Bar-Ilan
# 21.07.20
#
# in general:
# mapping fq files to ref-seq <bwa> -> keep mapped reads, sort, index <samtools> -> create consensus sequence <samtools & bcftools>
# -> align consensuses to ref-seq <mafft>
# + produce report

# requirements: samtools v1.10, bcftools v1.9, ivar v. 1.2.2, mafft v7.215
trap "kill 0" EXIT
eval "$(conda shell.bash hook)"
conda activate CoronaPipeline


if [ -d BAM/ ]; then # not first run, rm existing files to avoid mix ups
  rm BAM/* 2> /dev/null #  if file not found it's ok. no need to print error to screen.
fi
if [ -d CNS/ ]; then
  rm CNS/* 2> /dev/null
fi
if [ -d CNS_5/ ]; then
  rm CNS_5/* 2> /dev/null
fi
mkdir -p BAM CNS alignment results  # -p: only if directory does not exist. else continue without errors.






path=`dirname "${0}"`

function initialize_globals() {
  dirs_flag=false
  threads=20
  num_processes=1
  input_path=""
  single_end=false
  newNextclade=false
  nextclade25=false
  nextalign25=false
  spike=false
  SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

}

# parse input with flags
function usage() {
    cat <<EOF
Usage: $0 [options]

required: [-h | -d | -i AND -r]
-h| --help                      print this usage message and exit. Ignore the rest
-d|--create_dirs                create all project's directories in current working directory
-i              [fastq.gz/path] input path to fastq.gz files location
-r|--refseq     [refseq/path/]  user defined reference. required: refseq/path/ - path to reference fasta file

optional:
--threads       [int]           number of threads for each sample. default: 32
-p|--processes  [int]      number of processes (samples) to run in parallel. default: 1
-s|--single-end                 single end sequencing
-n| --newNextclade              use nextclade version 1.3.0
--spike                         spike sequencing only

EOF
exit 0
}

function get_user_input() {
  while (( "$#" )); do
    case "$1" in
      -d|--create_dirs)
        dirs_flag=true
        shift
        ;;
      -r|--refseq)
        shift
        refseq="$1"
        shift
        ;;
      -p|--processes)
        shift
        num_processes="$1"
        shift
        ;;
      --threads)
        shift
        threads="$1"
        shift
        ;;
      -i)
        shift
        input_path="$1"
        shift
        ;;
      -h|--help)
        usage
        shift
        ;;
      -s|--single-end)
        single_end=true
        shift
        ;;
      -n| --newNextclade)
        newNextclade=true
        shift
        ;;
      --nextclade2.5)
        nextclade25=true
        shift
        ;;     
      --nextalign2.5)
        nextalign25=true
        shift
        ;;     
      -q| --quasispecies)
        q=true
        shift
        ;;       
      --spike)
        spike=true
        shift
        ;;
      --noRecombinants)
        noRecombinants=true
        shift
        ;;
      -*|--*=)
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
    exit 1
  fi

  if [ -z "$input_path" ]; then
    echo "Please provide -i <input_path> to your fastq.gz files location." >&2
    usage
    exit 1
  fi

  if [ -z "$refseq" ]; then
    mkdir -p refs
    echo please provide reference sequence --refseq! >&2
    exit 1
  fi

  if [[ "$input_path" != */ ]]; then
    input_path="$input_path"/
  fi
}


function single_end_mapping() {
  bwa index "$refseq"
  if [ -d BAM/ ]; then # not first run, rm BAM files to avoid mixups
    rm BAM/* 2> /dev/null # if file not found it's ok. no need to see on screen
  fi
  if [ -d CNS/ ]; then
    rm CNS/* 2> /dev/null
  fi
  if [ -d CNS_5/ ]; then
    rm CNS_5/* 2> /dev/null
  fi
  mkdir -p BAM CNS alignment results

  for file in "$input_path"*.fastq*; do
    if [[ $file == *Undetermined*.fastq* || $file == *unpaired*.fastq*  || $file == *singletons* ]]; then # ignore undetermined and singletons
        continue
    fi

    if [[ $file == *.gz ]]; then
      output=`basename .fastq.gz`
    else
      output=`basename "$file" .fastq`
    fi
    # replace to short name:
    if [[ $output == Sh_* ]]; then
      output=${output/Sh_/}
    fi

    output=BAM/$output.bam

    bwa mem -v1 -t"$threads" "$refseq" "$file" | samtools view -@ "$threads" -b - > $output
  done
}

function bwa_mem(){
    if [[ $r1 == *Undetermined*.fastq* || $r1 == *unpaired*.fastq*  || $r1 == *singletons* ]]; then # ignore undetermined and singletons
        continue
    fi
    r2=${r1/R1/R2}
    output=${r1/_R1/}
    output=${output/_paired/}
    if [[ $r1 == *.gz ]]; then
      output=`basename "$output" .fastq.gz`
    else
      output=`basename "$output" .fastq`
    fi
    # replace to short name:
    if [[ $output == Sh_* ]]; then
      output=${output/Sh_/}
    fi
    output=$( echo "$output" | cut -d'_' -f 1 )  # SHORT NAME
    output=BAM/$output.bam
    bwa mem -v1 -t"$threads" "$refseq" "$r1" "$r2" | samtools view -@ "$threads" -b - > $output
}
# paired end mapping: map reads to reference. straight to bam format.
function map_to_ref() {
  # index reference
  bwa index "$refseq"


  for r1 in "$input_path"*R1*.fastq*; do
    ((i=i%num_processes)); ((i++==0)) && wait
    bwa_mem “$r1” &
  done
}

function keep_mapped_reads() {
  for file in BAM/*.bam; do
    ((i=i%num_processes)); ((i++==0)) && wait
    samtools view -@ "$threads" -b -F 260 $file > BAM/`basename $file .bam`.mapped.bam &
  done
}

function sort_bam(){
  for file in BAM/*.mapped.bam; do
    ((i=i%num_processes)); ((i++==0)) && wait
    samtools sort -@ "$threads" "$file" -o BAM/`basename "$file" .mapped.bam`.mapped.sorted.bam &
  done
}

function index_bam() {
  for file in BAM/*.mapped.bam; do
    ((i=i%num_processes)); ((i++==0)) && wait
    samtools index -@ "$threads" BAM/`basename "$file" .mapped.bam`.mapped.sorted.bam &
  done
}

function pileup(){
  for file in BAM/*.mapped.sorted.bam; do
    ((i=i%num_processes)); ((i++==0)) && wait
    samtools mpileup "$file" -f "$refseq" > pileup/`basename "$file" .mapped.sorted.bam`.pileup &
  done
}

function depth() {
    mkdir -p QC/depth/
    # create depth files from each mapped sorted and indexed bam file. -a: include all sample's positions including 0 depth.
    for file in BAM/*.mapped.sorted.bam; do
      ((i=i%num_processes)); ((i++==0)) && wait
      samtools depth -a "$file" > QC/depth/`basename "$file" .mapped.sorted.bam`.txt &
    done
}


function consensus() {
  mkdir -p CNS CNS_5
  for file in BAM/*.mapped.sorted.bam; do
    # ivar instead of bcftools:
    ((i=i%num_processes)); ((i++==0)) && wait
    file_name=`basename "$file" .mapped.sorted.bam`
    if [ "$single_end" == FALSE ]; then  # if paired end leave quality threshold 20
      # CNS1
      samtools mpileup -A "$file" | ivar consensus -t 0.6 -m 1 -p CNS/"$file_name" &
      # CNS5
      samtools mpileup -A "$file" | ivar consensus -t 0.6 -m 5 -p CNS_5/"$file_name" &
    else  # of single end allow low quality bases (MinIon)
      # CNS1
      samtools mpileup -A "$file" | ivar consensus -t 0.6 -m 1 -p CNS/"$file_name" -q 10 &
      # CNS5
      samtools mpileup -A "$file" | ivar consensus -t 0.6 -m 5 -p CNS_5/"$file_name" -q 10 &
    fi
  done
  wait
  # remove qual files:
  rm CNS/*.qual.txt CNS_5/*.qual.txt
}

# change fasta header from ref to sample name - not needed when using ivar
function change_fasta_header() {
  for file in CNS/*.fa*; do
    # change header to sample name:
    name=`basename $file`
    sed -i "s/>.*/>${name%%.*}/" "$file"
  done

  for file in CNS_5/*.fa*; do
    # change header to sample name:
    name=${file/CNS_5\//} # ${var/find/replace} => remove 'CNS/' prefix
    sed -i "s/>.*/>${name%%.*}/" "$file"
  done
}

function mafft_alignment() {
  # align with MAFFT
  # https://towardsdatascience.com/how-to-perform-sequence-alignment-on-2019-ncov-with-mafft-96c1944da8c6
  cat CNS_5/*.fa* > alignment/all_not_aligned.fasta

  conda activate nextstrain
  if [ "$nextalign25" == true ]; then
    nextalign run -r "$refseq" alignment/all_not_aligned.fasta -o alignment/all_aligned.fasta --output-insertions alignment/insertions.csv
  else 
    nextalign -i alignment/all_not_aligned.fasta -r "$refseq" --output-fasta alignment/all_aligned.fasta --output-insertions alignment/insertions.csv 
  fi
  conda deactivate
  
  

}


function muttable() {

    if [ "$spike" == true ]; then
       python "$path"/variants_spike.py alignment/all_aligned.fasta results/variants.csv "$path"/mutationsTable.xlsx QC/report.txt
    else
    # run pangolin
    conda deactivate
    
    echo "Run Pangolin" 1>&3
    conda activate pangolin
    pangolin alignment/all_not_aligned.fasta --outfile results/pangolinClades.csv
    conda deactivate

  echo "Run Nextclade" 1>&3
    conda activate nextstrain
    if [ "$newNextclade" == true ]; then
      nextclade --input-fasta alignment/all_not_aligned.fasta  --input-dataset $SCRIPT_DIR/nextclade --output-dir nextclade/ --output-tsv results/nextclade.tsv
    elif [ "$nextclade25" == true ]; then
      nextclade run --input-dataset $SCRIPT_DIR/nextclade --output-tsv results/nextclade.tsv alignment/all_not_aligned.fasta
    else 
        nextclade -i alignment/all_not_aligned.fasta -t results/nextclade.tsv
    fi
    conda deactivate

    conda activate CoronaPipeline
    echo "Variants analysis" 1>&3

          python "$path"/MutTable.py alignment/all_aligned.fasta results/nuc_muttable.xlsx  "$path"/mutationsTable.xlsx
          #python "$path"/translated_table.py alignment/all_aligned.fasta results/AA_muttable.xlsx "$path"/regions.csv "$path"/mutationsTable.xlsx
          python "$path"/variants.py alignment/all_aligned.fasta results/variants.csv results/pangolinClades.csv results/nextclade.tsv "$path"/mutationsTable.xlsx "$noRecombinants" QC/report.txt
          echo "Extra mutations report" 1>&3
          python "$path"/translate_extras.py "$path"/covid19_regions.csv results/non_variant_mutations.csv "$refseq" alignment/all_aligned.fasta


    fi
}

function over_50() {
  # create another multifasta with only sequences of <50% N's.
  python "$path"/filterNPercent.py alignment/all_not_aligned.fasta alignment/all_not_aligned_over50.fasta
  python "$path"/filterNPercent.py alignment/all_aligned.fasta alignment/all_aligned_over50.fasta
}

function result() {

  sample_name=`basename $file .mapped.sorted.bam`
  if [[ $sample_name == Sh_* ]]; then
    sample_name=${sample_name/Sh_/}
  fi

  if [ "$single_end" == false ]; then
    sample_name=$( echo "$sample_name" | cut -d'_' -f 1 ) # SHORT NAME
  fi

  original_bam=${file/.mapped.sorted.bam/.bam} # {var/find/replace}
  tot_reads=$(samtools view -c "$original_bam") # do not use -@n when capturing output in variable
  coverage_stats=( $(samtools coverage -H "$file" | cut -f4,5,6) ) # number of mapped reads, covered bases, coverage
  breadth_cns5=$(cut -f3 QC/depth/`basename $file .mapped.sorted.bam`.txt | awk '$1>5{c++} END{print c+0}')
  genome_size=$(cat QC/depth/`basename $file .mapped.sorted.bam`.txt | wc -l)
  coverage_cns5=$(echo "$breadth_cns5/$genome_size"*100 | bc -l)
  mapped_num=${coverage_stats[0]}
  percentage_mapped=$(awk -v m="$mapped_num" -v t="$tot_reads" 'BEGIN {print (m/t)*100}')
  depths=$(awk '{if($3==0){next}; if(min==""){min=max=$3}; if($3>max) {max=$3}; if($3< min) {min=$3}; total+=$3; count+=1} END {print total/count"\t"max"\t"min}' QC/depth/`basename $file .mapped.sorted.bam`.txt)
  echo -e "${sample_name}\t${percentage_mapped}\t${mapped_num}\t${tot_reads}\t${coverage_stats[1]}\t${coverage_stats[2]}\t${coverage_cns5}\t${depths}" >> "$report"
}


function results_report() {
  report=QC/report.txt
  # samtools coverage headers: 1#rname  2startpos  3endpos  4numreads  5covbases  6coverage  7meandepth  8meanbaseq  9meanmapq
  echo -e "sample\tmapped%\tmappedreads\ttotreads\tcovbases\tcoverage%\tcoverageCNS_5%\tmeandepth\tmaxdepth\tmindepth" > "$report"
  for file in BAM/*.mapped.sorted.bam; do
    ((i=i%num_processes)); ((i++==0)) && wait
    result “$file” &
  done
}

########################### MAIN ###############################
# call all functions
# user input:
initialize_globals
get_user_input "$@"
check_flags

touch results/pipeline.log
exec 3>&1 1>>"results/pipeline.log" 2>&1

echo "Starting Pipeline" 1>&3
echo "Mapping reads to reference" 1>&3
if [ "$single_end" == true ]; then
  single_end_mapping
else
  # start workflow:
  map_to_ref
  wait  
fi

echo "Filtering mapped reads" 1>&3
keep_mapped_reads
wait
sort_bam
wait 
index_bam
wait


echo "Calculating depth" 1>&3
depth
wait

echo "Determining consensus sequences" 1>&3
consensus
if [ "$single_end" == false ]; then
  change_fasta_header
fi
echo "Align sequences to the reference sequence" 1>&3
mafft_alignment

over_50
echo "Generate report" 1>&3
results_report
wait
muttable
wait

if [ "$q" == true ]; then
  echo "Quasispecies analysis" 1>&3
  mkdir pileup
  pileup
  wait
  python "$path"/quasispecies.py
fi

conda deactivate
echo "pipeline finished! (:" 1>&3
echo "Pipeline log can be found in results/pipeline.log" 1>&3