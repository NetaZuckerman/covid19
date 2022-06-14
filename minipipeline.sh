#!/bin/bash
# 18.05.2021 Dana Bar Ilan

# respond do ctrl+c to abort run:
trap "kill 0" EXIT
# allow conda activation from scrpit
eval "$(conda shell.bash hook)"
#source ~/miniconda3/etc/profile.d/conda.sh

# help function: (usage)
function usage() {
    cat <<EOF
  usage: $0 [-h] [-i | --suquences] [-r | --reference-sequence]

  -i | --sequences            [MULTI FASTA]
  -r | --reference-sequence   [FASTA]

  optional:
  -n| --newNextclade              use nextclade version 1.3.0
  --dontAlign                     input aligned fasta and skip alignment stage
EOF
exit 0
}

function get_user_input() {
  while (( "$#" )); do
    case "$1" in
      -i|--sequences)
        shift
        sequences="$1"
        shift
        ;;
      --dontAlign)
        dontAlign=true
        shift
        ;;
      --noRecombinants)
        noRecombinants=true
        shift
        ;;
      -q| --quasispecies)
        q=true
        shift
        ;;     
      -p|--processes)
        shift
        num_processes="$1"
        shift
        ;;      
      -n| --newNextclade)
        newNextclade=true
        shift
        ;;
      -r|--reference-sequence)
        shift
        refseq="$1"
        shift
        ;;
      -h|--help)
        usage
        shift
        ;;
      -*|--*)  # decline any other flag
        echo "Error: Unsupported flag $1" >&2
        exit 1
        ;;
    esac
  done
}

### MAIN ###
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
path=$(dirname "${0}")
get_user_input "$@"
# check input:
echo "checking input.."
[[ -z "$sequences" ]] && echo "Please provide multi-fasta sequence (-i|--sequences)" && exit 0
[[ -z "$refseq" ]] && [ -z "$dontAlign" ] && echo "Please provide reference sequence (-r|--reference-sequence)" && exit 0
echo "all input provided. continuing."
mkdir -p {alignment,results} # -p: create only if doesnt already exist

touch results/minipipeline.log
exec 3>&1 1>>"results/minipipeline.log" 2>&1
echo "Starting Minipipeline" 1>&3

function pileup(){
  for file in BAM/*.mapped.sorted.bam; do
    ((i=i%num_processes)); ((i++==0)) && wait
    samtools mpileup "$file" -f "$refseq" > pileup/`basename "$file" .mapped.sorted.bam`.pileup &
  done
}

conda activate nextstrain
if  [ -z "$dontAlign" ] ; then
echo "Align multifasta to reference sequence" 1>&3
# align multifasta to reference sequence using augur align:
#  augur align \
#  --sequences "$sequences" \
#  --reference-sequence "$refseq" \
#  --output alignment/all_aligned.fasta
  nextalign -i "$sequences" -r "$refseq" --output-fasta alignment/all_aligned.fasta --output-insertions alignment/insertions.csv
  aligned="alignment/all_aligned.fasta"
else 
  aligned=$sequences
fi

echo "Run Nextclade" 1>&3
conda activate nextstrain
if [ "$newNextclade" == true ]; then
  nextclade --input-fasta "$aligned"  --input-dataset $SCRIPT_DIR/nextclade --output-dir nextclade/ --output-tsv results/nextclade.tsv
else 
  nextclade -i "$aligned" -t results/nextclade.tsv
fi        
conda deactivate

conda activate pangolin
echo "Run Pangolin" 1>&3
pangolin "$sequences" --outfile results/pangolinClades.csv
conda deactivate


conda activate CoronaPipeline
echo "Variants analysis" 1>&3
#python "$path"/translated_table.py "$aligned" results/AA_muttable.xlsx "$path"/regions.csv "$path"/mutationsTable.xlsx
python "$path"/variants.py "$aligned" results/variants.csv results/pangolinClades.csv results/nextclade.tsv "$path"/mutationsTable.xlsx "$noRecombinants"

echo "Extra mutations report" 1>&3
python "$path"/translate_extras.py "$path"/covid19_regions.csv results/non_variant_mutations.csv "$refseq" "$aligned"

echo "Nucleotides mutation table report" 1>&3
python "$path"/MutTable.py "$aligned" results/nuc_muttable.xlsx "$path"/mutationsTable.xlsx

if [ "$q" == true ]; then
  echo "Quasispecies analysis" 1>&3
  mkdir pileup
  pileup
  wait
  python "$path"/quasispecies.py
fi
echo "Finished minipipeline (:" 1>&3
echo "Minipipline log can be found in results/minipipline.log" 1>&3
