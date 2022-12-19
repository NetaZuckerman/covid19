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
      -r|--reference-sequence)
        shift
        refseq="$1"
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
      --invr)
        invr=true
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
[[ -z "$refseq" ]] && echo "Please provide reference sequence (-r|--reference-sequence)" && exit 0
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


echo "Run Nextclade" 1>&3
conda activate nextstrain
nextclade run --input-dataset $SCRIPT_DIR/nextclade --output-tsv results/nextclade.tsv --output-fasta alignment/all_aligned.fasta --output-insertions alignment/insersions.csv "$sequences"
#nextclade run "$sequences" --input-dataset $SCRIPT_DIR/nextclade --output-tsv results/nextclade.tsv --output-fasta alignment/all_aligned.fasta
conda deactivate

if [ "$invr" == true ]; then
  echo "INVR" 1>&3
  python "$path"/invr.py "$path"/covid19_regions.csv "$refseq" alignment/all_aligned.fasta results/nextclade.tsv results/invr.csv
fi

conda activate CoronaPipeline
echo "Variants analysis" 1>&3
#python "$path"/translated_table.py "alignment/all_aligned.fasta results/AA_muttable.xlsx "$path"/regions.csv "$path"/mutationsTable.xlsx
python "$path"/variants.py alignment/all_aligned.fasta results/variants.csv results/nextclade.tsv "$path"/mutationsTable.xlsx "$noRecombinants"

echo "Extra mutations report" 1>&3
python "$path"/translate_extras.py "$path"/covid19_regions.csv results/non_variant_mutations.csv "$refseq" alignment/all_aligned.fasta

echo "Nucleotides mutation table report" 1>&3
python "$path"/MutTable.py alignment/all_aligned.fasta results/nuc_muttable.xlsx "$path"/mutationsTable.xlsx

if [ "$q" == true ]; then
  echo "Quasispecies analysis" 1>&3
  mkdir pileup
  pileup
  wait
  python "$path"/quasispecies.py
fi
echo "Finished minipipeline (:" 1>&3
echo "Minipipline log can be found in results/minipipline.log" 1>&3
