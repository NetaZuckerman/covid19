#!/bin/bash
# 18.05.2021 Dana Bar Ilan

# respond do ctrl+c to abort run:
trap "kill 0" EXIT
# allow conda activation from scrit
eval "$(conda shell.bash hook)"

# help function: (usage)
function usage() {
    cat <<EOF
  usage: $0 [-h] [-i | --suquences] [-r | --reference-sequence]

  --unaligned                 [unaligned MULTI FASTA]
  --aligned                   [aligned MULTI FASTA
  -r | --reference-sequence   [FASTA]
EOF
exit 0
}

function get_user_input() {
  while (( "$#" )); do
    case "$1" in
      --unaligned)
        shift
        unaligned="$1"
        shift
        ;;
      --aligned)
        shift
        aligned="$1"
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
path=$(dirname "${0}")
get_user_input "$@"
# check input:

mkdir -p results_new # -p: create only if doesnt already exist

conda activate nextstrain
# nextclade (-t: tsv output)
nextclade -i "$unaligned" -t results_new/nextclade.tsv

conda deactivate

conda activate pangolin
pangolin "$unaligned" --outfile results_new/pangolinClades.csv
conda deactivate

conda activate CoronaPipeline
python "$path"/MutTable.py "$aligned" results_new/nuc_muttable.xlsx "$path"/mutationsTable.xlsx
python "$path"/translated_table.py "$aligned" results_new/AA_muttable.xlsx "$path"/regions.csv "$path"/mutationsTable.xlsx
python "$path"/variants.py "$aligned" results_new/variants.csv results_new/pangolinClades.csv results_new/nextclade.tsv "$path"/mutationsTable.xlsx QC/report.txt

echo "finished minipipeline (:"
