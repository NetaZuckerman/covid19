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

  -i | --sequences            [unaligned MULTI FASTA]
  -r | --reference-sequence   [FASTA]
EOF
exit 0
}

function get_user_input() {
  while (( "$#" )); do
    case "$1" in
      -i|--sequences)
        shift
        unaligned="$1"
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
echo "checking input.."
[[ -z "$unaligned" ]] && echo "Please provide multi-fasta sequence (-i|--sequences)" && exit 0
[[ -z "$refseq" ]] && echo "Please provide reference sequence (-r|--reference-sequence)" && exit 0
echo "all input provided. continuing."
mkdir -p {alignment,results} # -p: create only if doesnt already exist

conda activate nextstrain
# nextclade (-t: tsv output)
nextclade -i "$unaligned" -t results/nextclade.tsv > /dev/null 2>&1
# align multifasta to reference sequence using augur align:
augur align \
--sequences "$unaligned" \
--reference-sequence "$refseq" \
--output alignment/all_aligned.fasta

conda deactivate

conda activate pangolin
pangolin "$unaligned" --outfile results/pangolinClades.csv
conda deactivate

conda activate CoronaPipeline
python "$path"/MutTable.py alignment/all_aligned.fasta results/nuc_muttable.xlsx "$path"/mutationsTable.xlsx
python "$path"/translated_table.py alignment/all_aligned.fasta results/AA_muttable.xlsx "$path"/regions.csv "$path"/mutationsTable.xlsx
python "$path"/variants.py alignment/all_aligned.fasta results/variants.csv results/pangolinClades.csv results/nextclade.tsv "$path"/mutationsTable.xlsx QC/report.txt

echo "finished minipipeline (:"
