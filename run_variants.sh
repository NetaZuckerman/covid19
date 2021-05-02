eval "$(conda shell.bash hook)"

conda activate pangolin
pangolin alignment/all_aligned.fasta --outfile results/pangolinClades_02052021.csv
conda deactivate

conda activate nextstrain
nextclade -i alignment/all_aligned.fasta --output-json results/nextclade_02052021.json
conda deactivate

conda activate CoronaPipeline
python /data/projects/Dana/scripts/covid19/variants.py alignment/all_aligned.fasta results/variants_02052021.csv results/pangolinClades.csv results/nextclade_02052021.json
