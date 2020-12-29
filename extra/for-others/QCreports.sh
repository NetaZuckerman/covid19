#!/bin/bash

################# QC reports ####################
# Dana Bar-Ilan, Nov. 18th 2020
# Bioinformatics |  Central Virology Laboratory, Israel Ministry of Health

[ $# -eq 0 ] && { echo "Usage: $0 {input path} {output path}">&2; exit 1; }
# $1 location of fastq.gz files  $2 location to drop reports in
input=$1
output=$2

# fastqc
for file in "$input"*.fastq.gz; do
  fastqc "$file" --outdir="$output" -q
done

# multiqc
multiqc --interactive "$output" -o "$output"
