#!/bin/bash

input_dir=$1
output_dir=$2
fastqc $input_dir/?*.fastq.gz --outdir=$output_dir # NOTICE: EXPECTING fastq.gz POSTFIX

