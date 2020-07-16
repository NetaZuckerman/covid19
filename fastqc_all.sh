#!/bin/bash
out=$1
fastqc fastq/raw/*.fastq.gz --outdir="$out"

