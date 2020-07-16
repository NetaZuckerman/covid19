#!/bin/bash
in=$1
out=$2
fastqc "$in" --outdir="$out"
