#!/bin/bash

# for all mapped.sorted.bam files produce consensus report using ivar
# example command:  (source: https://andersen-lab.github.io/ivar/html/manualpage.html)
#  samtools mpileup -d 1000 -A -Q 0 test.bam | ivar consensus -p test -q 20 -t 0
# samtools mpileup
# -d max depth  -A count orphands (Do not skip anomalous read pairs in variant calling.
#                  Anomolous read pairs are those marked in the FLAG field as paired in
#                  sequencing but without the properly-paired flag set.)
# ivar consesus
# -p prefix output files (required)
# -q min quality score to count
# -t min frequency
# -m min depth to call consensus

mkdir -p ivarCNS5
for file in BAM/*.mapped.sorted.bam; do
  samtools mpileup -A "$file" | ivar consensus -m 5 -p ivarCNS5/`basename "$file" .mapped.sorted.bam`
done
