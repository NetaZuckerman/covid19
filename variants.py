"""
variants script similare to variants_softcode.py but with optional arguments
for example: to not include pangolin table
etc.
"""
import csv
from sys import argv
import pandas as pd
from Bio import SeqIO
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Create Variants Table")
    parser.add_argument('aligned', metavar='<aligned fasta>', type=str, nargs=1,
                        help='aligned multi-fasta file')
    parser.add_argument('-o', metavar='path', type=str, nargs=1,
                        help='output path file')


    args = parser.parse_args()