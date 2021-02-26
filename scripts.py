import csv
from sys import argv
import pandas as pd
from Bio import SeqIO
# craete s_mutations table from full table

ref = open("REF_NC_045512.2.fasta")
getRef = SeqIO.to_dict(SeqIO.parse(ref, "fasta"))

for id, record in getRef.items():
    print("28248-28254 in REF: ")
    print(record.seq[28247:28254])