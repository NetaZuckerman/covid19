import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
# TODO: make general to pipeline - accept user input? or approach to aligned file in server.
df = pd.read_csv("mutTable.csv")  # maybe pull from other source in the server?

# create Series, index: fasta header, value: sequence
# fastaDict = {record.id: record.seq for record in SeqIO.parse('ngs1/all_aligned.fasta', 'fasta')}

fastadict = SeqIO.to_dict(SeqIO.parse('ngs1/allAligned.fasta', 'fasta')) # fastadict -> {id: seq object}

fastadict.pop('NC_045512.2', None)  # remove refseq from dictionary (if does not exist, will do nothing - no Error)
fastadict.pop('REF_NC_045512.2', None)

for file, seqrecord in fastadict.items():
    seq = seqrecord.seq
    mutpositions = []
    for pos in df["nuc pos"]:
        x = seq[int(pos)-1]
        mutpositions.append(x)
    df[file] = mutpositions


# df['val'] = df.apply(lambda row: print(row), axis=1)
varcol = df.apply(lambda row: row[5:].unique(), axis=1)
df.insert(5, "var",varcol)
df.to_csv('updatedmuttable-allAligned.csv')