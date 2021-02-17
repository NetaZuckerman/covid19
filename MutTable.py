import pandas as pd
from Bio import SeqIO
from sys import argv
# import xlsxwriter
# import numpy as np

# df = pd.read_csv("novelMutTable.csv")
df = pd.read_csv("/data/projects/Dana/scripts/covid19/novelMutTable.csv"


# argv[1] = input multi-fasta file (aligned by augur!)
# argv[2] = output csv table of mutations in samples

fastadict = SeqIO.to_dict(SeqIO.parse(argv[1], 'fasta')) # fastadict -> {id: seq object}
fastadict.pop('NC_045512.2', None)  # remove refseq from dictionary (if does not exist, will do nothing - no Error)
fastadict.pop('REF_NC_045512.2', None)

samples = []
for file, seqrecord in fastadict.items():
    seq = seqrecord.seq
    samples.append(file)
    mutpositions = []
    for pos in df["pos"]:
        x = seq[int(pos)-1]
        mutpositions.append(x)
    df[file] = mutpositions
    # df[file] = [seq[int(pos)-1] for pos in df["nuc pos"]]

df["REF"] = df["REF"].str.upper()

# df['val'] = df.apply(lambda row: print(row), axis=1)
varcol = df.apply(lambda row: row[5:].unique(), axis=1)
df.insert(7, "var",varcol)
df = df.sort_values(by=["gene", "Lineage"], ascending=[False, True]) # check!


df.to_csv(argv[2])
#
# def function(x):
#     df1 = pd.DataFrame('background-color: ', index=x.index, columns=x.columns)
#     for s in samples:
#         m1 = x[s] != x["REF"]
#         df1[s] = np.where(m1,'background-color: {}'.format('yellow'), df1[s])
#
# df.style.apply(function,)