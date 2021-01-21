import pandas as pd
from Bio import SeqIO
from sys import argv

df = pd.read_csv("/data/projects/Dana/scripts/covid19/novelMutTable.csv")  # maybe pull from other source in the server?


# argv[1] = input multi-fasta file (aligned by augur!)
# argv[2] = output csv table of mutations in samples

fastadict = SeqIO.to_dict(SeqIO.parse(argv[1], 'fasta')) # fastadict -> {id: seq object}
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
df.insert(7, "var",varcol)
df.to_csv(argv[2])