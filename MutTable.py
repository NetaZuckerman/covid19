import pandas as pd
from Bio import SeqIO  # biopython
from sys import argv   # no need for special download
"""
produce report of nucleotide monitored mutations, based on alignment to reference sequence. 
does not include insertions (for now). 
"""

# input: multi-fasta, name of output file


def highlight_row(row):
    """
    highlight mutations in yellow and refseq in silver.
    :param row: row of dataframe.
    :return: parallel row of colors.
    """
    mut = row.mut
    color_list = [""] * 7 + ["background-color: silver"] * 2 + ["background-color: yellow" if i == mut and i != 'N' else
                                                                "" for i in row[9:]]
    return color_list


# df = pd.read_csv("novelMutTable.csv")
df = pd.read_csv("/data/projects/Dana/scripts/covid19/novelMutTable.csv")  # name of mutations table

# argv[1] : input multi-fasta file (aligned by augur!)
# argv[2] : output csv table of mutations in samples

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

df["REF"] = df["REF"].str.upper()  # make sure all reference nucleotides from table are upper case.
df = df[df.type != 'Insertion']
varcol = df.apply(lambda row: row[8:].unique(), axis=1)
df.insert(6, "var", varcol)
df = df.sort_values(by=["gene", "lineage"], ascending=[False, True]) # check!

df = df.rename(columns={'pos': 'nuc pos', 'nucleotide': 'nuc name', 'AA': 'name'})

# change order of columns
sorted_cols = ['nuc pos', 'nuc name', 'type', 'gene', 'var', 'name', 'lineage', 'REF', 'mut']
df = df[sorted_cols + [col for col in df.columns if col not in sorted_cols]]
# write to file
df.style.apply(highlight_row, axis=1).to_excel(argv[2], index=False)
