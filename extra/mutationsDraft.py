from Bio import SeqIO
import csv
import pandas as pd
from sys import argv


# TODO: change from Hard Code! maybe duble mutations rows to get uniq values of lineages
muts = {
    "UK": "B.1.1.7 - UK",
    "SA": "B.1.351 - SA",
    "Manaus": "P.1 - Manaus",
    "Rio": "P.2- Rio de jeneiro",
    "Mink": "Mink cluster 5",
    "Aus": "Australian 501 variant"
}

mutTable = pd.read_csv("novelMutTable.csv")
mutTable["REF"] = mutTable["REF"].apply(lambda x: x.upper())


def getMuts(seq, muttable = mutTable):
    mutations_list = []
    S_mutations = []
    for line in muttable.iterrows():
        pos = int(line[1][5]) - 1
        alt = seq[pos]
        ref = line[1][6]
        gene = line[1][2]
        mutname = line[1][1]
        if alt != ref and alt != 'N':
            mutations_list.append(mutname)
        elif alt == 'N' and gene == "S":
            S_mutations.append(mutname)

    return mutations_list, S_mutations


def mutations_col(seq, muttable = mutTable):
    muts_list, s_muts_list = getMuts(seq, muttable)
    mutations = ','.join([str(x) for x in muts_list])
    s_mutations = ','.join([str(x) for x in s_muts_list])
    return mutations + ' /S: ' + s_mutations


def isVar(seq, mut, muttable = mutTable, mutdict = muts):
    lineage = mutdict[mut]
    is_var = True
    for line in muttable[muttable["Lineage"].str.contains(lineage)].iterrows():
        pos = int(line[1][5])-1
        alt = seq[pos]
        ref = line[1][6]
        if alt == ref or alt == 'N':
            is_var = False
            break
    return is_var


def knownVar(seq, muttable = mutTable):
    """
    Full Variant! Has all mutations of variant.
    :param seq: sample aligned sequence
    :param muttable: mutations table
    :return: Variant name
    """
    vars = []
    if isVar(seq, "UK", muttable):
        vars.append("UK")
    if isVar(seq, "SA", muttable):
        vars.append("SA")
    if isVar(seq, "Manaus", muttable):
        vars.append("Manaus")
    if isVar(seq, "Mink", muttable):
        vars.append("Mink")
    if isVar(seq, "Aus", muttable):
        vars.append("Australian")
    return vars if vars else "-"

# argv[1] = NGS29_aligned.fasta
# argv[2] = output csv table of mutations in samples

fastadict = SeqIO.to_dict(SeqIO.parse(argv[1], 'fasta')) # fastadict -> {id: seq object}
fastadict.pop('NC_045512.2', None)  # remove refseq from dictionary (if does not exist, will do nothing - no Error)
fastadict.pop('REF_NC_045512.2', None)

samples = []
for fa in fastadict:
    samples.append(fa)

table = pd.DataFrame()
table["sample"] = fastadict.keys()
table["seq"] = table["sample"].map(fastadict)  # to delete afterwards
table["Known Variant"] = table.apply(lambda row: knownVar(row["seq"]), axis=1)
table["Muts"] = table.apply(lambda row: mutations_col(row["seq"]), axis=1)



table.to_csv("NGS29MUTS.csv")