from Bio import SeqIO
import pandas as pd
from sys import argv
# get uk fasta
#
# uk_headers = open("../uk_names.txt")
# uk_names = [line.strip() for line in uk_headers]
# uk_names_ngs = [x+"ngs" for x in uk_names]
# uk_names_ngs2 = [x + "-ngs" for x in uk_names]
#
# all_fasta = SeqIO.to_dict(SeqIO.parse("../NGS29_aligned.fasta", 'fasta'))
#
# found = []cd
# with open("uk_all_aligned.fasta", 'w') as outfile:
#     for name, seqrecord in all_fasta.items():
#         if name in uk_names:
#             SeqIO.write(seqrecord, outfile, "fasta")
#             found.append(name)
#         elif name in uk_names_ngs:
#             SeqIO.write(seqrecord, outfile, "fasta")
#             found.append(name)
#         elif name in uk_names_ngs2:
#             SeqIO.write(seqrecord, outfile, "fasta")
#             found.append(name)
#
# with open("brits_headers.txt", 'w') as headers:
#     for f in found:
#         headers.write(f'{f}\n')
#
# uk_headers.close()

############################################
# muts - adding another columns to table in order to classify mutations
# df = pd.read_csv("../novelMutTable.csv")
#
# # argv[1] = input multi-fasta file (aligned by augur!)
# # argv[2] = output csv table of mutations in samples
#
# fastadict = SeqIO.to_dict(SeqIO.parse(argv[1], 'fasta')) # fastadict -> {id: seq object}
# fastadict.pop('NC_045512.2', None)  # remove refseq from dictionary (if does not exist, will do nothing - no Error)
# fastadict.pop('REF_NC_045512.2', None)
#
# for file, seqrecord in fastadict.items():
#     seq = seqrecord.seq
#     mutpositions = []
#     for pos in df["nuc pos"]:
#         x = seq[int(pos)-1]
#         mutpositions.append(x)
#     df[file] = mutpositions
#
# df["REF"] = df["REF"].str.upper()
#
# # df['val'] = df.apply(lambda row: print(row), axis=1)
# varcol = df.apply(lambda row: row[5:].unique(), axis=1)
# df.insert(7, "alts", varcol)
# df = df.sort_values("Lineage")
#
# # column -> Var! to determine variant
# # gather all muts of lineage
# # problem: lineage need to be parsed!
# british_var = pd.unique([x for x in df["Lineage"] if "B.1.1.7 - UK" in x])
# sa_var = pd.unique([x for x in df["Lineage"] if "B.1.351 - SA" in x])
# barzil_rio_var = pd.unique([x for x in df["Lineage"] if "P.2- Rio de jeneiro" in x])
# brazil_manus_var = pd.unique([x for x in df["Lineage"] if "P.1 - Manaus" in x])
# mink_var = pd.unique([x for x in df["Lineage"] if "Mink cluster 5" in x])
# australian_var = pd.unique([x for x in df["Lineage"] if "Australian 501 variant" in x])

# df.to_csv(argv[2])
#
# with open("../gisaid_Adi_May25.fasta") as adistern, open("adinames.txt",'w') as outfile:
#     records = SeqIO.parse(adistern, "fasta")
#     for recorn in records:
#         outfile.write(f'{recorn.id}\n')
muttable = pd.read_csv("../novelMutTable.csv")
print(pd.unique(muttable.lineage))


