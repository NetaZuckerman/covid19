from collections import namedtuple
from sys import argv
import pandas as pd
from Bio import SeqIO
import os
import subprocess


class NucleotidesDistribution:
    def __init__(self, a_nuc, c_nuc, g_nuc, t_nuc, n_nuc):
        self.a = a_nuc
        self.c = c_nuc
        self.g = g_nuc
        self.t = t_nuc
        self.n = n_nuc

    def __str__(self):
        return 'a:{a},c:{c},g:{g},t:{t},n:{n}'.format(a=self.a, c=self.c, g=self.g, t=self.t,
                                                      n=self.n)

# now: bam read-counter output file
# alternatively:
# bcftools mpileup --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
#     BAM/602_S16_L001_001.mapped.sorted.bam --no-reference | bcftools query --format '%CHROM\t%POS\t%ALT[\t%AD]\n' \
#     > somefile.vcf
# then parse vcf somehow

# each sample has a whole file
# iterate over all files and return a report
# consider writing a python script for whole pipeline

# df = pd.read_csv("/data/projects/Dana/scripts/covid19/novelMutTable.csv")  # maybe pull from other source in the server?
df = pd.read_csv("novelMutTable.csv")

# ReadCount = namedtuple('ReadCount', ['depth', 'ref', 'A', 'C', 'G', 'T', 'N'])  # needed?
# for now: one file
# later: iterate over all files

final_table = pd.DataFrame({
    'AA': df['AA'],
    'gene': df['gene'],
    'Lineage': df['Lineage'],
    'Type': df['Type'],
    'REF': df['REF'],
    'pos': df['nuc pos']
})

# iterate over all readcounts, file names are sample names
path = argv[1] if argv[1].endswith('/') else argv[1]+'/'
# positions = pd.DataFrame({
#     'pos': pd.unique(final_table.pos)
# })

for filename in os.listdir(path):
    with open(path+filename) as file:
        positions = {p: "" for p in pd.unique(final_table.pos)}
        for line in file:
            line = line.split('\t')
            pos = int(line[1])
            if pos in positions.keys():
                ref = line[2]
                total_depth = int(line[3])
                a = round((int(line[5].split(':')[1]) / total_depth) * 100, 2)
                c = round((int(line[6].split(':')[1]) / total_depth) * 100, 2)
                g = round((int(line[7].split(':')[1]) / total_depth) * 100, 2)
                t = round((int(line[8].split(':')[1]) / total_depth) * 100, 2)
                n = round((int(line[9].split(':')[1]) / total_depth) * 100, 2)
                positions[pos] = NucleotidesDistribution(a, c, g, t, n)

        print(filename.split('_')[0])
        final_table[filename.split('_')[0]] = final_table.apply(lambda row: positions[row.pos] if positions[row.pos] else 'None', axis=1)
print(final_table)

final_table.to_csv("readcountTable.csv")
# for now: iteration for one sample
# argv[1] = readcont file

# for p in df["nuc pos"]:
#     df["sample"] = positions[p] ######################## when multiple files - add sample name! TEMP!!!!
#

# output file: each sample is a column

