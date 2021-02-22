from sys import argv
from Bio import SeqIO
import pandas as pd
from Bio.Seq import translate

# user input:
aligned_fasta_path = argv[1]
outfile_path = argv[2]
regions_table_path = argv[3]

# get start and end of region
# don't forget -1 in position

# 1. load sequences and tables
regionsTable = pd.read_csv(regions_table_path)
multifasta = SeqIO.to_dict(SeqIO.parse(aligned_fasta_path, 'fasta'))
mutTable = pd.read_csv("novelMutTable.csv")  # TODO: change to server's path
# 2. keep only non-synonymous mutations
mutTable = mutTable[mutTable.type == 'non -Synonymous ']
finalTable = mutTable
# 3. iterate over mutations and create
for sample, record in multifasta.items():
    seq = record.seq
    sample_muts = []
    for mut in mutTable.iterrows():
        pos = mut[1][5]
        gene = mut[1][2]
        aa = mut[1][1]
        aa_number = int(aa[1:-1])  # strip letters of both sides
        aa_from = aa[0]
        aa_to = aa[-1]
        region_start = int(regionsTable[regionsTable.id == gene].start.values[0])
        region_end = int(regionsTable[regionsTable.id == gene].end.values[0])
        region_aa = seq[region_start-1:region_end].translate()
        if len(region_aa) < aa_number:
            print(aa)
        var = region_aa[aa_number-1]
        if var == aa_to:
            # mutation
            print('{}: {}'.format(sample, aa))
# TODO: TO FIX -> --- not translated. maybe translate with my own function??
# finalTable.REF = finalTable.REF.str.upper()


