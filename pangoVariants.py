from Bio import SeqIO
import csv
import pandas as pd
from sys import argv

alignment_file = argv[1]  # ?
output_file = argv[2]
pangolin_file = argv[3]
pangolinTable = pd.read_csv(pangolin_file)
# option for other mutations table in the future (like 'S' specific, etc..)
if len(argv) > 4:
    muttable_path = argv[4]
else:
    muttable_path = "/data/projects/Dana/scripts/covid19/novelMutTable.csv"
mutTable = pd.read_csv(muttable_path)
mutTable["REF"] = mutTable["REF"].apply(lambda x: x.upper())
mutTable["mut"] = mutTable["mut"].apply(lambda x: x.upper())
mutTable = mutTable[mutTable.type != 'Insertion']  # Ignore insertions for now

alignment = SeqIO.to_dict(SeqIO.parse(alignment_file, 'fasta'))  # fasta-dict -> {id: seq object}
alignment.pop('NC_045512.2', None)
alignment.pop('REF_NC_045512.2', None)

pangolinTable = pangolinTable.drop('REF_NC_045512.2')
