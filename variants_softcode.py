from Bio import SeqIO
import csv
import pandas as pd
from sys import argv

# TODO: Dictionary of dictionaries. maybe create mutations database more convenient
# Mayby create mutations as json??? and list all lineages in the same row. eah mutation can be a class instant that will
# TODO: 70% of certain mutation is enough but only if did not fully match other.
# mark partial match, sort it, and if top is > 70% choose it
# no hard code!

# Dict of lineages:

alignment_file = argv[1]
mutTable = pd.read_csv("novelMutTable_2.csv")
# mutTable = pd.read_csv("/data/projects/Dana/scripts/covid19/novelMutTable.csv")

mutTable["REF"] = mutTable["REF"].apply(lambda x: x.upper())
mutTable["mut"] = mutTable["mut"].apply(lambda x: x.upper())

mutations_by_lineage = {lin: [] for lin in pd.unique(mutTable["lineage"])}

alignment = SeqIO.to_dict(SeqIO.parse(alignment_file, 'fasta'))  # fasta-dict -> {id: seq object}
alignment.pop('NC_045512.2', None)
alignment.pop('REF_NC_045512.2', None)


