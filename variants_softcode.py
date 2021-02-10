from Bio import SeqIO
import csv
import pandas as pd
from sys import argv

# TODO: Dictionary of dictionaries. maybe create mutations database more convenient

# Dict of lineages:

mutTable = pd.read_csv("novelMutTable.csv")
# mutTable = pd.read_csv("/data/projects/Dana/scripts/covid19/novelMutTable.csv")
mutTable["REF"] = mutTable["REF"].apply(lambda x: x.upper())
mutTable["Mut"] = mutTable["REF"].apply(lambda x: x.upper())

lineages = pd.unique(mutTable["Lineage"])  # np. array

