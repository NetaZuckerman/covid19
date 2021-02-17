import csv
from sys import argv
import pandas as pd
# craete s_mutations table from full table

full_table_path = argv[1]
s_table_path = argv[2]

full_table = pd.read_csv(full_table_path)

s_table = full_table[full_table.gene == 'S']

s_table.to_csv(s_table_path)
