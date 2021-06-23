import pandas as pd
from Bio import SeqIO  # biopython
from sys import argv   # no need for special download
"""
produce report of nucleotide monitored mutations, based on alignment to reference sequence. 
does not include insertions (for now). 
"""

# input:
# argv[1] : input multi-fasta file (aligned by augur!)
# argv[2] : output xlsx table of mutations in samples. must end with xlsx.
# argv[3] : mutations table path



def highlight_row(row):
    """
    highlight mutations in yellow and refseq in silver.
    :param row: row of dataframe.
    :return: row of colors.
    """
    mut = row.mut
    color_list = [""] * 7 + ["background-color: silver"] * 2 + ["background-color: yellow" if i == mut and i != 'N' else
                                                                "" for i in row[9:]]
    return color_list


excel_table = pd.read_excel(argv[3], sheet_name=None, engine='openpyxl')  # load excel mutations table
for name in excel_table:
    excel_table[name]['lineage'] = name  # add lineage column to all tables before joining them

df = pd.concat(excel_table.values(), ignore_index=True)
df = df[['Position', 'Reference', 'Mutation', 'protein',
                     'variant', 'Mutation type', 'lineage', 'annotation']]  # select subset of columns
# compress identical mutations into one line and concat lineage names in the lineage column:
# df = df.groupby(  # to create compressed table:
#     ['Position', 'Reference', 'Mutation', 'protein', 'variant', 'Mutation type', 'annotation'], as_index=False).agg(
#     {'lineage': ';'.join}
# )

fastadict = SeqIO.to_dict(SeqIO.parse(argv[1], 'fasta'))  # fastadict -> {id: seq object}
fastadict.pop('NC_045512.2', None)  # remove refseq from dictionary (if does not exist, will do nothing - no Error)
fastadict.pop('REF_NC_045512.2', None)

samples = []  # list to aggregate a samples column to.
# iterate over fasta records of aligned fasta file and get values in mutations positions
df = df.dropna(how=6)
df.to_csv('concat_table.csv')
for file, seqrecord in fastadict.items():
    seq = seqrecord.seq
    samples.append(file)  # keep sample names in list
    mutpositions = []
    for pos in df["Position"]:  # for each mutations position get the value from fasta record (-1 bcs index starts from 0)
        x = seq[int(pos)-1]
        mutpositions.append(x)
    df[file] = mutpositions  # add sample as column

df = df[df['Mutation type'].str.lower() != 'insertion']  # ignore insertions # ignore case to avoid mistakes
varcol = df.apply(lambda row: row[9:].unique(), axis=1)  # add var columns- show which values in row (from col9 forward
# because 9 first columns are the muations table)
df.insert(6, "var", varcol)
df = df.sort_values(by=["lineage", "protein"], ascending=[True, False])  # sort by lineage and gene
df = df.rename(columns={'Position': 'nuc pos', 'Mutation type': 'type', 'protein': 'gene',
                        'variant': 'name', 'Reference': 'REF', 'Mutation': 'mut'})

# change order of columns
sorted_cols = ['nuc pos', 'type', 'gene', 'var', 'name', 'lineage', 'annotation', 'REF', 'mut'] # re-order columns
df = df[sorted_cols + [col for col in df.columns if col not in sorted_cols]]  # re-order columns

# write to file
# add highlights with designated function in the top of the page
try:
    df.style.apply(highlight_row, axis=1).to_excel(argv[2], index=False)
except ValueError:
    df.style.apply(highlight_row, axis=1).to_excel('MutTable.xlsx', index=False)
    print("output file is not valid - must end with .xlsx")
    print("output was written to MutTable.xlsx instead of %s" % (argv[2]))
