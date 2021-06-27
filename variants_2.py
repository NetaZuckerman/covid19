from Bio import SeqIO
import csv
import pandas as pd
from sys import argv
import os

"""
create variants table - for each sample in fasta multiple alignment file a covid variant is decided if found.
including pangolin and nextclades variants as well. 
"""

# get user input
alignment_file = argv[1]
output_file = argv[2]
pangolin_file = argv[3]
clades_path = argv[4]  # nextclade tsv
excel_path = argv[5]  # mutations table path
qc_report_path = argv[6]

try:
    qc = pd.read_csv(qc_report_path, sep='\t')
    qc['sample'] = qc['sample'].apply(str)
except FileNotFoundError:
    print("QC File does not exist")
    qc = pd.DataFrame()
# load pangolin + nextclade outputs, mutations table.
try:
    pangolinTable = pd.read_csv(pangolin_file)
    pangolinTable['taxon'] = pangolinTable['taxon'].apply(str)
except FileNotFoundError:
    print("Pangolin File Not Found")
    pangolinTable = pd.DataFrame()

clades_df = pd.read_csv(clades_path, sep='\t')

# excel muttable:
excel_mutTable = pd.read_excel(excel_path, sheet_name=None, engine='openpyxl')
for name in excel_mutTable:
    excel_mutTable[name]['lineage'] = name  # add a lineage column to all variant's tables

mutTable = pd.concat(excel_mutTable.values(), ignore_index=True)
mutTable['Mutation type'] = mutTable['Mutation type'].str.lower()
mutTable = mutTable[mutTable['Mutation type'] != 'insertion']
# prepare mutations table dataframe
mutTable["Reference"] = mutTable["Reference"].str.upper()  # make sure all upper case for safe comparisons
mutTable["Mutation"] = mutTable["Mutation"].str.upper()

# prepare multiple alignment dictionary (key: sample name, val: SeqIO record)
alignment = SeqIO.to_dict(SeqIO.parse(alignment_file, 'fasta'))
alignment.pop('NC_045512.2', None)
alignment.pop('REF_NC_045512.2', None)

# prepare nextclade dataframe
clades_df = clades_df[['seqName', 'aaSubstitutions', 'aaDeletions', 'clade', 'insertions']]
clades_df = clades_df.rename(columns={'seqName': 'sample'})
clades_df['sample'] = clades_df['sample'].apply(str)
clades_df = clades_df.fillna('')
# create dict of aasubstitutions and aadeletions. key=sample id, value=list of substitutions. to keep for final table.
aa_substitution_dict = {}
aa_deletions_dict = {}
insertions_dict = {}
for sample in clades_df['sample']:  # change appearance from nextclade format to: Mutation(Gene)
    aasubs = clades_df[clades_df['sample'] == sample].aaSubstitutions.values.tolist()
    aadels = clades_df[clades_df['sample'] == sample].aaDeletions.values.tolist()
    insertions = clades_df[clades_df['sample'] == sample].insertions.values.tolist()
    aa_substitution_dict[sample] = [f"{x.split(':')[1]}({x.split(':')[0]})" for x in aasubs[0].split(',')] \
        if (aasubs and aasubs != ['']) else ''
    aa_deletions_dict[sample] = [f"{x.split(':')[1]}({x.split(':')[0]})" for x in aadels[0].split(',')] \
        if (aadels and aadels != ['']) else ''
    insertions_dict[sample] = insertions[0].split(',')

samples_mutations = {id: [] for id in alignment}
# samples_s_not_covered = {id: [] for id in alignment}
samples_not_covered = {id: [] for id in alignment}
unexpected_mutations = {id: [] for id in alignment}

# iterate over all samples in multi-fasta and over all mutations in table, and check value of each mutation
mutTable = mutTable.dropna(thresh=10)  # lose empty rows that might get there by mistake in concat
for sample, record in alignment.items():
    for (idx, row) in mutTable.iterrows():
        pos = int(row.loc['Position']) - 1  # mutation position
        alt = record.seq[pos]  # fasta value in position
        ref = row.loc['Reference']  # reference in position
        table_mut = row.loc['Mutation']  # mutation  according to table
        gene = row.loc['protein']
        mutation_name = str(row.loc['variant'])
        if alt == table_mut:  # mutation exists in sequence
            samples_mutations[sample].append(mutation_name)
        elif alt == 'N':  # if position not covered in sequence
            samples_not_covered[sample].append(mutation_name)
        else:  # alt is not the expected mut. and is covered in sequencing (not N)
            unexpected_mutations[sample].append(mutation_name + "(alt:" + alt + ")")

mutations_by_lineage = {key: excel_mutTable[key].variant.tolist() for key in
                        excel_mutTable}  # mutations table with only lineage: mutation name

