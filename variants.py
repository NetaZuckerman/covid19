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
clades_path = argv[4]
excel_path = argv[5]  # mutations table path

# load pangolin + nextclade outputs, mutations table.
try:
    pangolinTable = pd.read_csv(pangolin_file)
except FileNotFoundError:
    print("Pangolin File Not Found")
    pangolinTable = pd.DataFrame()

clades_df = pd.read_csv(clades_path, sep='\t')

# excel muttable:
excel_mutTable = pd.read_excel(excel_path, sheet_name=None, engine='openpyxl')
mutTable_copy = excel_mutTable.copy()
for name, frame in mutTable_copy.items():
    frame['Mutation type'] = frame['Mutation type'].str.lower()  # make sure lower case to prevent mistakes
    excel_mutTable[name] = frame[frame['Mutation type'] != 'insertion']  # ignore insertions for now
    excel_mutTable[name]['lineage'] = name  # add a lineage column to all variant's tables

mutTable = pd.concat(excel_mutTable.values(), ignore_index=True)

# prepare mutations table dataframe
mutTable["Reference"] = mutTable["Reference"].str.upper()  # make sure all upper case for safe comparisons
mutTable["Mutation"] = mutTable["Mutation"].str.upper()

# prepare multiple alignment dictionary (key: sample name, val: SeqIO record)
alignment = SeqIO.to_dict(SeqIO.parse(alignment_file, 'fasta'))
alignment.pop('NC_045512.2', None)
alignment.pop('REF_NC_045512.2', None)

# prepare nextclade dataframe
clades_df = clades_df[['seqName', 'aaSubstitutions', 'aaDeletions', 'clade']]
clades_df = clades_df.rename(columns={'seqName': 'sample'})
clades_df['sample'] = clades_df['sample'].apply(str)
clades_df = clades_df.fillna('')
# create dict of aasubstitutions and aadeletions. key=sample id, value=list of substitutions. to keep for final table.
aa_substitution_dict = {}
aa_deletions_dict = {}
for sample in clades_df['sample']:  # change appearance from nextclade format to: Mutation(Gene)
    aasubs = clades_df[clades_df['sample'] == sample].aaSubstitutions.values.tolist()
    aadels = clades_df[clades_df['sample'] == sample].aaDeletions.values.tolist()
    aa_substitution_dict[sample] = [f"{x.split(':')[1]}({x.split(':')[0]})" for x in aasubs[0].split(',')] \
        if (aasubs and aasubs != ['']) else ''
    aa_deletions_dict[sample] = [f"{x.split(':')[1]}({x.split(':')[0]})" for x in aadels[0].split(',')] \
        if (aadels and aadels != ['']) else ''

samples_mutations = {id: [] for id in alignment}
# samples_s_not_covered = {id: [] for id in alignment}
samples_not_covered = {id: [] for id in alignment}
unexpected_mutations = {id: [] for id in alignment}
lineages_list = []

# iterate over all samples in multi-fasta and over all mutations in table, and check value of each mutation
mutTable = mutTable.dropna(thresh=10)  # lose empty rows that might get there by mistake in concat
for sample, record in alignment.items():
    for row in mutTable.iterrows():
        pos = int(row[1][0]) - 1  # mutation position
        alt = record.seq[pos]  # fasta value in position
        ref = row[1][1]  # reference in position
        table_mut = row[1][2]  # mutation  according to table
        gene = row[1][3] if row[1][3] else ''
        mutation_name = str(row[1][4])
        if alt != table_mut or alt == 'N':  # if position not covered in sequence OR the mutation is not the predictable
            if alt == 'N':  # mutation not covered in sequence - add to not covered dict under sample key
                samples_not_covered[sample].append(mutation_name)
            elif alt != ref and alt != 'N':  # alt is not the expected mut. and is covered in sequencing (not N)
                unexpected_mutations[sample].append(mutation_name + "(alt:" + alt + ")")  # add to unexpected mutations
        else:  # mutation exists in sequence
            samples_mutations[sample].append(mutation_name)  # accumulate all samples+mutations in dict:
            # {sampleName: mutationName}

unique_lineages = set(excel_mutTable.keys())  # get list of all unique lineages (= names of sheets of excel table)
mutations_by_lineage = {key: excel_mutTable[key].variant.tolist() for key in
                        excel_mutTable}  # mutations table with only lineage: mutation name

final_table = []
# iterate over all samples mutations and determine variants
for sample, sample_mutlist in samples_mutations.items():
    known_variant = ""
    more_muts = []
    notes = ''
    suspect = ''
    lin_percentages = {}
    lin_number = {}
    # iterate over mutations of each lineage, check how much of the lineage's mutations are covered by sample -
    # to decide on variant
    for lin, linmuts in mutations_by_lineage.items():
        temp = [m for m in linmuts]  # copy of mutations names into temp to use for calculations,
        # without harming original list
        temp_mutes = []
        for mut in sample_mutlist:  # iterate over  table mutations of sample
            if mut in temp:  # remove mutations of sample from temp list and keep in temp_mutes list.
                # eventually will check if all mutations are covered (if temp is empty)
                temp.remove(mut)
                temp_mutes.append(mut)
        if not temp:  # all mutations of that lineage exists in sample -> known variant=lineage name
            known_variant = lin
        elif len(linmuts) != len(temp):  # some of lineage's mutations exist but not all
            more_muts += temp_mutes  # keep found mutation in more_mutes
            # calculate percentage of lineage mutations found:
            lin_percentages[lin] = round(len(set(temp_mutes)) / len(set(linmuts)) * 100, 2)
            lin_number[lin] = (len(set(temp_mutes)), len(set(linmuts)))  # tuple:
            # (#lin_mutation_sample, #tot_lin_mutations)

    if known_variant:  # remove from more_mutes the mutations of the known variant to show only additional mutations
        more_muts = [x for x in more_muts if x not in mutations_by_lineage[known_variant]]
    else:  # did not find variant that has 100% mutations in sample
        max = 0
        var = ''
        flag = False  # True means known variant (at least 60% mutations)
        val = 0
        fraction = 0
        for key, tup in lin_number.items():  # pick variant that has highest muations % but at least 60 %
            # val = round((tup[0] / tup[1])*100, 2)
            val = lin_percentages[key]
            frac = f'({int(tup[0])}/{int(tup[1])})'  # keep as fraction to print to table
            if val >= 60:  # 60% -> known variant
                flag = True
            if val > max:  # make sure to choose best match and not the first over 60%
                max = val
                var = key
                fraction = frac
        if flag:  # found >60%
            known_variant = var
            # keep only additional mutations and not lineage mutations of known var.
            more_muts = [x for x in more_muts if x not in mutations_by_lineage[known_variant]]
        elif val >= 50 and var:  # > 50% mutations and < 60% because no known variant
            more_muts = [x for x in more_muts if x not in mutations_by_lineage[var]]

        if var and lin_number[var][0] >= 2:  # At least 2 mutations of lineage --> suspect variant
            suspect = f'suspect {var}: {str(lin_percentages[var])}% {fraction}'  # suspect_<lineage-name>): (%)(x/y)

    more_muts = set(more_muts)  # lose redundancies of mutations

    if not suspect and (more_muts or samples_not_covered[sample] or unexpected_mutations[sample]):
        # not specific suspect variant but some mutations exist \ not covered in sequencing - write as suspect
        suspect = 'suspect'

    # get pangolin info from table

    try:
        pangolin_clade = pangolinTable[pangolinTable.taxon == sample].lineage.values[0]
        pangolin_status = pangolinTable[pangolinTable.taxon == sample].status.values[0]
        pangolin_note = pangolinTable[pangolinTable.taxon == sample].note.values[0]
    except:
        pangolin_clade = '-'
        pangolin_status = ''
        pangolin_note = ''
    QCfail = True if pangolin_status == 'fail' else False

    # specific cases of mutations in the same location:
    temp = unexpected_mutations[sample].copy()
    for x in temp:
        short_name = x.split('(')[0]
        if 'P681' in short_name:
            if short_name == 'P681R' and 'P681H' in more_muts:
                unexpected_mutations[sample].remove(x)
            elif short_name == 'P681H' and 'P681R' in more_muts:
                unexpected_mutations[sample].remove(x)
        elif short_name == 'M234I':
            if 'M234I' in more_muts:
                unexpected_mutations[sample].remove(x)
        elif short_name == 'Q677H':
            if 'Q677H' in more_muts:
                unexpected_mutations[sample].remove(x)

    # get nextclade info from table
    nextclade = clades_df[clades_df['sample'] == sample].clade
    # create final table's line (sample line)

    more_muts_final = []
    for mut in more_muts:  # moremuts: keep only mutations in S and add (S) to mut. name
        if mutTable[mutTable.variant == mut].protein.values[0] == 'S':
            m_name = mut + "(S)"
            if m_name not in more_muts_final:  # make sure no duplicates
                more_muts_final.append(m_name)
        # else do not keep in list

    # if found known variant, write it down only if sample passed QC.
    if known_variant:
        known_variant = known_variant if not QCfail else "QC Fail"
    else:
        known_variant = "no monitored variant"


    not_covered_list = []
    for mut in samples_not_covered[sample]:
        protein_values = mutTable[mutTable.variant == mut].protein.values.tolist()
        if QCfail:
            break
        if protein_values:
            not_covered_list.append(mut + "(" + protein_values[0] + ")")
        else:
            not_covered_list.append(mut + "()")

    not_covered_list = ";".join(set(not_covered_list)) if not_covered_list else ''

    line = {
        "Sample": sample,
        "Known Variant": known_variant,
        "Suspect": suspect,
        "More Mutations": ';'.join(more_muts_final),
        "Not Covered": not_covered_list,
        "all mutations": ';'.join(aa_substitution_dict[sample]) if aa_substitution_dict and sample in
                                                                   aa_substitution_dict else 'NA',
        "aaDeletions": ';'.join(aa_deletions_dict[sample] if aa_deletions_dict and sample
                                                             in aa_deletions_dict else 'NA'),
        "nextclade": nextclade.values[0] if not nextclade.empty else '',
        "pangolin_clade": pangolin_clade,
        "status": pangolin_status,
        "pangolin-note": pangolin_note
    }
    final_table.append(line)

with open(output_file, 'w') as outfile:
    filednames = ["Sample", "Known Variant", "Suspect", "More Mutations", "Not Covered",  # 'S Not Covered'
                  "all mutations", "aaDeletions", "nextclade", "pangolin_clade", "status", "pangolin-note"]
    writer = csv.DictWriter(outfile, filednames, lineterminator='\n')
    writer.writeheader()
    for line in final_table:
        writer.writerow(line)
