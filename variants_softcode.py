from Bio import SeqIO
import csv
import pandas as pd
from sys import argv
import operator

# TODO: Dictionary of dictionaries. maybe create mutations database more convenient
# Mayby create mutations as json??? and list all lineages in the same row. eah mutation can be a class instant that will
# TODO: 70% of certain mutation is enough but only if did not fully match other.
# mark partial match, sort it, and if top is > 70% choose it
# no hard code!

# Dict of lineages:

alignment_file = argv[1]
output_file = argv[2]
# mutTable = pd.read_csv("novelMutTable_2.csv")
mutTable = pd.read_csv("/data/projects/Dana/scripts/covid19/novelMutTable_2.csv") # TODO: change to novelMutTable, change table itself

mutTable["REF"] = mutTable["REF"].apply(lambda x: x.upper())
mutTable["mut"] = mutTable["mut"].apply(lambda x: x.upper())
mutTable = mutTable[mutTable.type != 'Insertion'] # Ignore insertions for now


alignment = SeqIO.to_dict(SeqIO.parse(alignment_file, 'fasta'))  # fasta-dict -> {id: seq object}
alignment.pop('NC_045512.2', None)
alignment.pop('REF_NC_045512.2', None)


samples_mutations = {id: [] for id in alignment}
samples_s_not_covered = {id: [] for id in alignment}
unexpected_mutations = {id: [] for id in alignment}
lineages_list = []
for sample, record in alignment.items():
    for row in mutTable.iterrows():
        pos = int(row[1][5])-1
        alt = record.seq[pos]
        ref = row[1][6]
        table_mut = row[1][7]
        gene = row[1][2] if row[1][2] else ''
        lineage = row[1][3]
        lineages_list += [x.strip() for x in lineage.split(',')]
        mutation_name = row[1][1] if row[1][1] else row[1][0]  # if no AA name take nucleotide name
        if alt != table_mut or alt == 'N':
            if alt == 'N' and gene == 'S':
                samples_s_not_covered[sample].append(mutation_name)
            elif alt != ref and alt != 'N':  # alt is not the expected mut. and is covered in sequencing
                unexpected_mutations[sample].append(mutation_name + "(alt:" + alt + ")")
        else:  # some variant mutation
            samples_mutations[sample].append(mutation_name)


unique_lineages = set(lineages_list)

mutations_by_lineage = {x: mutTable[mutTable.lineage.str.contains(x)].AA.tolist() for x in unique_lineages}

final_table = []

for sample, sample_mutlist in samples_mutations.items():
    known_variant = ""
    more_muts = []
    notes = ''
    suspect = ''
    lin_percentages = {}
    for lin, linmuts in mutations_by_lineage.items():
        temp = [x for x in linmuts]
        temp_mutes = []
        for mut in sample_mutlist:
            if mut in temp:
                temp.remove(mut)
                temp_mutes.append(mut)
        if not temp:
            # all mutations of that lineage exists.
            known_variant = lin
        elif len(linmuts) != len(temp):  # some mutations do exist
            more_muts += temp_mutes
            lin_percentages[lin] = len(temp_mutes) / len(linmuts) * 100

    if known_variant:
        more_muts = [x for x in more_muts if x not in mutations_by_lineage[known_variant]]
    else:
        max=0
        var = ''
        flag = False
        for key, val in lin_percentages.items():
            if val >= 70:
                flag = True
            if val > max:
                max = val
                var = key
        if flag:
            known_variant = var
            suspect = 'suspect_' + var + ": " + str(lin_percentages[var]) + "%"

    more_muts = set(more_muts)
    if not suspect and (more_muts or samples_s_not_covered[sample] or unexpected_mutations[sample]):
        suspect = 'suspect'

    line = {
        "Sample": sample,
        "Known Variant": known_variant if known_variant else 'no variant',
        "Suspect": suspect if more_muts or samples_s_not_covered[sample] or unexpected_mutations[sample] else '',
        "More Mutations": ';'.join(set([x + "(" + mutTable[mutTable.AA == x].gene.values[0] + ")" for x in more_muts])),
        "S Not Covered": ';'.join(samples_s_not_covered[sample]),
        "non-Table Mutations": ';'.join(unexpected_mutations[sample])
    }
    final_table.append(line)

with open(output_file, 'w') as outfile:
    filednames = ["Sample", "Known Variant", "Suspect", "More Mutations", "S Not Covered",
                  "non-Table Mutations"]
    writer = csv.DictWriter(outfile, filednames, lineterminator='\n')
    writer.writeheader()
    for line in final_table:
        writer.writerow(line)