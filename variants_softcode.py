from Bio import SeqIO
import csv
import pandas as pd
from sys import argv

# TODO: specify mutations that share location: ignore and do not add to 'non table mutations'!
# TODO: create more functions based code, with main
alignment_file = argv[1]
output_file = argv[2]
pangolin_file = argv[3]

pangolinTable = pd.read_csv(pangolin_file)
# mutTable = pd.read_csv("novelMutTable.csv")
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


def specific_cases(unexpected_muts_dict, sample, variant):
    """
    remove double mutations from unexpected mutations dictionary.
    double mutations:
        1. P681H(UK) P981R(Uganda) pos 23604
        2. M234I G-T(Rio), G-A(NY) pos 28975

    :param more_muts_list: extra mutations, more than variant
    :param unexpected_muts_list: mutations that are different from ref but also from mutations table
    :param variant: which variant already chosen
    :return: more_muts_list and unexpected_muts_list updated.
    """
    if variant not in ["B.1.1.7 - UK", "A.23.1 Uganda", "P.2- Rio de jeneiro", "B.1.526 New york",
                       "VOI-18.02-WHO"]:
        return unexpected_muts_dict

    new_unexpected = unexpected_muts_dict.copy()  # create shallow copy to avoid changing the original

    for x in unexpected_muts_dict[sample]:
        if "P681R" in x and variant == "B.1.1.7 - UK":
            new_unexpected[sample].remove(x)
        elif "P681H" in x and variant == "A.23.1 Uganda":
            new_unexpected[sample].remove(x)
        elif "M234I" in x and variant in ["B.1.526 New york", "P.2- Rio de jeneiro"]:
            new_unexpected[sample].reomve(x)
        elif "Q677H" in x and variant in ["VOI-18.02-WHO", "B.1.1.7 - UK"]:
            new_unexpected[sample].remove(x)

    return new_unexpected


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
    lin_number = {}
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
            lin_percentages[lin] = round(len(temp_mutes) / len(linmuts) * 100, 2)
            lin_number[lin] = (len(temp_mutes), len(linmuts))  # tuple: (#mutation_sample, #tot_lin_mutations)

    if known_variant:
        more_muts = [x for x in more_muts if x not in mutations_by_lineage[known_variant]]
    else:
        max=0
        var = ''
        flag = False
        val = 0
        # for key, val in lin_percentages.items():
        #     if val >= 70:
        #         flag = True
        #     if val > max:
        #         max = val
        #         var = key
        for key, tup in lin_number.items():
            val = round((tup[0] / tup[1])*100, 2)
            if val >= 70:
                flag = True
            if val > max:
                max = val
                var = key
        if flag:
            known_variant = var
            more_muts = [x for x in more_muts if x not in mutations_by_lineage[known_variant]]
        elif val >= 50 and var:  # and < 70
            more_muts = [x for x in more_muts if x not in mutations_by_lineage[var]]

        if var and lin_number[var][0] >= 2:
            suspect = 'suspect_' + var + ": " + str(lin_percentages[var]) + "%"
    unexpected_mutations = specific_cases(unexpected_mutations, sample, known_variant)

    more_muts = set(more_muts)
    if not suspect and (more_muts or samples_s_not_covered[sample] or unexpected_mutations[sample]):
        suspect = 'suspect'

    line = {
        "Sample": sample,
        "Known Variant": known_variant if known_variant else 'no variant',
        "Suspect": suspect,
        "More Mutations": ';'.join(set([x + "(" + mutTable[mutTable.AA == x].gene.values[0] + ")" for x in more_muts])),
        "S Not Covered": ';'.join(samples_s_not_covered[sample]),
        "non-Table Mutations": ';'.join(unexpected_mutations[sample]),
        "pangolin_clade": pangolinTable[pangolinTable.taxon == sample].lineage.values[0],
        "status": pangolinTable[pangolinTable.taxon == sample].status.values[0],
        "pangolin-note": pangolinTable[pangolinTable.taxon == sample].note.values[0]
    }
    final_table.append(line)

with open(output_file, 'w') as outfile:
    filednames = ["Sample", "Known Variant", "Suspect", "More Mutations", "S Not Covered",
                  "non-Table Mutations", "pangolin_clade", "status", "pangolin-note"]
    writer = csv.DictWriter(outfile, filednames, lineterminator='\n')
    writer.writeheader()
    for line in final_table:
        writer.writerow(line)

