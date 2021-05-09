from Bio import SeqIO
import csv
import pandas as pd
from sys import argv
"""
create variants table - for each sample in fasta multiple alignment file a covid variant is decided if found.
including pangolin and nextclades variants as well. 
"""

# get user input
alignment_file = argv[1]
output_file = argv[2]
pangolin_file = argv[3]
clades_path = argv[4]


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
    if variant not in ["B.1.1.7", "A.23.1", "P.2", "B.1.526"]:
        return unexpected_muts_dict

    new_unexpected = unexpected_muts_dict.copy()  # create shallow copy to avoid changing the original

    for x in unexpected_muts_dict[sample]:
        if "P681R" in x and variant == "B.1.1.7":
            new_unexpected[sample].remove(x)
        elif "P681H" in x and variant == "A.23.1":
            new_unexpected[sample].remove(x)
        elif "M234I" in x and variant in ["B.1.526 ", "P.2"]:
            new_unexpected[sample].remove(x)
        elif "Q677H" in x and variant == "B.1.1.7":
            new_unexpected[sample].remove(x)

    return new_unexpected


# load pangolin + nextclade outputs, mutations table.
pangolinTable = pd.read_csv(pangolin_file)
clades_df = pd.read_csv(clades_path, sep='\t')
if len(argv) > 5:
    muttable_path = argv[5]
else:
    muttable_path = "/data/projects/Dana/scripts/covid19/novelMutTable.csv"
mutTable = pd.read_csv(muttable_path)

# prepare mutations table dataframe
mutTable["REF"] = mutTable["REF"].apply(lambda x: x.upper())
mutTable["mut"] = mutTable["mut"].apply(lambda x: x.upper())
mutTable = mutTable[mutTable.type != 'Insertion']  # Ignore insertions for now

# prepare multiple alignment dictionary (key: sample name, val: SeqIO record)
alignment = SeqIO.to_dict(SeqIO.parse(alignment_file, 'fasta'))
alignment.pop('NC_045512.2', None)
alignment.pop('REF_NC_045512.2', None)

# prepare nextclade dataframe
clades_df = clades_df[['seqName', 'aaSubstitutions', 'aaDeletions', 'clade']]  # or aaDeletions instead?? # TODO
clades_df = clades_df.rename(columns={'seqName': 'sample'})
clades_df['sample'] = clades_df['sample'].apply(str)
clades_df = clades_df.fillna('')
# create dict of aasubstitutions and aadeletions. key=sample id, value=list of substitutions..
aa_substitution_dict = {}
aa_deletions_dict = {}
for sample in clades_df['sample']: # for each sample aggregate
    aasubs = clades_df[clades_df['sample'] == sample].aaSubstitutions.values.tolist()
    aadels = clades_df[clades_df['sample'] == sample].aaDeletions.values.tolist()
    aa_substitution_dict[sample] = [f"{x.split(':')[1]}({x.split(':')[0]})" for x in aasubs[0].split(',')] \
        if (aasubs and aasubs != ['']) else ''
    aa_deletions_dict[sample] = [f"{x.split(':')[1]}({x.split(':')[0]})" for x in aadels[0].split(',')] \
        if (aadels and aadels != ['']) else ''

# some variables:
samples_mutations = {id: [] for id in alignment}
# samples_s_not_covered = {id: [] for id in alignment}
samples_not_covered = {id: [] for id in alignment}
unexpected_mutations = {id: [] for id in alignment}
lineages_list = []

# iterate over all samples in multifasta and over all mutations in table, and check value of each mutation
for sample, record in alignment.items():
    for row in mutTable.iterrows():
        pos = int(row[1][5])-1  # mutation position
        alt = record.seq[pos]  # fasta value in position
        ref = row[1][6]  # reference in position
        table_mut = row[1][7]  # mutation  according to table
        gene = row[1][2] if row[1][2] else ''
        lineage = row[1][3]
        lineages_list += [x.strip() for x in lineage.split(',')]  # accumulate all lineages to list - for later
        mutation_name = row[1][1] if row[1][1] else row[1][0]  # if no AA name take nucleotide name
        if alt != table_mut or alt == 'N':
            if alt == 'N':  # and gene == 'S':
                # samples_s_not_covered[sample].append(mutation_name)
                samples_not_covered[sample].append(mutation_name)
            elif alt != ref and alt != 'N':  # alt is not the expected mut. and is covered in sequencing
                unexpected_mutations[sample].append(mutation_name + "(alt:" + alt + ")")
        else:  # some variant mutation
            samples_mutations[sample].append(mutation_name)  # accumulate all samples+mutations in dict:
            # {sampleName: mutationName}

unique_lineages = set(lineages_list)
# create a dictionary of mutations by lineage -> {key=lineage: val=dataframe with all of lineage's mutations}
mutations_by_lineage = {x: mutTable[mutTable.lineage.str.contains(x, regex=False)].AA.tolist() for x in unique_lineages}

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
        temp = [x for x in linmuts]
        temp_mutes = []
        for mut in sample_mutlist:  # iterate over  table mutations of sample
            if mut in temp:  # remove mutations of sample from temp list and keep in temp_mutes list
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
            if val >= 60:  # 60% -> knwon variant
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
            suspect = f'suspect_{var}: {str(lin_percentages[var])}{fraction}'  # suspect_<lineagename>): (%)(x/y)
    unexpected_mutations = specific_cases(unexpected_mutations, sample, known_variant)

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
    line = {
        "Sample": sample,
        "Known Variant": known_variant if known_variant and not QCfail else 'no monitored variant',
        "Suspect": suspect,
        "More Mutations": ';'.join(set([x + "(" + mutTable[mutTable.AA == x].gene.values[0] + ")" for x in more_muts])),
        # "S Not Covered": ';'.join(samples_s_not_covered[sample]),
        "Not Covered": ';'.join(set([x + "(" + mutTable[mutTable.AA == x].gene.values[0] +
                                     ")" for x in samples_not_covered[sample]])) if not QCfail else '',
        # "non-Table Mutations": ';'.join(unexpected_mutations[sample]),
        "all mutations": ';'.join(aa_substitution_dict[sample]) if aa_substitution_dict else '',
        "nextclade": nextclade.values[0] if not nextclade.empty else '',
        "pangolin_clade": pangolin_clade,
        "status": pangolin_status,
        "pangolin-note": pangolin_note
    }
    final_table.append(line)

with open(output_file, 'w') as outfile:
    filednames = ["Sample", "Known Variant", "Suspect", "More Mutations", "Not Covered",  # 'S Not Covered'
                  "all mutations", "nextclade", "pangolin_clade", "status", "pangolin-note"]
    writer = csv.DictWriter(outfile, filednames, lineterminator='\n')
    writer.writeheader()
    for line in final_table:
        writer.writerow(line)

