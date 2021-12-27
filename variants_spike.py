from Bio import SeqIO
import csv
import pandas as pd
from sys import argv
from pathlib import Path
from scipy.stats import binom
import operator
from datetime import datetime

"""
create variants table - for each sample in fasta multiple alignment file a covid variant is decided if found.
including pangolin and nextclades variants as well. 
"""


def format_rows(var):
    fraction = f'({int(lin_number[var][0])}/{int(lin_number[var][1])})'
    fraction_noN = f'({int(lin_number_noN[var][0])}/{int(lin_number_noN[var][1])})'
    row = f'suspect {var}: {str(lin_percentages[var])}% {fraction};  ' \
          f'noN: {lin_percentages_noN[var]}% {fraction_noN} ;' \
          f'p={lin_binomial[var]}'
    return row


def sort_variants_per_sample(sortby, ascending=False):
    s = pd.Series(sortby)
    s = s.sort_values(ascending=ascending)
    s = s.index.map(format_rows).tolist()
    return s


def format_line(variant):
    suspect_info = f"{variant}: {lin_percentages[variant]}% " \
                   f"({lin_number[variant][0]}/{lin_number[variant][1]}); " \
                   f"noN: {lin_percentages_noN[variant]}% ({lin_number_noN[variant][0]}/" \
                   f"{lin_number_noN[variant][1]})"
    return suspect_info


def determine_variant():
    max_count = max(list(map(lambda x: x[0], lin_number_noN.values())))
    candidates = [var for var, count in lin_number_noN.items() if count[0] >= max_count - 1]
    candidates_freq = {candidate: lin_percentages_noN[candidate] for candidate in candidates}
    variant = max(candidates_freq.items(), key=operator.itemgetter(1))[0]
    suspect_info = format_line(variant)
    return suspect_info


def calculate_coverage(fasta_seq):
    """
    calculate percentage of coverage of fasta sequence. % of no Ns # TODO: include '-'?
    :param ref_length: length of coverage sequence
    :param fasta_seq: fasta sequence to calculate coverage as SeqIO
    :return: % coverage of fasta sequence
    """
    ref_length = len(fasta_seq)
    nCount = fasta_seq.upper().count('N')
    return ((ref_length - nCount) / ref_length) * 100


def remove_prefix(text, prefix):
    """
    remove prefix from string
    :param text: string to remove prefix from
    :param prefix: prefix to remove
    :return: new string wihtout prefix
    """
    if text.startswith(prefix):
        return text[len(prefix):]
    return text


# get user input
alignment_file = argv[1]
output_file = argv[2]
excel_path = argv[3]  # mutations table path

red_flags_path = excel_path.replace('mutationsTable.xlsx', 'red_flags.csv')
output_path = Path(output_file).parent
out_fname = datetime.now().strftime('%Y%m%d') + '_variants.csv'
output_file = output_path / out_fname

red_flags_df = pd.read_csv(red_flags_path)

if len(argv) > 4:
    qc_report_path = argv[4]
    try:
        qc = pd.read_csv(qc_report_path, sep='\t')
        qc['sample'] = qc['sample'].apply(str)
    except FileNotFoundError:
        print("QC File does not exist")
        qc = pd.DataFrame()
else:
    qc_report_path = ''


# excel muttable:
excel_mutTable = pd.read_excel(excel_path, sheet_name=None, engine='openpyxl')
for name in excel_mutTable:
    excel_mutTable[name] = excel_mutTable[name][excel_mutTable[name]['annotation'] == 'Spike']
    excel_mutTable[name]['lineage'] = name  # add a lineage column to all variant's tables
    excel_mutTable[name]['variant'] = excel_mutTable[name]['variant'].astype(str)

mutTable = pd.concat(excel_mutTable.values(), ignore_index=True)
mutTable['Mutation type'] = mutTable['Mutation type'].str.lower()
mutTable = mutTable[mutTable['Mutation type'] != 'insertion']
# prepare mutations table dataframe
mutTable["Reference"] = mutTable["Reference"].str.upper()  # make sure all upper case for safe comparisons
mutTable["Mutation"] = mutTable["Mutation"].str.upper()
mutTable = mutTable.dropna(thresh=3)

# prepare multiple alignment dictionary (key: sample name, val: SeqIO record)
alignment = SeqIO.to_dict(SeqIO.parse(alignment_file, 'fasta'))
alignment.pop('Spike', None)




samples_mutations = {id: [] for id in alignment}
# samples_s_not_covered = {id: [] for id in alignment}
samples_not_covered = {id: [] for id in alignment}
unexpected_mutations = {id: [] for id in alignment}
lineages_list = []

# iterate over all samples in multi-fasta and over all mutations in table, and check value of each mutation
with open("mutations.log", 'w') as log:
    for sample, record in alignment.items():
        for (idx, row) in mutTable.iterrows():
            if pd.isna(row.loc['Position']):
                log.write(f"NaN: {row}")
                continue
            pos = int(row.loc['Position'])   # mutation position
            alt = record.seq[pos - 21563 ] # fasta value in position
            ref = row.loc['Reference']  # reference in position
            table_mut = row.loc['Mutation']  # mutation  according to table
            gene = row.loc['protein']
            mutation_name = str(row.loc['variant'])
            if alt == table_mut:  # mutation exists in sequence
                samples_mutations[sample].append(mutation_name)
            elif alt == 'N':  # if position not covered in sequence
                samples_not_covered[sample].append(mutation_name)
            elif alt != ref:  # alt is not the expected mut. and is covered in sequencing (not N)
                unexpected_mutations[sample].append(mutation_name + "(alt:" + alt + ")")

unique_lineages = excel_mutTable.keys()  # get list of all unique lineages (= names of sheets of excel table)
mutations_by_lineage = {key: excel_mutTable[key].variant.tolist() for key in
                        excel_mutTable}  # mutations table with only lineage: mutation name

final_table = []
ranked_variants_dict = dict()
# iterate over all samples mutations and determine variants
n = len(samples_mutations)
c = 0
for sample, sample_mutlist in samples_mutations.items():
    c += 1
    # print(c/n)
    known_variant = ""
    notes = ''
    suspect_info = ''
    lin_percentages = {}
    lin_percentages_noN = {}
    lin_number = {}
    lin_number_noN = {}
    lin_not_covered = {}
    lin_binomial = {}
    # iterate over mutations of each lineage, check how much of the lineage's mutations are covered by sample -
    # to decide on variant
    for lin, linmuts in mutations_by_lineage.items():
        if linmuts:
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

            # elif len(linmuts) != len(temp):  # some of lineage's mutations exist but not all
            # calculate percentage of lineage mutations found:
            lin_percentages[lin] = round(len(set(temp_mutes)) / len(set(linmuts)) * 100, 2)
            no_N_mutations = [x for x in linmuts if (x not in samples_not_covered[sample]) or (x in temp_mutes)]
            lin_percentages_noN[lin] = round(len(set(temp_mutes)) / len(set(no_N_mutations)) * 100,
                                             2) if no_N_mutations else 0.0
            # also keep as fraction x/y
            lin_number[lin] = (len(set(temp_mutes)), len(set(linmuts)))  # tuple:
            #keep covered mutations
            # (#lin_mutation_sample, #tot_lin_mutations)
            lin_number_noN[lin] = (len(set(temp_mutes)), len(set(no_N_mutations))) if no_N_mutations else (0, 0)
            lin_binomial[lin] = binom.pmf(
                k=lin_number_noN[lin][0],
                n=lin_number_noN[lin][1],
                p=0.25)

    ranked_variants_dict[sample] = sort_variants_per_sample(lin_binomial, ascending=True)

    if not known_variant:  # did not find variant that has 100% mutations in sample
        var = max(lin_percentages_noN.items(), key=operator.itemgetter(1))[0]

        val = lin_percentages[var]
        val_noN = lin_percentages_noN[var]
        fraction = f'({int(lin_number[var][0])}/{int(lin_number[var][1])})'
        fraction_noN = f'({int(lin_number_noN[var][0])}/{int(lin_number_noN[var][1])})'

        over_60 = True if val >= 60 else False
        known_variant = var if over_60 else known_variant

        if var and lin_number[var][0] >= 2:  # At least 2 mutations of lineage --> suspect variant
            # suspect <lineage>: (%)(x/y); noN(x/#covered_mutations); SNP(#); SNP_silent(#) :
            suspect_info = \
                f'suspect {var}: {str(lin_percentages[var])}% {fraction};  ' \
                f'noN: {val_noN}% {fraction_noN}'

    else:  # variant has 100% mutations
        suspect_info = f"{known_variant}: {lin_percentages[known_variant]}% " \
                       f"({lin_number[known_variant][0]}/{lin_number[known_variant][1]}); " \
                       f"noN: {lin_percentages_noN[known_variant]}% ({lin_number_noN[known_variant][0]}/" \
                       f"{lin_number_noN[known_variant][1]})"


    if not suspect_info and (samples_not_covered[sample] or unexpected_mutations[sample]):
        # not specific suspect variant but some mutations exist \ not covered in sequencing - write as suspect
        suspect_info = 'suspect'

    suspect_info = determine_variant()
    sus_variant_name = suspect_info.split(':')[0] if suspect_info else ''



    # get coverage of sample from qc report.txt
    if qc_report_path:
        coverage = qc[qc['sample'] == sample]['coverageCNS_5%'].values[0].round(2)
    else:
        coverage = str(calculate_coverage(alignment[sample].seq))
        # coverage = ''

    # add protein value
    # not covered list = all the mutations in the bodek that the sequence did not covered (N's)
    not_covered_list = []
    for mut in samples_not_covered[sample]:
        protein_values = mutTable[mutTable.variant == mut].protein.values.tolist()

        if protein_values:
            not_covered_list.append(mut + "(" + protein_values[0] + ")")
        else:
            not_covered_list.append(mut + "()")

    not_covered_list = ";".join(set(not_covered_list)) if not_covered_list else ''

    # add protein value
    # covered list = all the mutations in the bodek that the were found in the sequence
    covered_list = []
    for mut in samples_mutations[sample]:
        protein_values = mutTable[mutTable.variant == mut].protein.values.tolist()
        if protein_values:
            covered_list.append(mut + "(" + protein_values[0] + ")")
        else:
            covered_list.append(mut + "()")

    covered_list = ";".join(set(covered_list)) if covered_list else ''

    #not covered in the suspected variant
    not_covered_in_suspect = [value for value in mutations_by_lineage[sus_variant_name] if value.split("(")[0] in not_covered_list]
    not_covered_in_suspect = ";".join(set(not_covered_in_suspect)) if not_covered_in_suspect else ''

    #covered in the suspected variant
    covered_in_suspect = [value for value in mutations_by_lineage[sus_variant_name] if value.split("(")[0] in covered_list]
    covered_in_suspect = ";".join(set(covered_in_suspect)) if covered_in_suspect else ''

    line = {
        "Sample": sample,
        "suspected variant": sus_variant_name,
        "suspect info": suspect_info,  # TODO add more info
        "covered mutations": covered_list,
        "mutations not covered": not_covered_list,
        "not covered in suspect": not_covered_in_suspect,
        "covered in suspect": covered_in_suspect,
        # mutations that are not part of variant list incase there is a known variant
        "% coverage": coverage,

    }
    final_table.append(line)

with open(output_file, 'w') as outfile:
    fieldnames = ["Sample", "suspected variant", "suspect info","covered mutations","mutations not covered", "non variant mutations", "not covered in suspect", "covered in suspect" ,"% coverage"]

    writer = csv.DictWriter(outfile, fieldnames, lineterminator='\n')
    writer.writeheader()
    for line in final_table:
        writer.writerow(line)

ranked_variants_file = output_path / 'ranked_variants.csv'
ranked_variants_df = pd.DataFrame(ranked_variants_dict)
ranked_variants_df.to_csv(ranked_variants_file)
# write a report of all the variants (that are not 0%) with percentages for each sample
