from Bio import SeqIO
import csv
import pandas as pd
from sys import argv
import operator

"""
create variants table - for each sample in fasta multiple alignment file a covid variant is decided if found.
including pangolin and nextclades variants as well. 
"""


def calculate_coverage(fasta_seq):
    """
    calculate percentage of coverage of fasta sequence. % of no Ns # TODO: include '-'?
    :param ref_length: length of coverage sequence
    :param fasta_seq: fasta sequence to calculate coverage as SeqIO
    :return: % coverage of fasta sequence
    """
    ref_length = len(fasta_seq)
    nCount = fasta_seq.upper().count('N')
    return ((ref_length - nCount)/ref_length) * 100


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
pangolin_file = argv[3]
clades_path = argv[4]  # nextclade tsv
excel_path = argv[5]  # mutations table path
if len(argv) > 6:
    qc_report_path = argv[6]
    try:
        qc = pd.read_csv(qc_report_path, sep='\t')
        qc['sample'] = qc['sample'].apply(str)
    except FileNotFoundError:
        print("QC File does not exist")
        qc = pd.DataFrame()
else:
    qc_report_path = ''

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
    # ';' instead of ',' as sep. for Ari's tables -> important
    aa_substitution_dict[sample] = [f"{x.split(':')[1]}({x.split(':')[0]})" for x in aasubs[0].split(',')] \
        if (aasubs and aasubs != ['']) else ''
    aa_deletions_dict[sample] = [f"{x.split(':')[1]}({x.split(':')[0]})" for x in aadels[0].split(',')] \
        if (aadels and aadels != ['']) else ''
    insertions_dict[sample] = insertions[0].split(',')

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
            elif alt != ref:  # alt is not the expected mut. and is covered in sequencing (not N)
                unexpected_mutations[sample].append(mutation_name + "(alt:" + alt + ")")


unique_lineages = excel_mutTable.keys()  # get list of all unique lineages (= names of sheets of excel table)
mutations_by_lineage = {key: excel_mutTable[key].variant.tolist() for key in
                        excel_mutTable}  # mutations table with only lineage: mutation name

final_table = []
# iterate over all samples mutations and determine variants
for sample, sample_mutlist in samples_mutations.items():
    known_variant = ""
    notes = ''
    suspect_info = ''
    lin_percentages = {}
    lin_percentages_noN = {}
    lin_number = {}
    lin_number_noN = {}
    lin_not_covered = {}
    extra_subs = []
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

        # elif len(linmuts) != len(temp):  # some of lineage's mutations exist but not all
        # calculate percentage of lineage mutations found:
        lin_percentages[lin] = round(len(set(temp_mutes)) / len(set(linmuts)) * 100, 2)
        no_N_mutations = [x for x in linmuts if (x not in samples_not_covered[sample]) or (x in temp_mutes)]
        lin_percentages_noN[lin] = round(len(set(temp_mutes)) / len(set(no_N_mutations)) * 100, 2) if no_N_mutations else 0.0
        # also keep as fraction x/y
        lin_number[lin] = (len(set(temp_mutes)), len(set(linmuts)))  # tuple:
        # (#lin_mutation_sample, #tot_lin_mutations)
        lin_number_noN[lin] = (len(set(temp_mutes)), len(set(no_N_mutations))) if no_N_mutations else (0, 0)

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

    if known_variant:  # there is a variant at least X% covered
        # then list extra mutations (= mutations that exist in nextclade list but not part of the variant)
        subs_list = [x.split('(')[0] for x in aa_substitution_dict[sample]]  # TODO - CHECK
        # TODO check protein as well (not just the name)
        extra_subs = [x for x in subs_list if x not in mutations_by_lineage[known_variant]]

    if not suspect_info and (samples_not_covered[sample] or unexpected_mutations[sample]):
        # not specific suspect variant but some mutations exist \ not covered in sequencing - write as suspect
        suspect_info = 'suspect'

    # get coverage of sample from qc report.txt
    if qc_report_path:
        coverage = qc[qc['sample'] == sample]['coverageCNS_5%'].values[0].round(2)
    else:
        coverage = str(calculate_coverage(alignment[sample].seq))
        # coverage = ''

    # get pangolin info from table
    try:
        pangolin_clade = pangolinTable[pangolinTable['taxon'] == sample].lineage.values[0]
        pangolin_status = pangolinTable[pangolinTable['taxon'] == sample].status.values[0]
        pangolin_scorpio = pangolinTable[pangolinTable['taxon'] == sample].scorpio_call.values[0]
    except:
        pangolin_clade = '-'
        pangolin_status = ''
        pangolin_scorpio = ''
    QCfail = True if pangolin_status == 'fail' else False

    # get nextclade info from table
    nextclade = clades_df[clades_df['sample'] == sample].clade

    not_covered_list = []
    for mut in samples_not_covered[sample]:
        protein_values = mutTable[mutTable.variant == mut].protein.values.tolist()
        if QCfail:  # do not output list of uncovered mutations in case of QCFail. just write QCFail.
            # probably most of the mutations will be N's in that case.
            break
        if protein_values:
            not_covered_list.append(mut + "(" + protein_values[0] + ")")
        else:
            not_covered_list.append(mut + "()")

    not_covered_list = ";".join(set(not_covered_list)) if not_covered_list else ''

    line = {
        "Sample": sample,
        "Variant": pangolin_clade if not QCfail else "QC fail",
        "suspect": None,
        "suspected variant": remove_prefix(suspect_info.split(':')[0], 'suspect ') if suspect_info else '',
        "suspect info": suspect_info,  # TODO add more info
        "AA substitutions": ';'.join(aa_substitution_dict[sample]) if aa_substitution_dict and
                                                                      sample in aa_substitution_dict else 'NA',
        "AA deletions": ';'.join(aa_deletions_dict[sample] if aa_deletions_dict and sample
                                                             in aa_deletions_dict else 'NA'),
        "Insertions": ';'.join(insertions_dict[sample]) if insertions_dict and sample in insertions_dict else 'NA',
        "mutations not covered": not_covered_list,
        "non variant mutations": ';'.join(extra_subs) if extra_subs else '',  # mutations that are not part of variant list incase there is a known variant
        "% coverage": coverage,
        "pangolin clade": pangolin_clade,
        "pangolin scorpio": pangolin_scorpio,
        "nextstrain clade": nextclade.values[0] if not nextclade.empty else ''
    }
    final_table.append(line)

with open(output_file, 'w') as outfile:
    fieldnames = ["Sample", "Variant", "suspect", "suspected variant", "suspect info", "AA substitutions",
                  "AA deletions", "Insertions", "mutations not covered", "non variant mutations", "% coverage",
                  "pangolin clade", "pangolin scorpio", "nextstrain clade"]
    writer = csv.DictWriter(outfile, fieldnames, lineterminator='\n')
    writer.writeheader()
    for line in final_table:
        writer.writerow(line)

# write a report of all the variants (that are not 0%) with percentages for each sample
