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

#sort the variants by the noN % field. (descending)
def sort_variants_per_sample(sortby, ascending=False):
    s = pd.Series(sortby)
    s = s.sort_values(ascending=False)
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
    candidates = [var for var, count in lin_number_noN.items() if count[0] >= max_count-1]
    candidates_freq = {candidate : lin_percentages_noN[candidate] for candidate in candidates}
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

def recombinant_suspect(nt_substitutions_list, suspect_info):
    if not suspect_info.startswith(("BA.1", "BA.2", "B.1.617.2 signature d")) : 
        return ""
    else:
        ba1_muts = [value for value in unique_muts[unique_muts["lineage"] == "BA.1"].substitution.tolist() if value in nt_substitutions_list]
        ba2_muts = [value for value in unique_muts[unique_muts["lineage"] == "BA.2"].substitution.tolist() if value in nt_substitutions_list]
        deltaD_muts = [value for value in unique_muts[unique_muts["lineage"] == "B.1.617.2 signature d"].substitution.tolist() if value in nt_substitutions_list]
        
        suspect_rec_muts = unique_muts[unique_muts.substitution.isin([*ba1_muts, *ba2_muts ,*deltaD_muts])].sort_values(by=['position']).reset_index()
        swap = 0 #count the swaps between variants
        mut_counter = 0 #count the mutations of the current variant
        brkpnt =[] #list of swaps position
        for index, row in suspect_rec_muts.iterrows():
            if index == 0:
                var = row['lineage']
            mut_counter += 1
            temp = var
            var = row['lineage']
            if var != temp:
                swap += 1
                brkpnt.append(str(suspect_rec_muts.iloc[index-1]["position"]) + "-" + str(suspect_rec_muts.iloc[index]["position"]))  #position of variant swap

        rec = {
            "Sample": sample,
            "swaps": swap,
            "BA.1 mutations": ';'.join(ba1_muts),
            "BA.2 mutations": ';'.join(ba2_muts),
            "deltaD mutations": ';'.join(deltaD_muts),
            "breakpoints": ';'.join(brkpnt)
            }
        recombinants.append(rec)
        if swap == 1 and ( (len(ba1_muts)>2 and len(ba2_muts)>2) or (len(ba1_muts)>2 and len(deltaD_muts)>2) or (len(deltaD_muts)>2 and len(ba2_muts)>2) ):
            return "X"
        return ""

DEBUG = False
if DEBUG:
    debug_path = Path('/home/omera/Code/sandbox/ar/128')
    alignment_file = debug_path / 'alignment' / 'all_aligned.fasta'
    output_file = debug_path / 'results' / 'variants_debug.csv'
    pangolin_file = debug_path / 'results' / 'pangolinClades.csv'
    clades_path = debug_path / 'results' / 'nextclade.tsv'
    excel_path = '/home/omera/Code/covid19/mutationsTable.xlsx'
    
else:
    # get user input
    alignment_file = argv[1]
    output_file = argv[2]
    pangolin_file = argv[3]
    clades_path = argv[4]  # nextclade tsv
    excel_path = argv[5]  # mutations table path
    
    
red_flags_path = excel_path.replace('mutationsTable.xlsx', 'red_flags.csv')
output_path = Path(output_file).parent
out_fname = datetime.now().strftime('%Y%m%d') + '_variants.csv'
output_file = output_path / out_fname

red_flags_df = pd.read_csv(red_flags_path)

unique_muts = pd.read_excel(excel_path.replace("mutationsTable.xlsx","unique_BA.1_BA.2_DeltaD_muts.xlsx"), engine='openpyxl')


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
mutTable = pd.concat(excel_mutTable.values(), ignore_index=True)
mutTable['mutation type'] = mutTable['mutation type'].str.lower()
mutTable = mutTable[mutTable['mutation type'] != 'insertion']
# prepare mutations table dataframe
mutTable["reference"] = mutTable["reference"].str.upper()  # make sure all upper case for safe comparisons
mutTable["mutation"] = mutTable["mutation"].str.upper()
mutTable = mutTable.dropna(thresh=3)

# prepare multiple alignment dictionary (key: sample name, val: SeqIO record)
alignment = SeqIO.to_dict(SeqIO.parse(alignment_file, 'fasta'))
alignment.pop('NC_045512.2', None)
alignment.pop('REF_NC_045512.2', None)

# prepare nextclade dataframe
clades_df = clades_df[['seqName', 'aaSubstitutions', 'aaDeletions', 'clade', 'insertions', 'substitutions']]
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

samples_mutations = {id: [] for id in aa_substitution_dict if id in alignment}
# samples_s_not_covered = {id: [] for id in alignment}
samples_not_covered = {id: [] for id in aa_substitution_dict if id in alignment}
unexpected_mutations = {id: [] for id in aa_substitution_dict if id in alignment}
lineages_list = []

# iterate over all samples in multi-fasta and over all mutations in table, and check value of each mutation
with open("mutations.log", 'w') as log:
    for sample, record in alignment.items():
        if sample in aa_substitution_dict.keys():
            for (idx, row) in mutTable.iterrows():
                if pd.isna(row.loc['position']):
                    log.write(f"NaN: {row}")
                    continue
                pos = int(row.loc['position']) - 1  # mutation position
                alt = record.seq[pos]  # fasta value in position
                ref = row.loc['reference']  # reference in position
                table_mut = row.loc['mutation']  # mutation  according to table
                gene = row.loc['protein']
                mutation_name = str(row.loc['variant'])
                if alt == table_mut:  # mutation exists in sequence
                    samples_mutations[sample].append(mutation_name)
                elif alt == 'N':  # if position not covered in sequence
                    samples_not_covered[sample].append(mutation_name)
                elif alt != ref:  # alt is not the expected mut. and is covered in sequencing (not N)
                    unexpected_mutations[sample].append(mutation_name + "(alt:" + alt + ")")


mutations_by_lineage = mutTable.groupby('variantname')['variant'].apply(list).to_dict()

final_table = []
recombinants = []
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
        lin_binomial[lin] = binom.pmf(
            k=lin_number_noN[lin][0],
            n=lin_number_noN[lin][1],
            p=0.25)
    
    ranked_variants_dict[sample] = sort_variants_per_sample(lin_number_noN, ascending=True)

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

    #if known_variant:  # there is a variant at least X% covered
     #   non_variant_mut = [value for value in aa_substitution_dict[sample] if value.split("(")[0] not in mutations_by_lineage[known_variant]]


    if not suspect_info and (samples_not_covered[sample] or unexpected_mutations[sample]):
        # not specific suspect variant but some mutations exist \ not covered in sequencing - write as suspect
        suspect_info = 'suspect'


    suspect_info = determine_variant()
    sus_variant_name = suspect_info.split(': ')[0] if suspect_info else ''


    # get coverage of sample from qc report.txt
    if qc_report_path:
        coverage = qc[qc['sample'] == sample]['coverageCNS_5%'].values[0].round(2)
    else:
        coverage = str(calculate_coverage(alignment[sample].seq))

        #coverage = ''

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

    #add protein value
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
    
    
    nt_substitutions = clades_df.loc[clades_df['sample'].eq(sample), 'substitutions'].str.split(',')
    nt_substitutions_list = ';'.join(nt_substitutions.values[0]).split(";")
    
    non_variant_mut = [value for value in aa_substitution_dict[sample] if value.split("(")[0] not in mutations_by_lineage[sus_variant_name]]
    non_variant_mut = [value for value in non_variant_mut if value.split("(")[1] not in ["ORF1a)","ORF1b)"]] #temporery - until the bodek's numecleture will match nextclade's
    non_variant_mut = ";".join(set(non_variant_mut))
    

    if nt_substitutions_list:
        red_flags = red_flags_df.loc[red_flags_df['SNP'].isin(nt_substitutions), 'SNP']
        red_flags_str = ';'.join(red_flags)
        nt_substitutions_str = ';'.join(nt_substitutions_list)
    else:
        red_flags_str = nt_substitutions_str = ''

    line = {
        "Sample": sample,
        "Variant": pangolin_clade if not QCfail else "QC fail",
        "suspect": None,
        "suspected variant": remove_prefix(suspect_info.split(':')[0], 'suspect').lstrip(' ') if int(suspect_info.split('noN:')[1].split(".")[0]) > 60 else 'QC fail', #if the suspected variant < 60% => QC fail
        "suspect info": suspect_info,  # TODO add more info
        'nt substitutions' : nt_substitutions_str,
        'red_flags' : red_flags_str,
        "AA substitutions": ';'.join(aa_substitution_dict[sample]) if aa_substitution_dict and
                                                                      sample in aa_substitution_dict else 'NA',
        "AA deletions": ';'.join(aa_deletions_dict[sample] if aa_deletions_dict and sample
                                                             in aa_deletions_dict else 'NA'),
        "Insertions": ';'.join(insertions_dict[sample]) if insertions_dict and sample in insertions_dict else 'NA',
        "mutations not covered": not_covered_list,
        "non variant mutations": non_variant_mut,  # mutations that are not part of variant list incase there is a known variant
        "% coverage": coverage,
        "pangolin clade": pangolin_clade,
        "pangolin scorpio": pangolin_scorpio,
        "nextstrain clade": nextclade.values[0] if not nextclade.empty else '',
        "recombinant suspect": recombinant_suspect(nt_substitutions_list, suspect_info) 

    }
    final_table.append(line)

with open(output_file, 'w') as outfile:
    fieldnames = ["Sample", "Variant", "suspect", "suspected variant", "suspect info", "AA substitutions",
                  "AA deletions", "Insertions", "mutations not covered", "non variant mutations", "% coverage",
                  "pangolin clade", "pangolin scorpio", "nextstrain clade", 'nt substitutions', 'red_flags',"recombinant suspect"]
    writer = csv.DictWriter(outfile, fieldnames, lineterminator='\n')
    writer.writeheader()
    for line in final_table:
        writer.writerow(line)


#append low quel
    low_quel = {id: [] for id in aa_substitution_dict if id not in alignment}
    for sample in low_quel:
        if qc_report_path:
            coverage = qc[qc['sample'] == sample]['coverageCNS_5%'].values[0].round(2)
        else:
            coverage = ""
        line = {
        "Sample": sample,
        "% coverage": coverage,
        "Variant": "QC fail",
        "suspected variant": "QC fail"
        }
        writer.writerow(line)


ranked_variants_file = output_path / 'ranked_variants.csv'
ranked_variants_df = pd.DataFrame(ranked_variants_dict)
ranked_variants_df.to_csv(ranked_variants_file)
# write a report of all the variants (that are not 0%) with percentages for each sample


#recombinants file
with open("results/recombinants.csv", 'w') as recombinant_file:
    names = [ "Sample", "swaps","BA.1 mutations", "BA.2 mutations", "deltaD mutations","breakpoints"]
    writer = csv.DictWriter(recombinant_file, names, lineterminator='\n')
    writer.writeheader()
    for rec in recombinants:
        writer.writerow(rec)
    