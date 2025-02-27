import numpy as np
from Bio import SeqIO
import pandas as pd
from sys import argv
from pathlib import Path
from datetime import datetime
# import resist

import time
start = time.time()



def format_line(info):
    '''
    format the suspected variant info for suspect info column in variants.csv.
    :return: example - BA.2: 73.85% (48/65); noN: 78.69% (48/61)
    '''
    return info["lineage"] + ": "  + str(info["mutation_precentage"]) + "% " + info["mutations_fraction"] + \
    "; noN: " + str(info["noN_mutation_precentage"]) + "% " + info["noN_mutations_fraction"]


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


def extra_mutations(sus_mutations, bodek_mutations, aa=None):
    '''
    find mutations that the samples had but are not in the suspected variant.
    works for amino acid and for nucleotide(default) mutations
    '''
    if aa:
        return [value for value in sus_mutations if value.split("(")[0] not in bodek_mutations 
                and value.split("(")[1] not in ["ORF1a)","ORF1b)"]] #temporery - until the bodek's numecleture will match nextclade's
    return [value for value in sus_mutations if value not in bodek_mutations]

def set(l):
    '''
    @overrides
    '''
    return list(dict.fromkeys(l))

def var_mut_counter(muts):
    '''
    get a list of the variant mutations and calculate unique mutations.
    this function was written to count a consecutive deletions as a single event.
    '''

    # del_lists is a list that will contain lists of deletions
    # each list will contain mutations of one deletion event.
    del_lists = []
    all_deletions = [int(x[1:-1]) for x in muts if x.endswith("-")]
    all_deletions.sort()
    dels = []
    prev_mut = 0 # initiate
    for mut in all_deletions:
        cur_mut = int(mut)#int(mut[1:-1])  # get only the Position of the mutation.
        if not cur_mut == prev_mut + 1:  # meaning its not the same deletion event
            del_lists.append(dels) if not prev_mut == 0 else None  # dont append the first list - its empty
            dels = []  # reset list - new deletion event
            dels.append(mut)  # add first mutation in the new event
        else:
            dels.append(mut)
        prev_mut = cur_mut
    del_lists.append(dels)
    muts_sum = len(muts) - len(all_deletions) + len(del_lists)
    return del_lists, muts_sum

def sample_mut_counter(muts,del_list):
    '''
    get a list of the samples' mutations and variant deletions, returns unique mutations.
    this function was written to count a consecutive deletions as a single event.
    :param muts: a list of the samples mutations.
    :param del_list: a list of a variant deletions lists. each deletions list contains mutations of one deletion event.
    :returns mut_sum: unique mutations calculation of the intersection between the variants' and samples' mutations.
    '''
    del_count = 0
    samples_del = [x for x in muts if x.endswith("-")]
    for del_event in del_list:  # itetates variants' deletions event
        shared_del = [x for x in del_event if x in samples_del]  # variant and sample mutations intersection.
        if len(shared_del) > 0:
            del_count += 1
    muts_sum = len(muts) - len(samples_del) + del_count
    return muts_sum

def rank_variants(mutations_list, mutations_not_covered, compare_with, lin_mut_count, lin_del_list):
    '''
    generated rank variants table for one sample by comparing the sample's mutations with the mutationTable(Ha-Bodek).
    :param mutations_list: list of the sample's mutation
    :param mutations_not_covered: list of unknown mutations = Position is N
    :param compare_with:
    :param lin_mut_count: dict {lineage : mutations count}. mutations count(int) is the sum of unique mutations of a variant
    :param lin_del_list: dict {lineage : deletions lists}. deletions list(list of lists) each list is one deletion event
    :return: suspected_lineage: the variant that had the best match.
    highest covered (noN) mutations, highest noN mutations percentage
    :return: ranked_variant: a DataFrame of all variant and their match to the sample
    '''

    ranked_variants = pd.DataFrame(columns=["lineage", "mutations_fraction", "mutation_precentage",
                                            "noN_mutations_fraction","noN_mutation_precentage", "sort_by"])

    #create list of the mutations in each variant 
    for lin, lin_muts in compare_with.items():

        shared_mut = set( [mut for mut in mutations_list if mut.split("(")[0] in lin_muts] )
        noN_mutations = [x for x in lin_muts if (x not in mutations_not_covered) or (x in shared_mut)]
        if not shared_mut:
            continue
        shared_mut_count = sample_mut_counter(shared_mut, lin_del_list[lin])
        noN_mut_count = sample_mut_counter(noN_mutations, lin_del_list[lin])

        line = {        
            "lineage" : lin,
            "mutations_fraction": "("+str(shared_mut_count) + "/" + str(lin_mut_count[lin])+")",
            "mutation_precentage": round(shared_mut_count / lin_mut_count[lin] * 100, 2),
            "noN_mutations_fraction": "("+str(shared_mut_count) + "/" + str(noN_mut_count)+")",
            "noN_mutation_precentage": round(shared_mut_count / noN_mut_count * 100, 2) if noN_mut_count else 0,
            "sort_by": shared_mut_count
        }

        line['suspect'] = format_line(line)
        line_df = pd.DataFrame(data=line, index = [0])
        ranked_variants = pd.concat([ranked_variants, line_df], axis = 0) 
    
    if len(ranked_variants) == 0:
        return "", ranked_variants
    
    ranked_variants = ranked_variants.sort_values(by=['sort_by'], ascending=False,ignore_index=True)
    # determine suspect
    suspected_candidates = ranked_variants[ranked_variants['sort_by'] >= ranked_variants.iloc[0]['sort_by'] - 1]
    suspected_lineage = suspected_candidates.sort_values(by=['noN_mutation_precentage'], ascending=False,ignore_index=True).iloc[0]['lineage'] 

    return suspected_lineage, ranked_variants


def swap_finder(suspected_variant, suspected_recombinant, sample_mut):
    '''
    recombinant part function. get 2 variants names and the samples' mutations list (nucleotides)
    and find the recombination details.
    only consider mutations that are not in both variants.
    return values will be used for the recombinants files.
    :param suspected_variant: suspected variant name
    :param suspected_recombinant: suspected recombinant name
    :param sample_mut: a list of the samples mutations.
    :return var_muts: samples' mutations that are unique for the suspected variant.
    :return rec_muts: samples' mutations that are unique for the suspected recombinant.
    :return swap: number of recombination swaps between the two
    :return brkpnt: list of swaps Positions
    :return is_suspect: (bool) is suspected recombinant has less than 3 swaps and each swaps contains at list 2 mutations
     => is suspect = X
    '''

    is_suspect = ""  # if suspected recombinant is suspected after this function -> suspect = X
    
    # get unique mutations of suspected variant and its suspected recombinant (in relative to each other) - bodek.
    uniq_var_mut = [value for value in mutations_by_lineage_nt[suspected_variant]
                    if value not in mutations_by_lineage_nt[suspected_recombinant]]
    uniq_rec_mut = [value for value in mutations_by_lineage_nt[suspected_recombinant]
                    if value not in mutations_by_lineage_nt[suspected_variant]]

    var_muts = [value for value in uniq_var_mut if value in sample_mut]
    rec_muts = [value for value in uniq_rec_mut if value in sample_mut]
    
    merged_sorted_muts = mutTable[mutTable.Var_name.isin([suspected_variant,suspected_recombinant])]
    merged_sorted_muts = merged_sorted_muts[merged_sorted_muts.Nuc_sub.isin([*var_muts, *rec_muts])].sort_values(by=['Position']).reset_index()
    
    swap = 0  # count the swaps between variants
    muts_found = []  # mutation that were found already. to avoid repeats of following deletions.
    # following deletion have the same "AA_sub" in the BODEK
    mut_counter = -1  # count the mutations of the current variantdf
    mut_counter_list = []
    brkpnt =[]  # list of swaps Position
    for index, row in merged_sorted_muts.iterrows():
        if index == 0:
            var = row['Var_name']
            
        if row['AA_sub'] not in muts_found :
            mut_counter += 1
            muts_found.append(row['AA_sub'])
        temp = var
        var = row['Var_name']
        if var != temp:
            swap += 1
            mut_counter_list.append(mut_counter)
            mut_counter = 0
            brkpnt.append(str(merged_sorted_muts.iloc[index-1]["Position"])
                          + "-" + str(merged_sorted_muts.iloc[index]["Position"]))  # Position of variant swap
        # is_suspect test - suspected recombinant must have 1 or 2 swaps.
        # the variant and the recombinant must have at least 2 mutations in a row.
    mut_counter_list.append(mut_counter+1)  # append the last swap
    if swap in [1,2] and 1 not in mut_counter_list:
        is_suspect = "X" 
    return var_muts, rec_muts, swap, brkpnt, is_suspect

# get user input
alignment_file = argv[1]
output_file = argv[2]
pangolin_file = argv[3]
clades_path = argv[4]  # nextclade tsv
mutTable_path = argv[5]  # mutations table path
no_rec = argv[6]



output_path = Path(output_file).parent 
out_fname = datetime.now().strftime('%Y%m%d') + '_variants.csv'
output_file = output_path / out_fname



if len(argv) > 7:
    qc_report_path = argv[7]
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
    if "qc_status" in pangolinTable:
        pangolinTable = pangolinTable.rename(columns={"qc_status": "status"})
except FileNotFoundError:
    print("Pangolin File Not Found")
    pangolinTable = pd.DataFrame()

clades_df = pd.read_csv(clades_path, sep='\t')
red_flag = pd.read_csv(mutTable_path.replace("mutationsTable", "red_flag"))['mutation'].to_list()


# excel muttable:
mutTable = pd.read_csv(mutTable_path)

# prepare mutations table dataframe
mutTable["Reference"] = mutTable["Reference"].str.upper()  # make sure all upper case for safe comparisons
mutTable["Mutation"] = mutTable["Mutation"].str.upper()
mutTable = mutTable.dropna(thresh=3)

# prepare multiple alignment dictionary (key: sample name, val: SeqIO record)
alignment = SeqIO.to_dict(SeqIO.parse(alignment_file, 'fasta'))
alignment.pop('NC_045512.2', None)
alignment.pop('REF_NC_045512.2', None)

# prepare nextclade dataframe
clades_df = clades_df.rename(columns={'seqName': 'Sample'})
clades_df['Sample'] = clades_df['Sample'].apply(str)
clades_df = clades_df.fillna('')
# create dict of aasubstitutions and aadeletions. key=sample id, value=list of substitutions. to keep for final table.
aa_substitution_dict = {}
aa_deletions_dict = {}
insertions_dict = {}
for sample in clades_df['Sample']:  # change appearance from nextclade format to: Mutation(Gene)
    aasubs = clades_df[clades_df['Sample'] == sample].aaSubstitutions.values.tolist()
    aadels = clades_df[clades_df['Sample'] == sample].aaDeletions.values.tolist()
    insertions = clades_df[clades_df['Sample'] == sample].insertions.values.tolist()
    
    # ';' instead of ',' as sep. for Ari's tables -> important
    aa_substitution_dict[sample] = [f"{x.split(':')[1]}({x.split(':')[0]})" for x in aasubs[0].split(',')] \
        if (aasubs and aasubs != ['']) else ''
    aa_deletions_dict[sample] = [f"{x.split(':')[1]}({x.split(':')[0]})" for x in aadels[0].split(',')] \
        if (aadels and aadels != ['']) else ''
    insertions_dict[sample] = insertions[0].split(',')

mutations_by_lineage_nt = mutTable.groupby('Var_name')['Nuc_sub'].apply(list).to_dict()

# contains the variant mutations sum, consecutive deletions will be considered as one event
lin_mut_count = {}
# contains a list of lists of the deletions event. consecutive deletions will be stored in the same list
lin_del_list = {}

for lin, muts in mutations_by_lineage_nt.items():
    lin_del_list[lin], lin_mut_count[lin] = var_mut_counter(muts)


uniq_muts = mutTable.drop_duplicates(subset='Nuc_sub', keep="first")

#######################################################################

# DataFrames for csv outputs
final_table = pd.DataFrame(columns=["Sample", "variant (catalog)", "suspect","variant (bodek)", "bodek info", 
                                    "dictionary","% coverage","pango_usher", "pango_nextclade","AA substitutions",
                  "AA deletions", "Insertions", "mutations not covered", 
                  'nt substitutions', "recombinant suspect", "red_flag"])
chosen_recombinants = pd.DataFrame(columns=["Sample","suspected variant","suspected recombinant", "suspected variant mutations",
                                            "suspected recombinant mutations","swaps","breakpoints"])
ranked_variants_df = pd.DataFrame(columns=["sample", "suspect"])
all_sus_rec_df = pd.DataFrame(columns=["sample", "suspect","lineage", "mutations_fraction", "mutation_precentage",
                                            "noN_mutations_fraction","noN_mutation_precentage"])


extra_df = pd.DataFrame(columns=["sample", "mutations"])

# iterate over all samples in multi-fasta and over all mutations in table, and check value of each mutation
with open("results/mutations.log", 'w') as log:
    for sample, record in alignment.items():
        samples_mutations_nt = []
        samples_not_covered_nt = []
        if sample in aa_substitution_dict.keys():
            for (idx, row) in uniq_muts.iterrows():
                if pd.isna(row.loc['Position']):
                    log.write(f"NaN: {row}")
                    continue
                pos = int(row.loc['Position']) - 1  # mutation Position
                alt = record.seq[pos]  # fasta value in Position
                ref = row.loc['Reference']  # reference in Position
                table_mut = row.loc['Mutation']  # mutation  according to table
                if alt == table_mut:  # mutation exists in sequence
                    samples_mutations_nt.append(str(row.loc['Nuc_sub']))
                elif alt == 'N':  # if Position not covered in sequence
                    samples_not_covered_nt.append(str(row.loc['Nuc_sub']))

            
            sus_variant_name, rank_variant = rank_variants(samples_mutations_nt, samples_not_covered_nt,
                                                           mutations_by_lineage_nt, lin_mut_count,
                                                           lin_del_list)
            if not sus_variant_name:
                QCfail = True
            if len(rank_variant) > 0:
                rank_variant['sample'] = sample
                ranked_variants_df = pd.concat([ranked_variants_df,rank_variant], axis = 0)
                suspect_info = format_line(rank_variant.loc[rank_variant['lineage'] == sus_variant_name].iloc[0])
            else:
                suspect_info = ""
            if not suspect_info and samples_not_covered_nt:
            # not specific suspect variant but some mutations exist \ not covered in sequencing - write as suspect
                suspect_info = 'suspect'
            
            # get coverage of sample from qc report.txt
            if qc_report_path:
                coverage = qc[qc['sample'] == sample]['coverageCNS_5%'].values[0].round(2)
            else:
                coverage = str(calculate_coverage(alignment[sample].seq))
        
            
            QCfail = True if not sus_variant_name else False
        
            # get pangolin info from table
    
            try:
                pangolin_clade = pangolinTable[pangolinTable['taxon'] == sample].lineage.values[0]
            except:
                pangolin_clade = '-'

            # get nextclade info from table
            nextclade_pango = clades_df[clades_df['Sample'] == sample].Nextclade_pango
            
            nt_substitutions = clades_df.loc[clades_df['Sample'].eq(sample), 'substitutions'].str.split(',')
            nt_substitutions_list = ';'.join(nt_substitutions.values[0]).split(";")
            
            non_variant_mut_nt = set(extra_mutations(nt_substitutions_list, mutations_by_lineage_nt[sus_variant_name])) if not QCfail else ""


            #red flags            
            red_flags_list = [x for x in nt_substitutions_list if x in red_flag]


            extra = pd.DataFrame()
            extra["mutations"] = np.array(non_variant_mut_nt)
            extra["sample"] = sample
            extra_df = pd.concat([extra_df,extra], axis = 0)


            # new recombinant part
            is_rec_suspect = ''
            if not no_rec:
                if not QCfail:
                    sus_recombinant,ranked_recombinant = rank_variants(non_variant_mut_nt, samples_not_covered_nt,
                                                                       mutations_by_lineage_nt, lin_mut_count,
                                                                       lin_del_list)
                    if len(ranked_recombinant) > 0:
                        ranked_recombinant['sample'] = sample
                        all_sus_rec_df = pd.concat([all_sus_rec_df, ranked_recombinant], ignore_index=True)
                        var_muts, rec_muts, swap, brkpnt, is_rec_suspect = swap_finder(sus_variant_name, sus_recombinant
                                                                                       , nt_substitutions_list)
                        
                        chosen_recombinants.loc[len(chosen_recombinants)] = {
                                                                                "Sample": sample,
                                                                                "suspected variant": sus_variant_name,
                                                                                "suspected recombinant": sus_recombinant,
                                                                                "suspected variant mutations": var_muts,
                                                                                "suspected recombinant mutations": rec_muts,
                                                                                "swaps": swap,
                                                                                "breakpoints": brkpnt
                                                                            }
                        
                        
            sus_dict = mutTable.loc[mutTable['Var_name'] == sus_variant_name, 'Aliases']
            sus_dict = sus_dict.iloc[0] if not sus_dict.empty else ""
            sus_variant_name = sus_variant_name.split(".m")[0] #  Neta wants to remove the mediatore name from the output but not from the analysis.
            line = pd.DataFrame(data = {
                "Sample": sample,
                "variant (catalog)": sus_variant_name if suspect_info and not QCfail else "QC fail",
                "suspect": None,
                "variant (bodek)": sus_variant_name if suspect_info and not suspect_info=='suspect' 
                and int(suspect_info.split('noN:')[1].split(".")[0].split("%")[0]) > 50 else 'QC fail', #if the suspected variant < 60% => QC fail
                "bodek info": suspect_info,
                "dictionary": sus_dict if not str(sus_dict) ==  "nan" else "",
                'nt substitutions' : ';'.join(nt_substitutions_list),
                "AA substitutions": ';'.join(aa_substitution_dict[sample]) if aa_substitution_dict and
                                                                              sample in aa_substitution_dict else 'NA',
                "AA deletions": ';'.join(aa_deletions_dict[sample] if aa_deletions_dict and sample
                                                                     in aa_deletions_dict else 'NA'),
                "Insertions": ';'.join(insertions_dict[sample]) if insertions_dict and sample in insertions_dict else 'NA',
                "mutations not covered": ";".join(set(samples_not_covered_nt)) if samples_not_covered_nt else '',
                "% coverage": coverage,
                "recombinant suspect": is_rec_suspect,
                "pango_usher": pangolin_clade,
                "pango_nextclade" : nextclade_pango.values[0] if not nextclade_pango.empty or not QCfail else "QC fail",
                "red_flag"  : ';'.join(red_flags_list) if len(red_flags_list) >= 2 else ""
                }, index=[0])
            final_table = pd.concat([final_table,line], axis = 0)

final_table = final_table.reset_index(drop=True)
        
# append low quelity samples
low_quel = {id: [] for id in aa_substitution_dict if id not in alignment}
for sample in low_quel:
    if qc_report_path:
        coverage = qc[qc['sample'] == sample]['coverageCNS_5%'].values[0].round(2)
    else:
        coverage = ""


    final_table.loc[len(final_table)] = {
    "Sample": sample,
    "% coverage": coverage,
    "Variant": "QC fail",
    "suspected variant": "QC fail"
}

# for Matrix(the company):
# cutting the accession_id when there is an 'EPI_ISL' start #$
for index, row in final_table.iterrows():
    if 'EPI_ISL' in row["Sample"]: #$
                row["Sample"] = ''.join([i for i in row["Sample"].split('/')[-1].split('|') if 'EPI_ISL' in i]) #$


# #add resistant mutations 
# resist = resist.run(final_table, alignment)
# final_table = final_table.merge(resist, how='left', on='Sample')

#concat nextclade's output for Matrix(the company)
final_table = final_table.merge(clades_df, on='Sample', how='left')


final_table.to_csv(output_file,index=False)
ranked_variants_df[ranked_variants_df.columns[0:2]].to_csv(output_path / 'ranked_variants.csv',index=False)
extra_df.to_csv("results/non_variant_mutations.csv", index=False)
# recombinants files
if not no_rec:
    chosen_recombinants.to_csv("results/suspected_recombinants.csv", index=False)
    all_sus_rec_df[all_sus_rec_df.columns[0:2]].to_csv("results/ranked_suspected_recombinants.csv",index=False)



end = time.time()
total_time = end - start
print("variants.py time:\n"+ str(total_time)+ 'sec')

