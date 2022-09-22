
import pandas as pd
from Bio import SeqIO
from math import floor
from sys import argv


# read files
region_file = argv[1]
ref_file = argv[2]
alignment_file = argv[3]
nextclade_path = argv[4]
output = argv[5]
regions = pd.read_csv (region_file)
new_df = pd.DataFrame(columns=["accession_id","position","referance","mutation","nuc_sub","nuc_var", \
                                   "genetic_region","genetic_region_orf","variant","mutation_type","variant_orf", \
                                       "mutation_type_orf","var_name","var_name_orf"])
# read sequences
reference = list(SeqIO.read(ref_file, 'fasta').seq)
alignment = SeqIO.to_dict(SeqIO.parse(alignment_file, 'fasta'))

for sample, record in alignment.items():
        alignment[sample] = list(str(record.seq).upper())

translate_table = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
}

def open_deletions(del_list):
    muts = []
    for deletion in del_list:
       pos = deletion.split("-")
       pos_list = list(range(int(pos[0]),int(pos[1])+1)) if len(pos) > 1 else pos
       pos[0] = int(pos[0])
       for pos in pos_list:
           ref_nuc = reference[pos]
           muts.append(ref_nuc + str(pos) + '-')
    return muts
    
def sort_by_pos(df):
    for index, row in df.iterrows():
        df.at[index, "pos"] = int(row["nuc_sub"][1:-1])
    df = df.sort_values(by=['pos'])
    
    return df[["accession_id", "nuc_sub"]]

def get_mut_df(nextclade_path):
    df = pd.DataFrame()    
        
    clades = pd.read_csv(nextclade_path,sep='\t')    
    
    for index, row in clades.iterrows():
    
        temp = pd.DataFrame()
        temp["nuc_sub"] = row["substitutions"].split(',')
        deletions = pd.DataFrame({"nuc_sub":open_deletions(row["deletions"].split(','))})
        temp = temp.append(deletions,ignore_index=True)
        
        temp["accession_id"] = row["seqName"]
        sorted_temp = sort_by_pos(temp)
        df = df.append(sorted_temp)
    return df
    

def get_gene(position):
    gene_list = []
    for index, row in regions.iterrows():
        if position in range(row["start"],row["end"]):
            gene_list.append([row["region"], row["start"], row["end"]])
    return gene_list

def get_codon_pos(position, start):
    pos_on_gene = position - start + 1
    aa_pos_on_gene = floor(pos_on_gene/3) + 1
    mod = pos_on_gene % 3
    if mod == 0:
        aa_pos_on_gene -= 1
        codon = (position - 2, position - 1, position)
    if mod == 1:
        codon = (position, position + 1, position + 2)
    if mod == 2:
        codon = (position - 1, position, position + 1)

    return codon, aa_pos_on_gene

def get_codon(seq, codon_pos):
    return seq[codon_pos[0]-1] + seq[codon_pos[1]-1] + seq[codon_pos[2]-1]

def get_aa(gene, start, nt_position):
    if "UTR" in gene:
        return ""
    codon_pos, aa_pos = get_codon_pos(nt_position, start)
    ref_codon = get_codon(reference, codon_pos)
    seq_codon = get_codon(alignment[str(row["accession_id"])], codon_pos)
    ref_aa = translate_table[ref_codon]
    if seq_codon == "---":
        seq_aa = '-'
    elif "-"  in seq_codon or "N"  in seq_codon:
        seq_aa = 'X'
    else:
        seq_aa = translate_table[seq_codon]
    return ref_aa + str(aa_pos) + seq_aa 

def mut_type(var): # var is a mutation format ref-position-alt example - A1232I
    if not var:
        return ""
    if var[-1] == var[0]:
        return "SNP_silent"
    if var[-1] == "-":
        return "Deletion"
    if not var[-1] == var[0]:
        return "SNP"
    
df = get_mut_df(nextclade_path)
for index, row in df.iterrows():
    nt_position = int(row["nuc_sub"][1:-1])
    gene_list = get_gene(nt_position)
    gene1 = gene2 = var1 = var2 = ""


    if len(gene_list) > 0:
        gene = gene_list[0][0]
        start = gene_list[0][1]
        var1 = get_aa(gene, start, nt_position)
        gene1 = gene

    if len(gene_list) == 2 :  # NSPs genes.
        gene = gene_list[1][0]
        start = gene_list[1][1]
        var2 = get_aa(gene, start, nt_position)
        gene2= gene

    row["nuc_var"] = row["nuc_sub"]
    row["position"] = nt_position
    row["referance"] = row["nuc_sub"][0]
    row["mutation"] = row["nuc_sub"][-1]

    if "ORF" in gene1: # case - (ORF7a, ORF7b) or  (ORF1 , NSP1) or (ORF, NUCAP)
        row["genetic_region"] = gene2
        row["genetic_region_orf"] = gene1
        row["variant"] = var2
        row["variant_orf"] = var1
        row["mutation_type"] = mut_type(var2)
        row["mutation_type_orf"] = mut_type(var1)
        row["var_name"] = row["genetic_region"] + ":" + row["variant"] if row["variant"] else ""
        row["var_name_orf"] = row["genetic_region_orf"] + ":" + row["variant_orf"] if row["variant_orf"] else ""
        
    if "ORF" in gene2 or ( len(gene_list) < 2 and "ORF" not in gene1): # case (NSP1, ORF1) or (UTR, _)
        row["genetic_region"] = gene1
        row["genetic_region_orf"] = gene2
        row["variant"] = var1
        row["variant_orf"] = var2
        row["mutation_type"] = mut_type(var2)
        row["mutation_type_orf"] = mut_type(var1)
        row["var_name"] = row["genetic_region"] + ":" + row["variant"] if row["variant"] else ""
        row["var_name_orf"] = row["genetic_region_orf"] + ":" + row["variant_orf"] if row["variant_orf"] else ""
        
        
    new_df = new_df.append(row)

# cutting the accession_id when there is an 'EPI_ISL' start #$
for index, row in new_df.iterrows():    
    if 'EPI_ISL' in row["accession_id"]: #$
        row["accession_id"] = ''.join([i for i in row["accession_id"].split('/')[-1].split('|') if 'EPI_ISL' in i]) #$
new_df.to_csv(output, index=False)
