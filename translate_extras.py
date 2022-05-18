import pandas as pd
from Bio import SeqIO
from math import floor
from sys import argv

# read files
region_file = argv[1]
results_file = argv[2]
ref_file = argv[3]
alignment_file = argv[4]

regions = pd.read_csv (region_file)
df = pd.read_csv(results_file)
new_df = pd.DataFrame(columns=["sample", "mutations", "aa_mut"])
# read sequences
reference = list(SeqIO.read(ref_file, 'fasta').seq._data)
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

def get_gene(position):
    for index, row in regions.iterrows():
        gene = None
        if position in range(row["START"],row["END"]):
            gene = row["GENE"]
            break
    return gene, row["START"], row["END"]

def get_codon_pos(position, start):
    pos_on_gene = position - start + 1
    aa_pos_on_gene = floor(pos_on_gene/3) + 1
    mod = pos_on_gene % 3
    if mod == 0:
        aa_pos_on_gene =- 1
        codon = (position - 2, position - 1, position)
    if mod == 1:
        codon = (position, position + 1, position + 2)
    if mod == 2:
        codon = (position - 1, position, position + 1)

    return codon, aa_pos_on_gene

def get_codon(seq, codon_pos):
    return seq[codon_pos[0]-1] + seq[codon_pos[1]-1] + seq[codon_pos[2]-1]



for index, row in df.iterrows():
    nt_position = int(row["mutations"][1:-1])
    gene, start, end = get_gene(nt_position)

    if gene == None:
        row["aa_mut"] = "UTR"
    else:
        codon_pos, aa_pos = get_codon_pos(nt_position, start)
        ref_codon = get_codon(reference,codon_pos)
        seq_codon = get_codon(alignment[row["sample"]], codon_pos)
        ref_aa = translate_table[ref_codon]
        seq_aa = translate_table[seq_codon]

        row["aa_mut"] = gene + ":" + ref_aa + str(aa_pos) + seq_aa if not ref_aa == seq_aa else gene + ": silent"
    new_df = new_df.append(row)

new_df.to_csv(results_file, index=False)
