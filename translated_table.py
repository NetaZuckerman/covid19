from sys import argv
from Bio import SeqIO, Seq, AlignIO
import pandas as pd

# user input:
aligned_fasta_path = argv[1]
outfile_path = argv[2]
regions_table_path = argv[3]  # tables of regions of the genome, to determine translation reading frame in translation.
excel_mutations_table_path = argv[4]  # TODO: pipeline - add as argument

def highlight_row(row):
    """
    Highlight the mutations cells in excel, row by row.
    :param row: row to return it's colors
    :return: colors list matching row indices.
    """
    colors_list = [""] * 7 + ["background-color: silver"] * 2  # color of the fixed part of the table
    # (the mutation table part, first 8 columns)
    mut = row["mut"]
    ref = row["REF"]

    for samp in row[9:]:  # now color all other cells (each column belongs to a different sample)
        if samp != ref:  # highlight cell if it has a mutation.
            if samp == mut:  # highlight only if the mutation is matching the mutation in table
                color = "background-color: yellow"
            else:  # sample == 'X' -> do not color.
                color = ''
        else:
            color = ''
        colors_list.append(color)
    return colors_list


codon_map = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}


def full_codon_gaps(sequence, start, end, gap='-'):
    """
    avoid partial gaps in codon and convert to whole gaps codon
    exmple: from A-- to ---
    :param sequence: Seq object (Biopython)
    :param start: start of reading frame
    :param end: end of reading frame
    :param gap: gap character. default: '-'
    :return: new sequence, with full codons
    """

    old_seq = str(sequence)
    new_seq = old_seq[:start]
    for i in range(start-1, end, 3):
        codon = old_seq[i: i+3]
        if '-' in codon:
            codon = '---'
        new_seq += codon
    new_seq += old_seq[end:]
    return Seq.Seq(new_seq)


def translate(sequence, start, end, codon_table):
    """
    translate nucleotides sequence in given region to amino acid sequence according to codons in region start->end
    :param sequence: nucleotides sequence as str #
    :param start: position of first nucleotide
    :param end: position of last nucleotide
    :param codon_table: dictionary of codons as keys and AA name as values, for translation.
    :return: translated sequence (aa seq)
    """
    tranlsated = []
    for i in range(start, end, 3):
        codon = sequence[i:i+3]
        if codon in codon_table:
            aa = codon_table[codon]  # get the codon's matching amino acid by codon table dictionary
        else:
            aa = 'X'  # ignore frameshifts
        tranlsated.append(aa)
    return tranlsated


# 1. load sequences and tables
regionsTable = pd.read_csv(regions_table_path)
multifasta = SeqIO.to_dict(SeqIO.parse(aligned_fasta_path, 'fasta'))
multifasta.pop('NC_045512.2', None)  # remove refseq sequence from alignment file if exists.
multifasta.pop('REF_NC_045512.2', None)

mutTable_excel = pd.read_excel(excel_mutations_table_path, sheet_name=None)

for name in mutTable_excel:
    mutTable_excel[name]['lineage'] = name  # add a lineage column to all variant's tables

mutTable = pd.concat(mutTable_excel.values(), ignore_index=True)
# select only part of the columns:
mutTable = mutTable[['Position', 'Reference', 'Mutation', 'protein',
                     'variant', 'Mutation type', 'lineage', 'annotation']]
# compress identical mutations into one line and concat lineage names in the lineage column:
# mutTable = mutTable.groupby(  # to create compressed table:
#     ['Position', 'Reference', 'Mutation', 'protein', 'variant', 'Mutation type', 'annotation'], as_index=False).agg(
#     {'lineage': ';'.join}
# )

# 2. keep only non-synonymous mutations
# comparing in lower case to avoid mistakes such as SNP_Stop != SNP_stop. to catch all cases.
mutTable = mutTable[(mutTable['Mutation type'].str.lower() == 'snp') | (mutTable['Mutation type'].str.lower() == 'snp_stop')]
finalTable = mutTable

# 3. iterate over mutations and create final table.
for sample, record in multifasta.items():
    # each sample found in fasta will create a column of final table next to mutations info.
    seq = record.seq  # 1 fasta sequence as string
    sample_muts = []  # will aggregate translated value of each mutation to this list that will be added as column.
    for mut in mutTable.iterrows():
        pos = mut[1][0]
        gene = mut[1][3]
        aa = mut[1][4]
        aa_number = int(aa[1:-1])  # strip letters of both sides (snp mutations example: Q57H) letter-number-letter
        aa_from = aa[0]  # original AA
        aa_to = aa[-1]  # mutated AA
        # get start and end regions sequence location from regions table.
        region_start = int(regionsTable[regionsTable.id == gene].start.values[0])
        region_end = int(regionsTable[regionsTable.id == gene].end.values[0])
        # translate region (with translate() function in the top of the page)
        region_translated = translate(str(seq), region_start-1, region_end, codon_map)
        alt = region_translated[aa_number-1]  # get the specific aa by its number (in its name)
        sample_muts.append(alt)  # add to mutations list of the specific mutation to create a column to the final table
    finalTable[sample] = sample_muts  # add the column of the sample.


varcol = finalTable.apply(lambda row: row[8:].unique(), axis=1)  # add a 'var' column to list unique values of row
finalTable.insert(6, 'var', varcol)  # insert var column
finalTable = finalTable.sort_values(by=["lineage", "protein"], ascending=[True, False])  # sort by lineage and then gene
finalTable = finalTable.rename(columns={'Position': 'nuc pos', 'Mutation type': 'type', 'protein': 'gene',
                                        'variant': 'name', 'Reference': 'REF', 'Mutation': 'mut'})  # rename columns (old: new)

sorted_cols = ['nuc pos', 'type', 'gene', 'var', 'name', 'lineage', 'annotation', 'REF', 'mut']  # re-order columns
finalTable = finalTable[sorted_cols + [col for col in finalTable.columns if col not in sorted_cols]]  # re-order columns

# write to file
# add highlights with designated function in the top of the page
finalTable.style.apply(lambda row: highlight_row(row), axis=1).to_excel(outfile_path, index=False)  # add highlights

