from sys import argv
from Bio import SeqIO, Seq, AlignIO
import pandas as pd
from Bio.Seq import translate
from Bio.Alphabet import Gapped, generic_dna, IUPAC
# import xlsxwriter
import Bio

# user input:
aligned_fasta_path = argv[1]
outfile_path = argv[2]
regions_table_path = argv[3]


def highlight_row(row):
    colors_list = [""] * 7 + ["background-color: grey"] * 2
    mut = row["mut"]
    ref = row["REF"]
    for samp in row[9:]:
        if samp != ref:
            if samp == mut:
                color = "background-color: yellow"
            elif samp != 'X':
                color = "background-color: lightyellow"
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
# get start and end of region
# don't forget -1 in position


def translate(sequence, start, end, codon_table):
    """
    translate nucleotides sequence in given region to amino acid sequence according to codons in region start-end
    :param sequence: nucleotides sequence as str #
    :param start: position of first nucleotide
    :param end: position of last nucleotide
    :return: translated sequence (aa seq)
    """
    tranlsated = []
    for i in range(start, end, 3):
        codon = sequence[i:i+3]
        if codon in codon_table:
            aa = codon_table[codon]
        else:
            aa = 'X' # ignore frameshift
        tranlsated.append(aa)
    return tranlsated


# 1. load sequences and tables
regionsTable = pd.read_csv(regions_table_path)
multifasta = SeqIO.to_dict(SeqIO.parse(aligned_fasta_path, 'fasta'))
multifasta.pop('NC_045512.2', None)
multifasta.pop('REF_NC_045512.2', None)

mutTable = pd.read_csv("/data/projects/Dana/scripts/covid19/novelMutTable.csv")
# 2. keep only non-synonymous mutations
mutTable = mutTable[mutTable.type == 'non -Synonymous ']
finalTable = mutTable
# 3. iterate over mutations and create
for sample, record in multifasta.items():
    seq = record.seq
    sample_muts = []
    for mut in mutTable.iterrows():
        pos = mut[1][5]
        gene = mut[1][2]
        aa = mut[1][1]
        aa_number = int(aa[1:-1])  # strip letters of both sides
        aa_from = aa[0]
        aa_to = aa[-1]
        region_start = int(regionsTable[regionsTable.id == gene].start.values[0])
        region_end = int(regionsTable[regionsTable.id == gene].end.values[0])
        region_translated = translate(str(seq), region_start-1, region_end, codon_map)
        alt = region_translated[aa_number-1]
        sample_muts.append(alt)
    finalTable[sample] = sample_muts


varcol = finalTable.apply(lambda row: row[8:].unique(), axis=1)
finalTable.insert(6, 'var', varcol)
finalTable = finalTable.sort_values(by=["gene", "lineage"], ascending=[False, True])
finalTable = finalTable.rename(columns={'pos': 'nuc pos', 'nucleotide': 'nuc name', 'AA': 'name'})

sorted_cols = ['nuc pos', 'nuc name', 'type', 'gene', 'var', 'name', 'lineage', 'REF', 'mut']
finalTable = finalTable[sorted_cols + [col for col in finalTable.columns if col not in sorted_cols]]

# write to file
finalTable["REF"] = finalTable.apply(lambda row: row["name"][0], axis=1)
finalTable["mut"] = finalTable.apply(lambda row: row["name"][-1], axis=1)
finalTable.style.apply(lambda row: highlight_row(row), axis=1).to_excel(outfile_path, index=False)
