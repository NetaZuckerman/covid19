from sys import argv
from Bio import SeqIO, Seq, AlignIO
import pandas as pd
from Bio.Seq import translate
from Bio.Alphabet import Gapped, generic_dna, IUPAC
import Bio

# user input:
aligned_fasta_path = argv[1]
outfile_path = argv[2]
regions_table_path = argv[3]

codon_map = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "STOP", "TAG": "STOP",
    "TGT": "C", "TGC": "C", "TGA": "STOP", "TGG": "W",
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
    translate nucleotides sequence to amino acid sequence according to codons in region start-end
    :param sequence: nucleotides sequence (Seq? or str?) # TODO decide
    :param start: position of first nucleotide
    :param end: position of last nucleotide
    :return: translated sequence (aa seq)
    """


# 1. load sequences and tables
regionsTable = pd.read_csv(regions_table_path)
multifasta = SeqIO.to_dict(SeqIO.parse(aligned_fasta_path, 'fasta'))

mutTable = pd.read_csv("novelMutTable.csv")  # TODO: change to server's path
# 2. keep only non-synonymous mutations
mutTable = mutTable[mutTable.type == 'non -Synonymous ']
finalTable = mutTable
# 3. iterate over mutations and create
for sample, record in multifasta.items():
    seq = record.seq
    sample_muts = []
    print(sample)
    try:
        for mut in mutTable.iterrows():
            pos = mut[1][5]
            gene = mut[1][2]
            aa = mut[1][1]
            aa_number = int(aa[1:-1])  # strip letters of both sides
            aa_from = aa[0]
            aa_to = aa[-1]
            region_start = int(regionsTable[regionsTable.id == gene].start.values[0])
            region_end = int(regionsTable[regionsTable.id == gene].end.values[0])
            new_seq = full_codon_gaps(seq, region_start, region_end)

            region_aa = new_seq[region_start-1:region_end].translate(gap='-')
    except Bio.Data.CodonTable.TranslationError:
        print('sample: {} - Unable to translate'.format(sample))
        # add to list?



# # TODO: TO FIX -> A-- not translated. maybe translate with my own function??
# # finalTable.REF = finalTable.REF.str.upper()
#
