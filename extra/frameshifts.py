import pandas as pd
from Bio import SeqIO
from sys import argv

multifasta_path = argv[1]

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

def translate(sequence, start, end, codon_table):
    """
    translate nucleotides sequence in given ragion to amino acid sequence according to codon talbe
    BUT: frameshift -> close gaps?
    :param sequence:
    :param start:
    :param end:
    :param codon_table:
    :return:
    """

# translate all samples by mutations table and write down all frameshifts
muttable = pd.read_csv("../novelMutTable.csv")  # TODO: ADD SERVER ADDRESS
muttable = muttable[muttable.type == 'non -Synonymous ']

multifasta = SeqIO.to_dict(SeqIO.parse(multifasta_path, 'fasta'))
multifasta.pop("NC_045512.2", None)
multifasta.pop("REF_NC_045512.2", None)

sample = []
for file, seqrecord in multifasta.items():
    seq = seqrecord.seq
    frameshifts = []
    # for pos in muttable["pos"]:
