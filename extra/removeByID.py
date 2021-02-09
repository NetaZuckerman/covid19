from Bio import SeqIO
from sys import argv

multifasta = argv[1]
newfasta = argv[2]
id_toremove = ['s6831','s6933', 's6912']


with open(multifasta) as original, open(newfasta) as corrected:
    records = SeqIO.parse(original, "fasta")
    for record in records:
        if record.id not in id_toremove:
            SeqIO.write(record, corrected, "fasta")

