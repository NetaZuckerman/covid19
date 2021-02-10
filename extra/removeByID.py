from Bio import SeqIO
from sys import argv

multifasta = argv[1]
newfasta = argv[2]
id_toremove = ['3464']


with open(multifasta) as original, open(newfasta, 'w') as corrected:
    records = SeqIO.parse(original, "fasta")
    for record in records:
        if record.id not in id_toremove:
            SeqIO.write(record, corrected, "fasta")

