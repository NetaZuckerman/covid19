from Bio import SeqIO
from sys import argv

# user input: 1. multi-fasta file  2. output multi-fasta file path
fasta_path = argv[1]
outfile = argv[2]

# check user input, if missing print message and exit.
if len(argv) < 3:
    print(f"missing input.\nusage: {argv[0]} <input fasta> <output fasta> <optional: n% threshold>")
    exit(1)

threshold = 50  # threshold default, if not given by user.
if len(argv) > 3:
    threshold = argv[3]  # argv[3] optional

fasta = open(fasta_path, 'r')

with open(outfile, 'w') as out:
    # iterate over all fasta records in file,
    # keep only high quality ones (determined by threshold of N%)
    for record in SeqIO.parse(fasta, 'fasta'):
        sequence = record.seq
        nCount = sequence.upper().count('N')
        length = len(sequence)
        if (nCount/length)*100 >= threshold:  # if N percentage is higher/equals threshold, add to output file.
            SeqIO.write(record, out, 'fasta')

fasta.close()

