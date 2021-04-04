from Bio import SeqIO

# iterate over
fasta = list(SeqIO.parse('REF_NC_045512.2.fasta', 'fasta'))[0]
nCount = fasta.seq.upper().count('N')
length = len(fasta.seq)
