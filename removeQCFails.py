from sys import argv
import pandas as pd
from Bio import SeqIO

# keep samples that has >50% coverage in QC report
if len(argv) < 4:
    print(f"usage: {argv[0]} <multi fasta> <qc_report> <output fasta>")

# user input:
fasta_path = argv[1]
qc_report_path = argv[2]
fasta_output = argv[3]


# upload QC report table
report = pd.read_csv(qc_report_path, sep='\t')
fasta_file = open(fasta_path, 'r')
for record in SeqIO.parse(fasta_file, 'fasta'):
    id = record.id
    # add try and except?
    line = report.loc[report['sample'] == id]['coverageCNS_5%']  # check value




