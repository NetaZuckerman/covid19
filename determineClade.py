from sys import argv
from Bio import SeqIO
import json

json_path = "/data/projects/Michal/nCoV2019/data/nextstrain_Israel/results/clades.json"  # change later for pipeline
with open(json_path) as f:
    data = json.load(f)

fastadict = SeqIO.to_dict(SeqIO.parse(argv[1], 'fasta'))
fastadict.pop('NC_045512.2', None)  # lose ref-seq. If not there will do nothing
fastadict.pop('REF_NC_045512.2', None)
samples_dict = {id: data['nodes'][id]['clade_membership'] for id, seqrecord in fastadict.items()}

with open('results/clades.csv', 'w') as file:
    file.write('id,clade\n')
    for key, val in samples_dict.items():
        file.write(f'{key},{val}\n')
