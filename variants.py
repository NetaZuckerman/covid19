from Bio import SeqIO
import csv
import pandas as pd
from sys import argv

# TODO: Less hard coded!
# TODO: add to pipeline!

input_file = argv[1]
output_file = argv[2]

muts = {
    "UK": "B.1.1.7 - UK",
    "SA": "B.1.351 - SA",
    "Manaus": "P.1 - Manaus",
    "Rio": "P.2- Rio de jeneiro",
    "Mink": "Mink cluster 5",
    "Aus": "Australian 501 variant"
}

mutTable = pd.read_csv("/data/projects/Dana/scripts/covid19/novelMutTable.csv")
mutTable = pd.read_csv("novelMutTable.csv")
mutTable["REF"] = mutTable["REF"].apply(lambda x: x.upper())

uk_muts = mutTable[mutTable["Lineage"] == "B.1.1.7 - UK"]["AA"].tolist()
sa_muts = mutTable[mutTable["Lineage"] == "B.1.351 - SA"]["AA"].tolist()
manaus_muts = mutTable[mutTable["Lineage"] == "P.1 - Manaus"]["AA"].tolist()
rio_muts = mutTable[mutTable["Lineage"] == "P.2- Rio de jeneiro"]["AA"].tolist()
mink_muts = mutTable[mutTable["Lineage"] == "Mink cluster 5"]["AA"].tolist()
aus_muts = mutTable[mutTable["Lineage"] == "Australian 501 variant"]["AA"].tolist()


def isVar(seq, mut, muttable = mutTable, mutdict = muts):
    lineage = mutdict[mut]
    is_var = True
    var_mutations = []
    s_mutations = []
    for line in muttable[muttable["Lineage"].str.contains(lineage)].iterrows():
        pos = int(line[1][5])-1
        alt = seq[pos]
        ref = line[1][6]
        mutname = line[1][1]
        gene = line[1][2]
        if alt == ref or alt == 'N':
            is_var = False
            if alt == 'N' and gene == 'S':
                s_mutations.append(mutname)
        else:  # mutation
            if gene == 'S':
                s_mutations.append(mutname)
            else:
                var_mutations.append(mutname)
    return is_var, var_mutations, s_mutations


fastadict = SeqIO.to_dict(SeqIO.parse(input_file, 'fasta')) # fastadict -> {id: seq object}
fastadict.pop('NC_045512.2', None)  # remove refseq from dictionary (if does not exist, will do nothing - no Error)
fastadict.pop('REF_NC_045512.2', None)

mutTable = pd.read_csv("novelMutTable.csv")
mutTable["REF"] = mutTable["REF"].apply(lambda x: x.upper())


finalTable = []
for id, seqrecord in fastadict.items():
    seq = seqrecord.seq

    ukVar, ukMuts, uk_s = isVar(seq, "UK")
    saVar, saMuts, sa_s = isVar(seq, "SA")
    manausVar, manausMuts, manaus_s = isVar(seq, "Manaus")
    rioVar, rioMuts, rio_s = isVar(seq, "Rio")
    minkVar, minkMuts, mink_s = isVar(seq, "Mink")
    ausVar, ausMuts, aus_s = isVar(seq, "Aus")

    known = ""
    extraMuts = ukMuts + saMuts + manausMuts + rioMuts + minkMuts + ausMuts
    sMuts = uk_s + sa_s + manaus_s + rio_s + mink_s + aus_s
    if ukVar:
        known = muts["UK"]
        extraMuts = [x for x in extraMuts if x not in uk_muts]
        sMuts = [x for x in sMuts if x not in uk_muts]
    elif saVar:
        known = muts["SA"]
        extraMuts = [x for x in extraMuts if x not in sa_muts]
        sMuts = [x for x in sMuts if x not in sa_muts]
    elif manausVar:
        known = muts["Manaus"]
        extraMuts = [x for x in extraMuts if x not in manaus_muts]
        sMuts = [x for x in sMuts if x not in manaus_muts]
    elif rioVar:
        known = muts["Rio"]
        extraMuts = [x for x in extraMuts if x not in rio_muts]
        sMuts = [x for x in sMuts if x not in rio_muts]
    elif minkVar:
        known = muts["Mink"]
        extraMuts = [x for x in extraMuts if x not in mink_muts]
        sMuts = [x for x in sMuts if x not in mink_muts]
    elif ausVar:
        known = muts["Aus"]
        extraMuts = [x for x in extraMuts if x not in aus_muts]
        sMuts = [x for x in sMuts if x not in aus_muts]

    extraMuts = set(extraMuts)  # unique mutations
    sMuts = set(sMuts)

    line = {
            "Sample": id,
            "Known Variant": known if known else "no variant",
            "Suspect": "suspect" if (extraMuts or sMuts) else "-",
            "More Mutations": ','.join(extraMuts),
            "S Mutations": ','.join(sMuts)
            }
    finalTable.append(line)


with open(output_file, 'w') as outflie:
    fieldsnames = ["Sample", "Known Variant", "Suspect", "More Mutations", "S Mutations"]
    writer = csv.DictWriter(outflie, fieldsnames,  lineterminator='\n')
    writer.writeheader()
    for line in finalTable:
        writer.writerow(line)
