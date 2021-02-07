from Bio import SeqIO
import csv
import pandas as pd
from sys import argv

# TODO: Less hard coded!
# CHECK RESULTS VS MUTATIONS TABLE?
# TEST NGS-29 by writing down all UK variants and their mutations


input_file = argv[1]
output_file = argv[2]

variantNames = {
    "UK": "B.1.1.7 - UK",
    "SA": "B.1.351 - SA",
    "Manaus": "P.1 - Manaus",
    "Rio": "P.2- Rio de jeneiro",
    "Mink": "Mink cluster 5",
    "Aus": "Australian 501 variant",
    "California": "B.1.429- California"
}

# mutTable = pd.read_csv("novelMutTable.csv")
mutTable = pd.read_csv("/data/projects/Dana/scripts/covid19/novelMutTable.csv") ############### TODO!

mutTable["REF"] = mutTable["REF"].apply(lambda x: x.upper())
mutTable["Mut"] = mutTable["Mut"].apply(lambda x: x.upper())


uk_muts = mutTable[mutTable["Lineage"] == "B.1.1.7 - UK"]["AA"].tolist()
sa_muts = mutTable[mutTable["Lineage"] == "B.1.351 - SA"]["AA"].tolist()
manaus_muts = mutTable[mutTable["Lineage"] == "P.1 - Manaus"]["AA"].tolist()
rio_muts = mutTable[mutTable["Lineage"] == "P.2- Rio de jeneiro"]["AA"].tolist()
mink_muts = mutTable[mutTable["Lineage"] == "Mink cluster 5"]["AA"].tolist()
aus_muts = mutTable[mutTable["Lineage"] == "Australian 501 variant"]["AA"].tolist()
cali_muts = mutTable[mutTable["Lineage"] == "B.1.429- California"]["AA"].tolist()


def isVar(seq, var, muttable=mutTable, vardict=variantNames):
    """
    Determine if sample is variant(var), and return all mutations types.

    :param seq: sample aligned fasta sequence
    :param var: variant short name (UK, SA, etc..) (vardict key)
    :param muttable: mutations table as pd.DataFrame
    :param vardict: dictionary of variant {short name: full name}
    :return: True or False (is Variant or not), variant mutations list,
    S uncovered mutations list and mutations that are different from table
    """
    lineage = vardict[var]
    is_var = True
    var_mutations = []  # variant mutation
    s_notCovered = []  # S mutation not covered (N nucleotide)
    diff_mutations = []  # mutation not by mutTable. !=ref && != mutation
    for line in muttable[muttable["Lineage"].str.contains(lineage)].iterrows():
        pos = int(line[1][5])-1  # genome position
        alt = seq[pos]    # sample's nucleotide in position
        ref = line[1][6]
        mutation = line[1][7]
        mutname = line[1][1]
        if "insertion" in mutname:
            continue
        gene = line[1][2]
        if alt != mutation or alt == 'N':
            is_var = False
            if alt == 'N' and gene == 'S':
                s_notCovered.append(mutname)
            elif alt != ref and alt != 'N':  # AND != Mutation (previous conditions)
                diff_mutations.append(mutname + "(alt:" + alt + ")")
        else:  # variant mutation
            var_mutations.append(mutname)
    return is_var, var_mutations, s_notCovered, diff_mutations


fastadict = SeqIO.to_dict(SeqIO.parse(input_file, 'fasta')) # fastadict -> {id: seq object}
fastadict.pop('NC_045512.2', None)  # remove refseq from dictionary (if does not exist, will do nothing - no Error)
fastadict.pop('REF_NC_045512.2', None)


mutTable["REF"] = mutTable["REF"].apply(lambda x: x.upper())


finalTable = []
for id, seqrecord in fastadict.items():
    seq = seqrecord.seq

    ukVar, ukMuts, uk_s, uk_diff = isVar(seq, "UK")
    saVar, saMuts, sa_s, sa_diff = isVar(seq, "SA")
    manausVar, manausMuts, manaus_s, manaus_diff = isVar(seq, "Manaus")
    rioVar, rioMuts, rio_s, rio_diff = isVar(seq, "Rio")
    minkVar, minkMuts, mink_s, mink_diff = isVar(seq, "Mink")
    ausVar, ausMuts, aus_s, aus_diff = isVar(seq, "Aus")
    caliVar, caliMuts, cali_s, cali_diff = isVar(seq, "California")

    known = ""
    extraMuts = ukMuts + saMuts + manausMuts + rioMuts + minkMuts + ausMuts + caliMuts
    sMuts = uk_s + sa_s + manaus_s + rio_s + mink_s + aus_s + cali_s
    different_muts = uk_diff + sa_diff + manaus_diff + rio_diff + mink_diff + aus_diff + cali_diff
    if ukVar:
        known = variantNames["UK"]
        extraMuts = [x for x in extraMuts if x not in uk_muts]
        sMuts = [x for x in sMuts if x not in uk_muts]
    elif saVar:
        known = variantNames["SA"]
        extraMuts = [x for x in extraMuts if x not in sa_muts]
        sMuts = [x for x in sMuts if x not in sa_muts]
    elif manausVar:
        known = variantNames["Manaus"]
        extraMuts = [x for x in extraMuts if x not in manaus_muts]
        sMuts = [x for x in sMuts if x not in manaus_muts]
    elif rioVar:
        known = variantNames["Rio"]
        extraMuts = [x for x in extraMuts if x not in rio_muts]
        sMuts = [x for x in sMuts if x not in rio_muts]
    elif minkVar:
        known = variantNames["Mink"]
        extraMuts = [x for x in extraMuts if x not in mink_muts]
        sMuts = [x for x in sMuts if x not in mink_muts]
    elif ausVar:
        known = variantNames["Aus"]
        extraMuts = [x for x in extraMuts if x not in aus_muts]
        sMuts = [x for x in sMuts if x not in aus_muts]
    elif caliVar:
        known = variantNames["California"]
        extraMuts = [x for x in extraMuts if x not in cali_muts]
        sMuts = [x for x in sMuts if x not in cali_muts]



    # extraMuts = set(extraMuts)  # unique mutations
    extraMuts = set([x + '(' + set(mutTable[mutTable["AA"] == x]["gene"]).pop() + ')' for x in extraMuts])
    sMuts = set(sMuts)


    line = {
            "Sample": id,
            "Known Variant": known if known else "no variant",
            "Suspect": "suspect" if (extraMuts or sMuts or different_muts) else "-",
            "More Mutations": ','.join(extraMuts),
            "S Not Covered": ','.join(sMuts),
            "non-Table Mutations": ','.join(different_muts)
            }
    finalTable.append(line)


with open(output_file, 'w') as outflie:
    fieldsnames = ["Sample", "Known Variant", "Suspect", "More Mutations", "S Not Covered", "non-Table Mutations"]
    writer = csv.DictWriter(outflie, fieldsnames,  lineterminator='\n')
    writer.writeheader()
    for line in finalTable:
        writer.writerow(line)
