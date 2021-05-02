import pandas as pd
from sys import argv
"""
create csv file from json nextclade output file.
filter only wanted information: clade, aa-substitutions
"""

def getAASubs(substitutions_list):
    """
    extract amino acid substitutions from 'substitutions' column.
    :param substitutions_list: list of substitutions (as dictionaries)
    :return: list of amino acid substitutions of the row.
    """
    row_list = []
    for dic in substitutions_list:
        if dic['aaSubstitutions']:
            codon = dic['aaSubstitutions'][0]['codon']
            gene = dic['aaSubstitutions'][0]['gene']
            refAA = dic['aaSubstitutions'][0]['refAA']
            queryAA = dic['aaSubstitutions'][0]['queryAA']
            sub = f"{refAA}{codon}{queryAA}({gene})"
            row_list.append(sub)
    return row_list


if __name__ == '__main__':
    # user input:
    clades_json_path = argv[1]
    output_path = argv[2]

    # open json file:
    clades_df = pd.read_json(clades_json_path)

    # get amino acid substitutions
    clades_df['aaSubstitutions'] = clades_df.apply(lambda row: getAASubs(row.substitutions), axis=1)

    # prepare output file
    output_df = clades_df[['seqName', 'aaSubstitutions', 'clade']]
    output_df.to_csv(output_path)
