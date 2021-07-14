import argparse
import pandas as pd


def initialize_pangolin(path):
    """
    Initialize pangolin table DataFrame
    :param path: path to pangolin csv table
    :return: Data Frame made from pangolin table
    """
    try:
        pangolin = pd.read_csv(path)
    except FileNotFoundError:
        print(f"Pangolin table not found at {path}.")
        pangolin = pd.DataFrame()
    else:
        pangolin = pd.DataFrame()

    return pangolin


def initialize(arguments):
    """
    Initialize essentials.
    Get data from user input and prepare dataframes for later use.
    :param args: User input from argparse (namedtuple)
    :return: variables
    """
    # pangolin
    try:
        pangolin_table = pd.read_csv(arguments.pangolin)
    except FileNotFoundError or AttributeError:
        print(f"Pangolin table not found. Ignoring")
        pangolin_table = pd.DataFrame()
    else:
        pangolin_table = pd.DataFrame()

    # nextstrain
    try:
        nextstrain = pd.read_csv(arguments.nextstrain, sep='\t')
    except:
        pass
    return pangolin_table

if __name__ == '__main__':
    # get user input #
    # argparse tutorial: https://docs.python.org/3/howto/argparse.html
    # will present a manual in command line,
    # input will be found in args as named tuple

    parser = argparse.ArgumentParser()
    parser.add_argument("alignment",  help="alignment fasta file")
    parser.add_argument("--outfile", help="output file path")
    parser.add_argument("--pangolin", help="pangolin csv file")
    parser.add_argument("--nextclade", help="nextclade tsv file")
    parser.add_argument("mutations-table", help="excel mutations table")
    parser.add_argument("--qc report", help="qc report txt file")

    args = parser.parse_args()



