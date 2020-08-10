
from collections import defaultdict
from Bio import SeqIO
import argparse


def mask_fasta(fasta_record, positions_list):
    """
    Mask fasta file to N in specified positions.
    :param fasta_record: fasta SeqIO record
    :param positions_list: list of positions to replace to 'N'
    :return: the masked fasta
    """
    mutable_seq = fasta_record.seq.tomutable()
    for position in positions_list:
        mutable_seq[position] = 'N'
    # new record, the same with new sequence
    return SeqIO.SeqRecord(mutable_seq, id=fasta_record.id, description=fasta_record.description,
                           name=fasta_record.name, dbxrefs=fasta_record.dbxrefs)


def parse_depth(depth_file, n=1):
    """
    Parse samtools depth file - each position -1 to match str index.
    :param depth_file: Position's dept, as samtools depth command's output.
    :param n: Masking threshold. Every position with depth<n will be changed to 'N'. default:1
    :return: default dict of {<chr>:[pos1,pos2,...] record name with list of positions to mask.
    """
    n_pos = defaultdict(list)  # chr1:[pos1,...], chr2:[...]
    with open(depth_file) as df:
        for line in df:  # line[0]:chr  line[1]:pos line[2]:depth
            line = line.split('\t')
            if int(line[2]) < n:
                n_pos[line[0]].append(int(line[1])-1) # one to the left - depth always starts with 1 instead of 0.
    return n_pos


def mask_multiple_records(fasta_file_name, depth_file, n=1):
    """

    :param fasta_file_name: fasta file to mask. contains 1 or more records.
    :param depth_file:depth of positions. as samtools depth output.
    :param n: Masking threshold. Every position with depth<n will be changed to 'N'. default:1
    :return: list of fasta records
    """
    depth_dict = parse_depth(depth_file, n)
    all_records = [mask_fasta(record, depth_dict[record.id]) for record in SeqIO.parse(fasta_file_name, "fasta")]
    return all_records


if __name__ == '__main__':
    # usage: python mask_fasta.py <fasta_file> <depth_file> <N masking-threshold>
    parser = argparse.ArgumentParser(description="Mask fasta file with Ns according to depth of each position. "
                                                 "Supports multiple fasta in one file as well.")
    parser.add_argument('fasta', type=str, help="fasta file to mask.")
    parser.add_argument('depth', type=str, help="depth file as samtools depth command's output.")
    parser.add_argument('-n', type=int, help="depth threshold. every position with depth lower than N will be "
                                                      "masked as 'N'. default: 1")
    args = parser.parse_args()
    N = 1
    if args.n:
        N = args.n
    all_fasta_records = mask_multiple_records(args.fasta, args.depth, N)
    SeqIO.write(all_fasta_records, args.fasta, "fasta")






