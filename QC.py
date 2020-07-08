import argparse

parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group()
parser.add_argument("input_dir", help="Provide input directory that contains all fastq.gz samples")
parser.add_argument("output_dir", help="Provide output directory")
group.add_argument("-t", "--trim", help="QC trimming. Provide csv file with trimmomatic flags as columns")
group.add_argument("-r", "--reports",  help="Produce fastqc reports of all fastq.gz files in input directory",
                    action="store_true")

args = parser.parse_args()
in_dir = args.input_dir
out_dir = args.output_dir

