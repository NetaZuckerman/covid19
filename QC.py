import argparse
import csv
import subprocess


def parse_flags_from_csv(flags_file):
    trimm_path='/data/software/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar'
    with csv.DictReader(open(flags_file)) as file:
        for line in file:
            sample = line['sample']
            ends = line['ends']
            phred = "-" + line['phred']
            if ends not in ('PE' or 'SE'):
                print("%s: ends must be 'PE' or 'SE'. Fix and run again.")
            # subprocess.call(['java', '-jar', trimm_path, ends, phred]) # TODO: decide how to run trimming command
            # MUST FINISH BY TUESDAY!




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



