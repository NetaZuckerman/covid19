import argparse  # for pretty help menu and parsing arguments
import csv  # for writing and reading csv files
import subprocess  # to run shell commands
import os  # to run os commands
import re  # regex
import fnmatch # regex
trimm_path='/data/software/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar'


# customised exception class
class InputError(Exception):
    def __init__(self, message=None):
        self.msg = message

    def __str__(self):
        if self.msg:
            return 'InputError: {0}'.format(self.msg)
        else:
            return 'InputError: Some required fields are missing in trim csv.'


# validate csv input (receives line in csv)
def validate_input(line_dict):
    """
    validate csv file (ends vs input provided)
    :param line_dict:
    """
    if not(line_dict['ends']):
        raise InputError

    if not line_dict['input_forward']:
        raise InputError("No input file provided. Please fix csv file and try again.")
    if line_dict['ends'] == 'PE' and not line_dict['input_reverse']:
        raise InputError("No reverse input file, please provide one or choose SE instead of PE.")

    elif line_dict['ends'] not in ['PE', 'SE']:
        raise InputError("ends may only be SE or PE.")
    return


def trim(flags_file, path=''):
    """
    construct the trimmomatic command and execute it
    :param flags_file: csv input file
    """
    global trimm_path
    prefix_fq = path + "fastq/raw/"
    prefix_out = path + "fastq/trimmed/"
    with open(flags_file, 'r') as csv_file:
        file = csv.DictReader(csv_file)  # reads csv as dict. header is keys. each line is a dict.
        #  iterate over all lines(dicts) and gather all parameters to one command.
        for line in file:
            validate_input(line)
            ends = line['ends']
            command = ['java', '-jar', trimm_path, ends]

            if ends == 'PE':  # paired ends
                # create input and output files
                in_f = prefix_fq + line['input_forward']
                in_r = prefix_fq + line['input_reverse']
                out_f = prefix_out + line['input_forward'].rstrip('.fastq.gz')
                out_f_paired = out_f + '_paired.fastq.gz'
                out_f_unpaired = out_f + '_unpaired.fastq.gz'
                out_r = prefix_out + line['input_reverse'].rstrip('.fastq.gz')
                out_r_paired = out_r + '_paired.fastq.gz'
                out_r_unpaired = out_r + '_unpaired.fastq.gz'
                # add input and output files to final command
                command += [in_f, in_r, out_f_paired, out_f_unpaired, out_r_paired, out_r_unpaired]
            else:  # single ends
                # create input and output files
                inp = prefix_fq + line['input_forward']
                out = prefix_out + line['input_forward'].rstrip('.fastq.gz')
                out_paired = out + '_paired.fastq.gz'
                out_unpaired = out + '_unpaired.fastq.gz'
                # add input and output files to final command
                command += [inp, out_paired, out_unpaired]

            if line['threads']:  # if threads argument is filled by user
                command += ['-threads', line['threads']]

            if line['phred']:  # if phred argument is filled by user
                command.append('-phred'+line['phred'])

            for key, val in line.items():  # iterate over all arguments. skip checked or empty arguments.
                # gather all additional arguments into the command. key is the trimmomatic flag and val is the value
                # filled by user.
                if (key in ['phred', 'threads', 'input', 'input_forward', 'input_reverse', 'ends']) or not val:
                    continue
                command += [key+val]

            print(command)
            with open('trimmomatic.log') as trim_log:
                subprocess.call(command, stderr=trim_log)  # run trimmomatic command. errors will be kept at trim_log


# produce a template csv file used in trimming
def template_csv(fq_path):
    """
    create template file with default arguments for trimmomatic
    :param path: path for all fastq.gz files.
    """

    fq_R1 = fnmatch.filter(os.listdir(fq_path), "*R1*.fastq.gz")  # Assuming reads are marked as "R1" and "R2"!!

    with open('template.csv', 'w') as output:
        headers = ['ends', 'input_forward', 'input_reverse', 'phred', 'threads', 'ILLUMINACLIP:', 'SLIDINGWINDOW:',
                   'LEADING:', 'TRAILING:', 'CROP:', 'HEADCROP:', 'MINLEN:']
        writer = csv.DictWriter(output, fieldnames=headers)
        writer.writeheader()
        for file in fq_R1:
            writer.writerow({
                'ends': 'PE',
                'input_forward': file,
                'input_reverse': file.replace('R1', 'R2'),
                'phred': '33',
                'threads': '32',
                'TRAILING:': '28',
            })

    print('finished template file. check it out at template.csv in your working directory.')


def fastqc_reports():  # TODO: maybe allow user to specify input and output
    # TODO: add first/second qc options, to produce reports after trimming as well
    try:
        fastqc_script_path='/home/dana/covid19/fastqc_all.sh'
        os.chmod(fastqc_script_path, 755)
        with open('fastqc.log', 'w') as log:
            subprocess.call(['bash', fastqc_script_path], stderr=log)
    except Exception:
        pass
        print('Problem executing fastqc')
    print('finished producing reports')


def multiqc_report():
    try:
        subprocess.call(['multiqc', 'QC/fastqc', '-o', 'QC/'])
    except:
        print("Problem executing multiqc")
    print('finished multiqc')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-t", "--trim", help="QC trimming. Provide csv file with trimmomatic flags as columns. "
                                            "base_path: path to base directory, where fastq/ and QC/ dirs are found.",
                       type=str, metavar="CSV, base_path", nargs="+")
    group.add_argument("--template", help="Produce trimmomatic auto-trimmig template csv file with default values",
                       dest='fq_path', const="fastq/raw/", nargs="?")
    group.add_argument("-r", "--reports",  help="Produce fastqc and multifastqc reports of all fastq.gz files in input"
                                                "directory", action='store_true')

    args = parser.parse_args()

    if args.trim:
        print('trim')
        print(args.trim)
        if len(args.trim) == 1:
            trim(args.trim[0])
        elif len(args.trim) == 2:
            trim(args.trim[0], args.trim[1])
        else:
            print("too many arguments for trim option. see -h for help.")
            exit(1)
        print('finished trimming. results are found in fastq/trimmed directory')

    elif args.reports:
        print('reports')
        fastqc_reports()
        multiqc_report()

    elif args.fq_path:
        # TODO: add PE SE options to create the right template
        # TODO: check SE as well
        print('producing template file')
        print(args.fq_path)
        template_csv(args.fq_path)

    else: # if user did not choose an option -> print help
        parser.print_usage()