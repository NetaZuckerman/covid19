import argparse  # for pretty help menu and parsing arguments
import csv  # for writing and reading csv files
import subprocess  # to run shell commands
import os  # to run os commands
import fnmatch  # regex
trimm_path='/data/software/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar'


class InputError(Exception):
    """
    customised exception class
    """
    def __init__(self, message=None):
        self.msg = message

    def __str__(self):
        if self.msg:
            return 'InputError: {0}'.format(self.msg)
        else:
            return 'InputError: Some required fields are missing in trim csv.'


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


def trim(flags_file, path='fastq/raw/', prefix_out="fastq/trimmed/"):
    """
    construct the trimmomatic command and execute it
    :param flags_file: csv input file
    :param path: path to raw fastq files
    :param prefix_out: Path to destination directory
    """
    global trimm_path
    if not os.path.exists(prefix_out):
        os.mkdir(prefix_out)

    with open(flags_file, 'r') as csv_file:
        file = csv.DictReader(csv_file)  # reads csv as dict. header is keys. each line is a dict.
        #  iterate over all lines(dicts) and gather all parameters to one command.
        for line in file:
            validate_input(line)
            ends = line['ends']
            command = ['java', '-jar', trimm_path, ends]

            if ends == 'PE':  # paired ends
                # create input and output files
                in_f = path + line['input_forward']
                in_r = path + line['input_reverse']
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
                inp = path + line['input_forward']
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

            subprocess.call(command)  # run trimmomatic command. errors will be kept at trim_log


# produce a template csv file used in trimming
def template_csv(fq_path='fastq/raw/', out_loc='template.csv'):
    """
    create template file with default arguments for trimmomatic
    :param path: path for all fastq.gz files.
    """

    fq_R1 = fnmatch.filter(os.listdir(fq_path), "*R1*.fastq.gz")  # Assuming reads are marked as "R1" and "R2"!!

    with open(out_loc, 'w') as output:
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

    print('finished template file. find it here: %s' %out_loc)


def fastqc_reports(out_dir="QC/", in_dir='fastq/raw/'):
    out_dest = out_dir + "fastqc/"
    if not os.path.exists(out_dest):
        os.mkdir(out_dest)
    for fqfile in os.listdir(in_dir):
        if not fqfile.endswith(".fastq.gz"):
            continue  # step over files that are not fastq.gz format
        print(in_dir + fqfile)
        with open("error.log", 'w') as log:
            subprocess.call(['fastqc', in_dir + fqfile, "--outdir=%s" % out_dest], stderr=log, stdout=log)
    print('Finished fastqc reports. Find them here: %s' % out_dest)


def multiqc_report(out_dir="QC/"):
    in_dir = out_dir + 'fastqc/'
    print('Start multiqc report')
    subprocess.call(['multiqc', in_dir, '-o', out_dir])
    print('Finished multiqc. Find it here: %s' % out_dir)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-t", "--trim", help="QC trimming. Provide csv file with trimmomatic flags as columns ",
                       type=str, metavar="[CSV]", nargs=1)
    group.add_argument("--template", help="produce trimmomatic auto-trimmig template csv file with default values",
                       action="store_true")
    group.add_argument("-r", "--reports",  help="produce fastqc and multifastqc reports of all fastq.gz files in input "
                                                "directory \nIf output path not provided (-o), default is QC/fastqc/",
                       action="store_true")
    parser.add_argument("-i", help="working directory of the program. \noptional, default is current directory",
                        nargs=1, type=str, dest='wd')
    parser.add_argument("-o", help="output path - optional.", nargs=1, type=str, dest='output_path')

    args = parser.parse_args()

    if args.wd:
        wd = args.wd[0]
    else:
        wd = ""

    out_path = ""
    if args.output_path:
        out_path = args.output_path[0]
        if not out_path.endswith('/'):
            out_path += '/'

    if args.trim:
        print('trim')
        print(args.trim[0])
        if wd and out_path:
            trim(args.trim[0], path=wd, prefix_out=out_path)
            print('finished trimming. results are found in %s provided' % out_path)
        elif wd or out_path:
            if wd:
                trim(args.trim[0], path=wd)
                print('finished trimming. results are found in fastq/trimmed directory')
            else:
                trim(args.trim[0], prefix_out=out_path)
                print('finished trimming. results are found in %s provided' % out_path)
        else:
            trim(args.trim[0])
            print('finished trimming. results are found in fastq/trimmed directory')

    elif args.reports:  # reports
        print('reports')
        if out_path:
            fastqc_reports(out_path)
            multiqc_report(out_path)
        else:
            fastqc_reports()
            multiqc_report()

    elif args.template:
        if wd and out_path:
            template_csv(fq_path=wd, out_loc=out_path)
        elif wd or out_path:
            if wd:
                template_csv(fq_path=wd)
            else:
                template_csv(out_loc=out_path)
        else:  # use defaults only
            template_csv()
        # TODO: add PE SE options to create the right template
        # TODO: check SE as well

    else:  # if user did not choose an option -> print help
        parser.print_help()
