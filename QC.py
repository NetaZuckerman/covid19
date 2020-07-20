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
    :param path: path to base project directory
    """
    global trimm_path
    if not os.path.exists("fastq/trimmed/"):
        os.mkdir("fastq/trimmed/")

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

            subprocess.call(command)  # run trimmomatic command. errors will be kept at trim_log


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


def fastqc_reports(out_dir, working_dir='fastq/raw/'):
    # TODO: add first/second qc options, to produce reports after trimming as well
    # fastqc_script_path = '/home/dana/covid19/fastqc_all.sh'
    # os.chmod(fastqc_script_path, 755)
    # with open('fastqc.log', 'w') as log:
    #     subprocess.call(['bash', fastqc_script_path, out_dir], stderr=log)

    for fqfile in os.listdir(working_dir):
        if not fqfile.endswith(".fastq.gz"):
            continue  # step over files that are not fastq.gz format
        print(working_dir+fqfile)
        with open("error.log", 'w') as log:
            subprocess.call(['fastqc', working_dir+fqfile, "--outdir=%s" % out_dir], stderr=log, stdout=log)

    print('finished producing reports')
    return out_dir


def multiqc_report(fq_out_dir, working_dir):
    subprocess.call(['multiqc', fq_out_dir, '-o', working_dir])
    print('finished multiqc')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-t", "--trim", help="QC trimming. Provide csv file with trimmomatic flags as columns. ",
                       type=str, metavar="[CSV]", nargs=1)
    group.add_argument("--template", help="Produce trimmomatic auto-trimmig template csv file with default values",
                       dest='template_fqpath', const="fastq/raw/", nargs="?")
    group.add_argument("-r", "--reports",  help="Produce fastqc and multifastqc reports of all fastq.gz files in input"
                                                "directory. Reports output path is optional, default is QC/fastqc/.",
                       type=str, dest='reports_outdir', nargs="?", const='QC/fastqc/')
    parser.add_argument("-w", "--working_dir", help="Working directory of the program. Default is current directory",
                        nargs=1, type=str, dest='wd')

    args = parser.parse_args()

    if args.wd:
        wd = args.wd
    else:
        wd = os.getcwd()

    if args.trim:
        print('trim')
        print(args.trim[0])
        if args.wd:
            print(args.wd[0])
            trim(args.trim[0], args.wd[0])
        else:
            print("not args wd")
            trim(args.trim[0])
        print('finished trimming. results are found in fastq/trimmed directory')

    elif args.reports_outdir:  # reports
        print('reports')
        fastqc_reports(args.reports_outdir, wd)
        multiqc_report(args.reports_outdir, wd)  # at the moment: multiqc produces report in pwd

    elif args.template_fqpath:
        if args.wd:
            print('No need to specify working directory when creating template. Exiting.')
            exit(1)
        # TODO: add PE SE options to create the right template
        # TODO: check SE as well
        print('producing template file')
        print(args.template_fqpath)
        template_csv(args.template_fqpath)

    else:  # if user did not choose an option -> print help
        parser.print_usage()
