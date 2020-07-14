import argparse
import csv
import subprocess
import os
trimm_path='/data/software/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar'

class InputError(Exception):
    def __init__(self, message=None):
        self.msg = message

    def __str__(self):
        if self.msg:
            return 'InputError: {0}'.format(self.msg)
        else:
            return 'InputError: Some required fields are missing in trim csv.'


def check_trim_input(line_dict):
    if not(line_dict['ends']):
        raise InputError
    if line_dict['ends'] == 'PE':
        if not (line_dict['input_forward'] and line_dict['input_reverse']):
            raise InputError()
    elif line_dict['ends'] == 'SE':
        if not line_dict['input']:
            raise InputError()
    else:
        raise InputError("ends may only be SE or PE")
    return


def parse_flags_from_csv(flags_file):
    global trimm_path
    with open(flags_file, 'r') as csv_file:
        file = csv.DictReader(csv_file)
        for line in file:
            check_trim_input(line)
            ends = line['ends']
            command = ['java', '-jar', trimm_path, ends]

            if ends == 'PE':  # paired ends
                in_f = 'fastq/raw/' + line['input_forward']
                in_r = 'fastq/raw/' + line['input_reverse']
                out_f = 'fastq/trimmed/' + line['input_forward'].rstrip('.fastq.gz')  # remove only from the end
                out_f_paired = out_f + '_paired.fastq.gz'
                out_f_unpaired = out_f + '_unpaired.fastq.gz'
                out_r = 'fastq/trimmed/' + line['input_reverse'].rstrip('.fastq.gz')  # remove only from the end
                out_r_paired = out_r + '_paired.fastq.gz'
                out_r_unpaired = out_r + '_unpaired.fastq.gz'
                command += [in_f, in_r, out_f_paired, out_f_unpaired, out_r_paired, out_r_unpaired]
            else:  # single ends
                inp = 'fastq/raw/' + line['input']
                out = 'fastq/trimmed/' + line['input'].rstrip('.fastq.gz')
                out_paired = out + '_paired.fastq.gz'
                out_unpaired = out + '_unpaired.fastq.gz'
                command += [inp, out_paired, out_unpaired]

            if line['threads']:
                command += ['-threads', line['threads']]

            if line['phred']:
                command.append('-phred'+line['phred'])

            for key, val in line.items():
                if (key in ['phred', 'threads', 'input', 'input_forward', 'input_reverse', 'ends']) or not val:
                    continue
                command += [key+val]

            print(command)
            subprocess.call(command)  # run trimmomatic command.


def template_csv():
    with open('template.csv', 'w') as output:
        headers = ['ends', 'input_forward', 'input_reverse', 'phred', 'threads', 'ILLUMINACLIP:', 'SLIDINGWINDOW:',
                   'LEADING:', 'TRAILING:', 'CROP:', 'HEADCROP:', 'MINLEN:']
        csv_writer = csv.writer(output)
        csv_writer.writerow(headers)


def fastqc_reports():  # TODO: maybe allow user to specify input and output
    # TODO: add first/second qc options, to produce reports after trimming as well
    try:
        fastqc_script_path='/home/dana/covid19/auto_QC.sh'
        os.chmod(fastqc_script_path, 755)
        with open('fastqc_error.log', 'w') as log:
            subprocess.call(['bash', '/home/dana/covid19/auto_QC.sh'], stderr=log)
    except Exception:
        pass
        print('Problem executing fastqc')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-t", "--trim", help="QC trimming. Provide csv file with trimmomatic flags as columns", type=str)
    group.add_argument("--template", help="Produce trimmomatic auto-trimmig template csv file", action='store_true')
    group.add_argument("-r", "--reports",  help="Produce fastqc reports of all fastq.gz files in input directory",
                       action='store_true')

    args = parser.parse_args()

    if args.trim:
        print('trim')
        print(args.trim)
        parse_flags_from_csv(args.trim)
        print('finished trimming. results are found in fastq/trimmed directory')

    elif args.reports:
        print('reports')
        fastqc_reports()

    elif args.template:
        # TODO: add PE SE options to create the right template
        # TODO: check SE as well
        print('producing template file')
        template_csv()
