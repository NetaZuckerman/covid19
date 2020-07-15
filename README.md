# covid19

`usage: QC.py [-h] [-t CSV, base_path [CSV, base_path ...] | --template
             [FQ_PATH] | -r]`

### Examples:
#### Template:
To produce template of csv file to later serve as trimming input, use --template.
--template: if no path provided, will automatically use default path of "fastq/raw/".If your fastq files (that needs trimming) are in a different location please provide a path.  
`python3 QC.py --template`\
`python3 QC.py --template /some/path/to/fastq/directory/`\
The output: template.csv file in your working directory.
#### Trimming:
To trim run the following command:\
`python3 QC.py -t path/to/file.csv ` or
`python3 QC.py --trim path/to/file.csv`\
To trim with the template csv as input, simply run `python3 QC.py -t template.csv`. Make sure you are in the base
directory of the project so the program could access fastq/raw, fastq/trim etc.
If you prefer running the program from another location, please add the path to the base location, like so:
`python3 QC.py -t path/to/file.csv path/to/base/dir`


For more info:
>Trimmomatic 0.39 manual: http://www.usadellab.org/cms/?page=trimmomatic



Dana Bar-Ilan.

last modified: 09.07.2020