# covid19

`usage: QC.py [-h] [-t CSV, base_path [CSV, base_path ...] | --template
             [FQ_PATH] | -r]`

### Examples:
#### Template:
To produce template of csv file to later serve as trimming input, use --template.\
If no path provided, will automatically use default path of "fastq/raw/".If your fastq files are in a different 
location, please provide a path.  
`python3 QC.py --template`\
`python3 QC.py --template /some/path/to/fastq/location/`\
The output: _template.csv_ file in your working directory.

Template example:\
| ends | input_forward | input_reverse | phred | threads | ILLUMINACLIP: |   |   |   |   |   |   |\
|------|---------------|---------------|-------|---------|---------------|---|---|---|---|---|---|
|   |   |   |   |   |   |   |   |   |   |   |   |
|   |   |   |   |   |   |   |   |   |   |   |   |
For more info about Trimmomatic's arguments, see Trimmomatic documentation:
>Trimmomatic 0.39 manual: http://www.usadellab.org/cms/?page=trimmomatic
>http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf


#### Trimming:
To trim run the following command:\
`python3 QC.py -t path/to/file.csv ` or
`python3 QC.py --trim path/to/file.csv`\
To trim with the template csv as input, simply run `python3 QC.py -t template.csv`. Make sure you are in the base
directory of the project so the program could access fastq/raw, fastq/trim etc.
If you prefer running the program from another location, please add the path to the base location, like so:
`python3 QC.py -t path/to/file.csv path/to/base/dir`



Dana Bar-Ilan.

last modified: 09.07.2020