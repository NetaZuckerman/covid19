# covid19

`usage: QC.py [-h] [-t CSV, base_path [CSV, base_path ...] | --template
             [FQ_PATH] | -r]`

### Usage:
#### Template:
To produce template of csv file to later serve as trimming input, use --template.\
If no path provided, will automatically use default path of "fastq/raw/".If your fastq files are in a different 
location, please provide a path.  \
The csv produced will include some default arguments, feel free to change them as you wish. 
`python3 QC.py --template`\
`python3 QC.py --template /some/path/to/fastq/location/`\
The output: _template.csv_ file in your working directory.


#### Trimming:
To trim run the following command:\
`python3 QC.py -t path/to/file.csv ` or
`python3 QC.py --trim path/to/file.csv`\
To trim with the template csv as input, simply run `python3 QC.py -t template.csv`. Make sure you are in the base
directory of the project so the program could access fastq/raw, fastq/trim etc.
If you prefer running the program from another location, please add the path to the base location, like so:
`python3 QC.py -t path/to/file.csv path/to/base/dir`\


### The CSV file:
Template csv header:
>ends   |	input_forward   |	input_reverse   |	phred   |	threads |   ILLUMINACLIP:   |	SLIDINGWINDOW:  
>|   LEADING:    |	TRAILING:   |	CROP:   |	HEADCROP:   |	MINLEN: 

Each line will represent a fastq file.

**Required fields**: ends, input_forward. Make sure those fields are filled in the csv file \
trimmomatic requires output files as well, but those are given automatically in the script. Don't include them in the csv.
* ends: `PE` or `SE` only. If you choose `PE` (paired end) -> be sure to include both input_forward **and**  input_reverse! 

**Optional fields**: 
Do not remove unwanted fields, just leave them empty and they will be ignored.
* -phred`[33 | 64]` -  -phred33 or -phred64 specifies the base quality encoding. If no quality encoding is specified,
it will be determined automatically 
* threads<int>: indicates the number of threads to use, which improves performance on multi-core
computers. If not specified, it will be chosen automatically. 
* ILLUMINACLIP:`<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip
threshold>:<simple clip threshold>` - Cut adapter and other illumina-specific sequences from the read. 
* SLIDINGWINDOW:`<windowSize>:<requiredQuality>` -  Perform a sliding window trimming, cutting once the average quality 
within the window falls below a threshold. 
* LEADING:/TRAILING:`<quality>` - Cut bases off the start/end of a read respectively, if below a threshold quality. 
quality: Specifies the minimum quality required to keep a base.
* CROP:`<length>`  - Removes bases regardless of quality from the end of the read, so that the read has maximally
the specified length after this step has been performed. length: The number of bases to keep, from the start of the read.
* HEADCROP:`<length>` - Removes the specified number of bases, regardless of quality, from the beginning of the read.
length: The number of bases to remove from the start of the read
* MINLEN:`<length>` - Removes reads that fall below the specified minimal length.  If required, it should
normally be after all other processing steps. 

NOTE: When adding new arguments to csv, add the argument as it will appear in the Trimmomatic command,
 including symbols like '-', ':', etc. 
 
>For more info about Trimmomatic's arguments, see Trimmomatic documentations:
>http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
>http://www.usadellab.org/cms/?page=trimmomatic


Dana Bar-Ilan.

last modified: 15.07.2020