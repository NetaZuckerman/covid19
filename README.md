# covid19
add name to 
## Suggested Workflow Using QC.py and pipeline.sh
1. Create project's directories - Optional. \
`bash pipeline.sh -d` \
The directories will be created in the current working directory, and it is recommended to run all following commands 
from that working directory.
* Fill fastq/raw/ with all raw fastq.gz files. 
* Fill refs/ with your desired reference sequence.
To see the full directories hierarchy scroll to the bottom of the page.

2. Produce reports for each fastq.gz file in fastq/raw (or any other location of your choice) \
`python3 QC.py --reports -i fastq/raw/`

3. Create template csv file for convenient trimming. See details at Template. \
`python3 QC.py --template -i fastq/raw/ ` \
Take a look at _template.csv_ produced as output and fix some trimming parameters according to fastqc/multiqc reports.
Scroll down to read more about the csv file.  

4. Trim fastq files according to the csv file created in (2). \
`python3 QC.py --trim template.csv -i /input/path/to/fastq.gz.files -o /trimmed/path` 

5. Produce second reports to make sure you're happy with the trimming results. \
`python3 QC.py --reports -i /trimmed/path -o /output/fastqc/2/path` \
Be sure to provide the path to the trimmed fastq.gz this time. Repeat steps 2-5 if more trimming is needed.

6. Run pipeline to produce consensus sequences for your samples, and get additional info such as number of reads, depth, 
coverage, etc. in results/report.txt. Run the pipeline on the trimmed fastq.gz files (the pipeline uses 
the paired files, ignoring singletons). 
To run the pipeline with raw fastq files omit the _--trimmed_fq_ flag. \
`bash pipeline.sh --trimmed_fq -i fastq/trimmed/ -r /refs/refseq.fa` 

# QC.py
### Usage:
`python QC.py [-h] [-t [CSV] | --template | -r] [-i WD] [-o OUTPUT_PATH]`
#### QC Reports:
##### -r | --reports
Produces fastqc and multiqc reports for all fastq.gz files. \
`python3 QC.py -r` \
You may also provide input (-i) and/or output (-o) paths: \
`python3 QC.py -r -i path/to/fastq/files -o path/for/output/reports`
default -i location: fastq/raw/ \
default -o location: QC/ 

#### Template:
##### --template 
Produces template of csv file with some default trimmomatic arguments.\
`python3 QC.py --template`\
You may also provide input (-i) and/or output (-o) paths: \
`python3 QC.py --template -i some/path/to/fastq/location/ -o some/path/for/output` \
default -i location: fastq/raw/ \
default -o location: current working directory
The output file is called template.csv, it will be located in current directory, or output path (-o) if provided.

#### Trimming:
##### -t|--trim <path/to/csv> <optional-path/to/base/dir>
To trim run the following command:\
`python3 QC.py -t path/to/file.csv `\ or
`python3 QC.py --trim path/to/file.csv`\
You may also provide input (-i) and/or output (-o) paths: \
`python3 QC.py -t path/to/file.csv -i path/to/raw/fastq/files -o path/to/trimmed/fatsq/output`
To trim with the template csv as input, simply run `python3 QC.py -t template.csv`.\
default -i location: fastq/raw/ \
default -o location: fastq/trimmed/

### The CSV file:
Template csv header:
>ends, input_forward, input_reverse, phred, threads, ILLUMINACLIP:, SLIDINGWINDOW:, LEADING:, TRAILING:, CROP:,
>HEADCROP:, MINLEN: 

Each line will represent a fastq file.

**Required fields**: ends, input_forward. Make sure those fields are filled in the csv file \
trimmomatic requires output files as well, but those are given automatically in the script. Don't include them in the csv.
* ends`[PE|SE]` - PE:paired end, SE: single end. If you choose `PE`, be sure to include both input_forward **and**  input_reverse! 
* input_forward `<sample_R1.fastq.gz>`: sample name
* input_reverse `<sample_R2.fastq.gz>`: sample name

**Optional fields**: 
Do not remove unwanted fields, just leave them empty and it will be ignored.
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
 including symbols like -, :, etc. 
 
>For more info about Trimmomatic's arguments, see Trimmomatic documentations:
>http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf,
>http://www.usadellab.org/cms/?page=trimmomatic
---------------
# pipeline.sh
### Usage:
 bash pipeline.sh [options]

#### Help:
##### -h | --help

#### Create project's directories:
##### -d|--create_dirs 
Create all the project directories in the current directory and exit. \
The directories: fastq/raw, fastq/trimmed, QC/fastqc, BAM, CNS, alignment, Trees, results
`bash pipeline.sh -d` OR `bash pipeline.sh --crate_dirs`

#### Run pipeline on samples to get consensus sequences and additional data:
#### Provide input path:
##### -i /input/path (required)
`/input/path` - the path to fastq.gz files. may be trimmed or not. If the files are trimmed by QC.py, add --trimmed_fq 
flag. See example in the Trimmed data section below. 

#### Provide path to reference sequence:
##### -r|--refseq <refseq/path/> (required)
Don't worry about indexing the fasta file, it happens automatically.

#### Trimmed data:
##### -t|--trimmed_fq (optional)
Run the pipeline with trimmed data (by trimmomatic).

### If you are still confused, here are some examples
To get the usage menu: \
`bash pipeline.sh -h` OR `bash pipeline.sh --help`

To run the whole pipeline on _trimmed_ fastq files found in _/input/path_ and refseq found in _references/refseq.fasta_:\
`bash pipeline.sh -i /input/path/ --refseq references/refseq.fasta --trimmed_fq` \
OR with abbreviated flags: \
`bash pipeline.sh -i /input/path/ -r references/refseq.fasta -t`

### Results:
The consensus sequences are found in CNS/, and further info in results/report.txt

### Directories Hierarchy:
If you choose to use ` bash pipeline.sh -d` to create directories, the hierarchy is as follows:

![alt text](https://github.com/ShebaVirals/covid19/blob/master/dirs_hierarchy.png?raw=true)

### Pipeline's Basic Commands
1. Index the reference sequence. \
`bwa index /refseq/path`

2. For each fastq sample - map to reference and keep as bam file. \
For each fastq.gz file run: \
    `bwa mem /refseq/path sample_R1 sample_R2 | samtools view -Sb - > BAM/sample_name.bam` \
R1 and R2: forward and reverse paired ends.
  
3. Keep only mapped reads of each sample. \
For each bam file run: \
    `samtools view -b -F 260 file_name.bam > BAM/file_name.mapped.bam` \

4. Sort and index each sample in BAM. \
For each mapped bam file run:
    `samtools sort file_name.mapped.bam file_name.mapped.sorted`

5. Create consensus files for each sample. \
For each mapped and sorted bam file run:
    `samtools mpileup -uf /refseq/path file_name.mapped.sorted.bam | bcftools call -mv -Oz -o CNS/sample_name.vcf.gz` \
    `bcftools index CNS/sample` \
    `bcftools consensus -f /refseq/pah CNS/samples_name.vcf.gz > CNS/sample_name.fasta` \
    `rm CNS/sample_name.vcf.gz CNS/sample_name.vcf.gz.csi`

6. Change consensus header name (instead of reference sequence name) \
For eah fasta consensus in CNS run: \
    `sed -i "s/>.*/>${new_name%%.*}/" "$CNS/sample_name.fasta"`

7. Align all consensus sequences and references sequence: \
Gather all fasta sequences (consensus + reference): \
    `cat CNS/*.fasta /refseq/path > alignment/all_not_aligned.fasta` \
Align with mafft and save output in alignment/ directory: \
    `mafft --clustalout alignment/all_not_aligned.fasta > alignment/all_aligned.clustalout` 
    `mafft alignment/all_not_aligned.fasta > alignment/all_aligned.fasta`

8. Report.txt.
    `samtools coverage -H BAM/sample_name.mapped.sorted.bam` for some coverage statistics and number of mapped reads. \
    `samtools view -c BAM/sample_name.bam` for total number of reads in sample. \
    `samtools depth BAM/sample_name.mapped.sorted.bam` for depth of each position in sample \
Check out pipeline.sh code for specifics.
---------------
Dana Bar-Ilan. 
21/07/2020