# covid19
### CoronaPipeline conda environment
* In order to install all dependencies using conda, run: \
    `conda env create -f <your/covid19/path>/env/CoronaPipeline.yml`
* Make sure to install pangolins conda environment as well as nextstrain, we use pangolin and nextclade (V0) in out script.
for pangolin see https://cov-lineages.org/resources/pangolin/installation.html
for nextclade (V0): https://www.npmjs.com/package/@nextstrain/nextclade
* Make sure you have samtools V 1.10 or above installed.
## Suggested Workflow Using QC.py and pipeline.sh
1. Create project's directories - Optional. \
`python3 QC.py --dirs` \
or \
`bash pipeline.sh -d` \
The directories will be created in the current working directory, and it is recommended to run all following commands 
from that working directory.
* **Fill fastq/raw/ with all raw fastq files.** 
* **Fill refs/ with your desired reference sequence.**
Check out the [directories hierarchy](#directories-hierarchy).

2. Produce reports for each fastq file in fastq/raw (or any other location of your choice) \
`conda activate CoronaPipeline`\
`python3 QC.py --reports -i fastq/raw/`

if trimming is not necessary skip to step 6. 

3. Create a template csv file for convenient trimming. See details at Template. \
`python3 QC.py --template -i fastq/raw/ ` \
Take a look at _template.csv_ produced as output and fix some trimming parameters according to fastqc/multiqc reports.
Scroll down to read more about the csv file.  

4. Trim fastq files according to the csv file created in (2). \
`python3 QC.py --trim template.csv -i /input/path/to/fastq/files -o /trimmed/path` 

5. Produce second reports to make sure you're happy with the trimming results. \
`python3 QC.py --reports -i /trimmed/path -o /output/fastqc/2/path` \
Be sure to provide the path to the trimmed fastq files this time. Repeat steps 2-5 if more trimming is needed.
   

6. Run pipeline to produce consensus sequences for your samples, and get additional info such as number of reads, depth, 
coverage, etc. in results/report.txt. Run the pipeline on the trimmed fastq files (the pipeline uses 
the paired files, ignoring singletons). \
-i: input fastq files path \
   * **Be sure to deactivate environment before running the pipeline.sh!** Conda is automatically activated in the pipeline as well as pangolins environment.
    
   `conda deactivate`\
`bash pipeline.sh -i fastq/trimmed/ -r /refs/refseq.fa` 



## Additional info: QC.py
### Usage:
`python QC.py [-h] [-t [CSV] | --template | -r] [-i WD] [-o OUTPUT_PATH]`
#### -r | --reports : QC Reports
Produces fastqc and multiqc reports for all fastq files. \
_-i_ : path to fastq files locations. default: fastq/raw/ \
_-o_ : path to drop the template.csv. default: QC/ \
`python3 QC.py -r -i path/to/fastq/files -o path/for/output/reports` \
To create QC reports with default input and output arguments:
`python3 QC.py -r` 

#### --template : Template of trimmomatic command for each sample
 _-i_ : path to fastq files location. default: fastq/raw/ \
 _-o_ : path to drop the template.csv. default: working directory \
`python3 QC.py --template -i some/path/to/fastq/location/ -o some/path/for/output` \
To create template with default input and output arguments: \
`python3 QC.py --template` 

#### -t|--trim <path/to/csv> <optional-path/to/base/dir> : Trim fastq files with trimmomatic, according to csv file.
_-i_ : path to fastq files location. default: fastq/raw/ \
_-o_ : path to drop the trimmed samples. default: fastq/trimmed \
`python3 QC.py -t path/to/file.csv -i path/to/raw/fastq/files -o path/to/trimmed/fatsq/output` \
To trim fastq files with default input and output arguments: \
`python3 QC.py -t path/to/file.csv ` or
`python3 QC.py --trim path/to/file.csv`\
To trim with the template csv as input, simply run `python3 QC.py -t template.csv`.

### The CSV file:
Template csv header:
>ends, input_forward, input_reverse, phred, threads, ILLUMINACLIP:, SLIDINGWINDOW:, LEADING:, TRAILING:, CROP:,
>HEADCROP:, MINLEN: 

Each line will represent a fastq file.

**Required fields**: ends, input_forward. Make sure those fields are filled in the csv file \
trimmomatic requires output files as well, but those are given automatically in the script. Don't include them in the csv.
* ends`[PE|SE]` - PE:paired end, SE: single end. If you choose `PE`, be sure to include both input_forward **and**  input_reverse! 
* input_forward `<sample_R1.fastq>`: sample name
* input_reverse `<sample_R2.fastq>`: sample name

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
# Additional Info: pipeline.sh
### Usage:
 `bash pipeline.sh [-h | -r [REFSEQ_PATH] [-i WD] [-o OUTPUT_PATH]`

### Run pipeline on samples to get consensus sequences and additional data:
#### -i /input/path : provide input path  (required)
`/input/path` - the path to fastq files.

#### -r|--refseq <refseq/path/> : provide refseq path (required)
Don't worry about indexing the fasta file, it happens automatically.

#### --spike (optional)
run the pipeline for hi-spike sequencing. \
 `bash pipeline.sh  -r path/to/spike_ref.fa -i path/to/fastq/location -o path/for/output/reports --spike `

#### -n|--newNextclade (optional)
use nextclade version 1.3.0. \
 `bash pipeline.sh  -r path/to/ref.fa -i path/to/fastq/location -o path/for/output/reports --newNextclade `
    
#### -p|--processes <int> (optional)
number of processes (samples) to run in parallel. default: 1. \
 `bash pipeline.sh  -r path/to/ref.fa -i path/to/fastq/location -o path/for/output/reports -p 5 `
    
#### -t|--threads <int> (optional)
number of threads for each sample. default: 32. \
 `bash pipeline.sh  -r path/to/ref.fa -i path/to/fastq/location -o path/for/output/reports --threads 32 `
        

### If you are still confused, here are some examples
To get the usage menu: \
`bash pipeline.sh -h` OR `bash pipeline.sh --help`

### Results:
The consensus sequences are found in CNS/ and CNS_5/, and further info in results/

### Pipeline's Basic Steps
`conda activate CoronaPipeline`
1. Index the reference sequence. \
`bwa index /refseq/path`

2. For each fastq sample - map to reference and keep as bam file. \
For each fastq file run: \
    `bwa mem /refseq/path sample_R1 sample_R2 | samtools view -b - > BAM/sample_name.bam` \
R1 and R2: forward and reverse paired ends.
  
3. Keep only mapped reads of each sample. \
For each bam file run: \
    `samtools view -b -F 260 file_name.bam > BAM/file_name.mapped.bam` 

4. Sort and index each sample in BAM. \
For each mapped bam file run:
    `samtools sort BAM/file_name.mapped.bam BAM/file_name.mapped.sorted.bam`

5. Create consensus files for each sample. \
For each mapped and sorted bam file run:
    `samtools mpileup -A  BAM/file_name.mapped.sorted.bam | ivar consensus -m 5 -p CNS5/file_name`

6. Align all consensus sequences and references sequence: \
Gather all fasta sequences (consensus + reference): \
    `cat CNS5/*.fasta /refseq/path > alignment/all_not_aligned.fasta` \
Align with augur (mafft based) and save output in alignment/ directory: \
  ` augur align \
  --sequences alignment/all_not_aligned.fasta \
  --reference-sequence /refseq/path \
  --output alignment/all_aligned.fasta`

7. Classify mutations with pangolin: \
`conda deactivate` \
`conda activate pangolin` \
   `pangolin alignment/all_not_aligned.fasta --outfile results/pangolinClades.csv` \
   `conda deactivate`

8. Check for mutations and variants: \
    `conda activate CoronaPipeline` \
    `python MutTable.py alignment/all_aligned.fasta results/nuc_muttable.xlsx` \
    `python translated_table.py alignment/all_aligned.fasta results/AA_muttable.xlsx regions.csv` \
    `python variants.py alignment/all_aligned.fasta results/variants.csv results/pangolinClades.csv`
    

9. Report.txt: \
    some coverage statistics and number of mapped reads: `samtools coverage -H BAM/sample_name.mapped.sorted.bam`  \
    total number of reads in sample: `samtools view -c BAM/sample_name.bam` \
    depth of each position in sample: `samtools depth BAM/sample_name.mapped.sorted.bam`  \
Check out pipeline.sh code for specifics.


## Minipipeline
### Script for analysis of fasta files (instead of fastq)
As in the regular pipeline, make sure to set the CoronaPipeline conda environment. No need to actiavte it 
before the run. \
### Usage:
`bash /path/to/minipipeline.sh -i fasta_path.fa -r covid_refseq.fasta` 
#### --dontAlign : run the minipipeline with aligned fasta input. (optional)
The first sequence in the alignment must be the reference sequence \
`bash minipipeline.sh -i alignment.fasta --dontAlign` 
#### --invr : generate INVR report containing information about the mutations. (optional)
`bash minipipeline.sh -i alignment.fasta -r ref.fasta --invr` 

---------------
Dana Bar-Ilan. \
last update: 20.02.2022 \
last updated by: Hagar Morad
