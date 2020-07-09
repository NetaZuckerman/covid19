# covid19

QC.py [--trim [str]| --template | --reports] for automated QC. 

>_--trim [csv file path]:_ Input is the csv file path. The file contains a constant header (please don't change it) and each line represents a diferent sample.
 Samples are found in fastq/raw, trimmed samples are found in fastq/trimmed. Make sure to run the program in the right working dir, and that all required directories exist (that might be taken care of in the future.)

---
For more info:
>Trimmomatic 0.39 manual: http://www.usadellab.org/cms/?page=trimmomatic



Dana Bar-Ilan.

last modified: 09.07.2020