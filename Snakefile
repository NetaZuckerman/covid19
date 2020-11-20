import subprocess
import csv

samples=['42958_S2_L001', 'eg_S1_L001'] # ADD SAMPLES NAMES HERE  * quoted, without .fastq.gz suffix
# for example: samples=['41326_S1_L001', '41357_S2_L001']
refseq="refs/REF_NC_045512.2.fasta"
threads=32
input_path='fastq/raw/'


rule all:
    input:
        "QC/report.txt",
        expand("CNS_5/{sample}_001.fa", sample=samples),
        "alignment/all_aligned.fasta", "alignment/all_aligned.clustalout"

rule map_to_refseq:
    input:
         fwd=input_path+"{sample}_R1_001.fastq.gz", rev=input_path+"{sample}_R2_001.fastq.gz"
    output:
         bam="BAM/{sample}_001.bam", mapped_bam="BAM/{sample}_001.mapped.bam"
    params:
        threads=threads
    shell:
         "bwa mem -v 1 -t {params.threads} {refseq} {input.fwd} {input.rev} | samtools view -@ {params.threads} -b - | tee {output.bam} | samtools view -@ {params.threads} -b -F 260 -o {output.mapped_bam} -"


rule sort_index:
    input:
         "BAM/{sample}_001.mapped.bam"
    output:
        samp="BAM/{sample}_001.mapped.sorted.bam", index="BAM/{sample}_001.mapped.sorted.bam.bai"
    run:
         shell("samtools sort {input} -o {output.samp}")
         shell("samtools index {output.samp}")


rule depth:
    input:
         sample="BAM/{sample}_001.mapped.sorted.bam", index="BAM/{sample}_001.mapped.sorted.bam.bai"
    output:
          "QC/depth/{sample}_001.txt"
    shell:
         "samtools depth -a {input.sample} -o {output}"


rule consensus:
    input:
       "BAM/{sample}_001.mapped.sorted.bam"
    output:
        cns_1="CNS/{sample}_001.fa", cns_5="CNS_5/{sample}_001.fa"
    shell:
        """samtools mpileup -A {input} | ivar consensus -m 1 -p {output.cns_1}
        samtools mpileup -A {input} | ivar consensus -m 5 -p {output.cns_5}"""

rule concat_all_cns:
    input:
        expand("CNS/{sample}_001.fa", sample=samples) # aggregate all samples together
    output:
        "alignment/all_not_aligned.fasta"
    shell:
        "cat {input} {refseq} > {output}"


rule mafft_clustal:
    input:
        "alignment/all_not_aligned.fasta"
    output:
        "alignment/all_aligned.clustalout"
    shell:
        "mafft --clustalout {input} > {output}"


rule mafft_fasta:
    input:
        "alignment/all_not_aligned.fasta"
    output:
        "alignment/all_aligned.fasta"
    shell:
        "mafft {input} > {output}"


rule produce_report:
    input:
         bams=expand('BAM/{samp}_001.bam', samp=samples), mapped_sorted_bams=expand('BAM/{samp}_001.mapped.sorted.bam', samp=samples),
         depths=expand('QC/depth/{samp}_001.txt', samp=samples), indexed_sorted_bams=expand('BAM/{samp}_001.mapped.sorted.bam.bai', samp=samples),
    output:
          report="QC/report.txt"
    run:
        header = ["sample", "mapped%", "mappedreads", "totreads", "covbases", "coverage%", "meandepth",
                                 "maxdepth", "mindepth"]
        with open(output.report, 'w') as report:
            csv_report = csv.DictWriter(report, delimiter='\t', fieldnames=header)
            csv_report.writeheader()
            for sample_bam in input.bams:
                mapped_sorted_bam = sample_bam.replace('.bam', '.mapped.sorted.bam')
                depth = sample_bam.replace('BAM/', 'QC/depth/').replace('.bam', '.txt')
                tot_reads = float(subprocess.Popen(['samtools', 'view', '-c', sample_bam], stdout=subprocess.PIPE).communicate()[0])
                coverage_stats = subprocess.Popen(['samtools', 'coverage', '-H', mapped_sorted_bam], stdout=subprocess.PIPE).communicate()[0].split()
                # print(coverage_stats)
                mapped_reads = float(coverage_stats[3])
                covbases = float(coverage_stats[4])
                percentage_cov = float(coverage_stats[5])
                percentage_mapped = round((mapped_reads / tot_reads)*100, 4)
                with open(depth, 'r') as depth_file:
                    depth_col = [float(row.split('\t')[-1]) for row in depth_file]
                col_no_zeroes=[i for i in depth_col if i>0]
                min_depth = int(min(col_no_zeroes))
                max_depth = int(max(col_no_zeroes))
                mean_depth = round(sum(col_no_zeroes) / len(col_no_zeroes), 4)
                sample_name=sample_bam.replace('BAM/', '').replace('.bam','')

                line = {'sample': sample_name, 'mapped%': percentage_mapped, 'mappedreads': int(mapped_reads),
                        'totreads': int(tot_reads), 'covbases': int(covbases), 'coverage%': percentage_cov,
                        'meandepth': mean_depth, 'maxdepth': max_depth, 'mindepth': min_depth}
                csv_report.writerow(line)
