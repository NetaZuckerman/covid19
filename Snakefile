import os
from Bio import SeqIO
import subprocess
import collections
import csv

samples=[] # ADD SAMPLES NAMES HERE  * cuoted, without .fastq.gz suffix
# for example: samples=['S1', 'S2', 'S3']
refseq="refs/REF_NC_045512.2.fasta"

def mask_fasta(fasta_record, positions_list):
    """
    Mask fasta file to N in specified positions.
    :param fasta_record: fasta SeqIO record
    :param positions_list: list of positions to replace to 'N'
    :return: the masked fasta
    """
    mutable_seq = fasta_record.seq.tomutable()
    for position in positions_list:
        mutable_seq[position] = 'N'
    # new record, the same with new sequence
    return SeqIO.SeqRecord(mutable_seq, id=fasta_record.id, description=fasta_record.description,
                           name=fasta_record.name, dbxrefs=fasta_record.dbxrefs)


def parse_depth(depth_file, n=1):
    """
    Parse samtools depth file - each position -1 to match str index.
    :param depth_file: Position's dept, as samtools depth command's output.
    :param n: Masking threshold. Every position with depth<n will be changed to 'N'. default:1
    :return: default dict of {<chr>:[pos1,pos2,...] record name with list of positions to mask.
    """
    n_pos = collections.defaultdict(list)  # chr1:[pos1,...], chr2:[...]
    with open(depth_file) as df:
        for line in df:  # line[0]:chr  line[1]:pos line[2]:depth
            line = line.split('\t')
            if int(line[2]) < n:
                n_pos[line[0]].append(int(line[1])-1)  # -1  - depth always starts with 1 instead of 0.
    return n_pos


rule all:
    input:
        "QC/report.txt",
        expand("CNS_5/{sample}_001.fasta", sample=samples),
        "alignment/all_aligned.fasta", "alignment/all_aligned.clustalout"


rule map_to_refseq:
    input:
         fwd="fastq/raw/{sample}_R1_001.fastq.gz", rev="fastq/raw/{sample}_R2_001.fastq.gz"
    output:
         bam="BAM/{sample}_001.bam", mapped_bam="BAM/{sample}_001.mapped.bam"
    params:
        threads='10'
    shell:
         "bwa mem -v 1 -t {params.threads} {refseq} {input.fwd} {input.rev} | samtools view -b -@ {params.threads} - | tee {output.bam} | samtools view -b -@ {params.threads} -F 260 -o {output.mapped_bam} -"


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


rule vcf:
    input:
         sample="BAM/{sample}_001.mapped.sorted.bam", index="BAM/{sample}_001.mapped.sorted.bam.bai",
         depths="QC/depth/{sample}_001.txt"
    output:
          sample="CNS/{sample}_001_calls.vcf.gz",
            index="CNS/{sample}_001_calls.vcf.gz.csi"
    run:
        shell("bcftools mpileup -f {refseq} {input.sample} | bcftools call -mv -Oz -o {output.sample}")
        shell("bcftools index {output.sample}")


rule consensus:
    input:
        samp="CNS/{sample}_001_calls.vcf.gz", index="CNS/{sample}_001_calls.vcf.gz.csi", depth="QC/depth/{sample}_001.txt"
    output:
        cns_1="CNS/{sample}_001.fasta", cns_5="CNS_5/{sample}_001.fasta"
    run:
        shell("bcftools consensus -f {refseq} {input.samp} -o {output.cns_1}")
        # shell("python /home/dana/covid19/mask_fasta.py {output.cns_1} {input.depth} -n 1 -o {output.cns_1}")
        # shell("python /home/dana/covid19/mask_fasta.py {output.cns_1} {input.depth} -n 5 -o {output.cns_5}")
        # cns 1:
        depth_dict = parse_depth(input.depth)
        all_records = [mask_fasta(record, depth_dict[record.id]) for record in SeqIO.parse(output.cns_1, "fasta")]
        SeqIO.write(all_records, output.cns_1, 'fasta')
        # cns 5:
        depth_dict_5 = parse_depth(input.depth, n=5)
        all_records_5 = [mask_fasta(record, depth_dict_5[record.id]) for record in SeqIO.parse(output.cns_1, "fasta")]
        SeqIO.write(all_records_5, output.cns_5, 'fasta')

        os.remove(input.index)
        os.remove(input.samp)
        # change fasta headers from refseq to sample:
        with open(output.cns_1,'r') as fasta1, open(output.cns_5, 'r') as fasta5:
            # CNS
            records=SeqIO.parse(fasta1, 'fasta')
            record = next(records)
            record.id = wildcards.sample + "_001"
            record.description = ""
            # CNS 5
            records_5=SeqIO.parse(fasta5, 'fasta')
            record_5 = next(records_5)
            record_5.id=wildcards.sample + "_001"
            record_5.description = ""
        with open(output.cns_1, 'w') as fasta1, open(output.cns_5, 'w') as fasta5:
            SeqIO.write(record, fasta1, 'fasta')
            SeqIO.write(record_5, fasta5, 'fasta')

rule concat_all_cns:
    input:
        expand("CNS/{sample}_001.fasta", sample=samples) # aggregate all samples together
    output:
        "alignment/all_not_aligned.fasta"
    shell:
        "cat {input} {refseq} > {output}"

rule mafft_alignment:
    input:
        "alignment/all_not_aligned.fasta"
    output:
        clustal="alignment/all_aligned.clustalout", alignment="alignment/all_aligned.fasta"
    run:
        shell("mafft --clustalout {input} > {output.clustal}")
        shell("mafft {input} > {output.alignment}")


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
