"""
Author: Katarzyna Kedzierska
Affiliation: WCHG, Unieversity of Oxford
Aim: A simple Snakemake workflow for processing ATAC-seq data.
Date: Fri 12 Apr - 
Run: snakemake -s Snakefile

While preparing this pipeline I used the following file as a reference,
especially for Snakemake syntax as this is my first snakemake pipeline.
https://github.com/porchard/ATACseq-Snakemake/blob/master/src/Snakefile
"""

#TODO: replace + with os.path.join!

from glob import glob
import re
import os 

configfile: "config.yaml"

# Helper functions
def fastq_to_bn(fastq_bn):
    return re.sub(r".fastq[.gz]*", "", fastq_bn)

# FASTQ files
READS = glob(config['READS'] + "/*.fastq[.gz]*") 
#SAMPLES = list(set([os.path.basename(read_file).split("_")[0] for read_file in READS]))
#FASTQ_BN = [fastq_to_bn(os.path.basename(fastq_file)) for fastq_file in READS]

SAMPLES=['SRR5876159', 'SRR5876662']

rule all:
    input:
        # FASTQC output for RAW reads
        expand(config['FASTQC'] + '{sample}_R{read}_fastqc.zip', 
               sample = SAMPLES, 
               read = ['1', '2']),

        # Trimming 
        expand(config['TRIMMED'] + '{sample}_R{read}_val_{read}.fq.gz', 
               sample = SAMPLES, 
               read = ['1', '2']),
        # Alignment 
        expand(os.path.join(config['ALIGNMENT'], '{sample}_sorted.bam'), 
               sample = SAMPLES), 

        # multiqc report
        "multiqc_report.html"

#rule peak_calling:
#
#rule mark_duplicates:
#    input:
#    
#    output:
#    shell:
#        
#
#rule filter:
#    input:
#    shell
#

    message:
        '\n#################### ATAC-seq pipeline #####################\n'
        'Running all necessary rules to produce complete output.\n'
        '############################################################'

rule alignment:
    input:
        config['REFERENCE'],
        os.path.join(config['TRIMMED'], '{sample}_R1_val_1.fq.gz'),
        os.path.join(config['TRIMMED'], '{sample}_R2_val_2.fq.gz')
    output:
        os.path.join(config['ALIGNMENT'], '{sample}_sorted.bam')
    params: 
        sort_tmp = os.path.join(config['ALIGNMENT'], '{sample}_sort_tmp')
        #rg = omitting this for now, could be useful if multiple lanes per
        #sample, would be added with -R \'{params.rg}\' to bwa mem
    log:
        bwa = os.path.join(config['LOGS'], 'bwa_{sample}.log'), 
        samtools = os.path.join(config['LOGS'], 'samtools_sort_{sample}.log')
    message: 
        '\n######################### Mapping ##########################\n'
        'Running bwa_mem follow by sorting on {sample}\n'
        '############################################################'
    shell:
        'bwa mem -t {threads} {input} 2> {log.bwa} |' 
        ' samtools sort -m 1g -@ {threads} -O bam -T {params.sort_tmp}'
        ' -o {output} - 2> {log.samtools}'

rule trimming:
    input:
        forw = config['READS'] + '{sample}_R1.fastq.gz',
        rev = config['READS'] + '{sample}_R2.fastq.gz' 
    output:
        config['TRIMMED'] + '{sample}_R1_val_1.fq.gz',
        config['TRIMMED'] + '{sample}_R2_val_2.fq.gz',
        config['TRIMMED'] + '{sample}_R1.fastq.gz_trimming_report.txt',
        config['TRIMMED'] + '{sample}_R2.fastq.gz_trimming_report.txt'
    params:
        qc_outdir = config["FASTQC"],
        outdir = config['TRIMMED']
    log:
        config["LOGS"] + 'trim_galore_{sample}.log' 
    message:
        '\n######################### Trimming #########################\n'
        'Running trim_galore on:\n'
        '{input.forw}\n{input.rev}\n'
        '############################################################'
    shell:
        'trim_galore --fastqc_args "--outdir {params.qc_outdir}"'
        ' --gzip -o {params.outdir} --paired {input.forw} {input.rev} &> {log}'

rule fastqc:
    input:
        forw = config['READS'] + '{sample}_R1.fastq.gz',
        rev = config['READS'] + '{sample}_R2.fastq.gz' 
    output:
        config['FASTQC'] + '{sample}_R1_fastqc.zip',
        config['FASTQC'] + '{sample}_R2_fastqc.zip'
    params:
        outdir = config["FASTQC"]
    log:
        config["LOGS"] + 'fastqc_{sample}.log'
    shell:
        'fastqc {input.forw} {input.rev} -o {params.outdir} &> {log}'

rule multiqc:
    input:
        # FASTQC output for RAW reads
        expand(config['FASTQC'] + '{sample}_R{read}_fastqc.zip', 
               sample = SAMPLES, 
               read = ['1', '2']),

        # Trimming reports
        expand(config['TRIMMED'] + '{sample}_R{read}.fastq.gz_trimming_report.txt', 
               sample = SAMPLES,
               read = ['1', '2'])
    output:
        "multiqc_report.html"
    message:
        '\n######################### Multiqc ##########################\n'
        'Running multiqc on all intermediate\n'
        'quality checks.\n'
        '############################################################'
    shell:
        "multiqc --force ."

