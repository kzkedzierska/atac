"""
Author: Katarzyna Kedzierska
Affiliation: WCHG, Unieversity of Oxford
Aim: A simple Snakemake workflow for processing ATAC-seq data.
Date: Fri 12 Apr - 
Run: snakemake -s Snakefile
"""
from glob import glob
import re
import os 

configfile: "config.yaml"

# Helper functions
def fastq_to_bn(fastq_bn):
    return re.sub(r".fastq[.gz]*", "", fastq_bn)

# FASTQ files
READS = glob(config['READS'] + "/*.fastq[.gz]*") 
SAMPLES = list(set([os.path.basename(read_file).split("_")[0] for read_file in READS]))
FASTQ_BN = [fastq_to_bn(os.path.basename(fastq_file)) for fastq_file in READS]



rule all:
    input:
        # FASTQC output for RAW reads
        expand(config['FASTQC_DIR'] + '{sample}_R{read}_fastqc.zip', 
               sample = SAMPLES, 
               read = ['1', '2']),

        # Trimming 
        expand(config['TRIMMED'] + '{sample}_R{read}_trimmed.fq.gz', 
               sample = SAMPLES, 
               read = ['1', '2']),

        # multiqc report
        "multiqc_report.html"

    message:
        '\n#################### ATAC-seq pipeline #####################\n'
        'Running all necessary rules to produce complete output.\n'
        '############################################################'

rule trimming:
    input:
        forw = config['READS'] + '{sample}_R1.fastq.gz',
        rev = config['READS'] + '{sample}_R2.fastq.gz' 
    output:
        config['TRIMMED'] + '{sample}_R1_trimmed.fq.gz',
        config['TRIMMED'] + '{sample}_R2_trimmed.fq.gz',
        config['TRIMMED'] + '{sample}_R1_fastq.gz_trimming_report.txt',
        config['TRIMMED'] + '{sample}_R2_fastq.gz_trimming_report.txt'
    params:
        qc_outdir = config["FASTQC_DIR"],
        outdir = config['TRIMMED']
    log:
        config["LOG_DIR"] + 'trim_galore_{sample}.log' 
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
        config['FASTQC_DIR'] + '{sample}_R1_fastqc.zip',
        config['FASTQC_DIR'] + '{sample}_R2_fastqc.zip'
    params:
        outdir = config["FASTQC_DIR"]
    log:
        config["LOG_DIR"] + 'fastqc_{sample}.log'
    shell:
        'fastqc {input.forw} {input.rev} -o {params.outdir} &> {log}'

rule multiqc:
    input:
        # FASTQC output for RAW reads
        expand(config['FASTQC_DIR'] + '{sample}_R{read}_fastqc.zip', 
               sample = SAMPLES, 
               read = ['1', '2']),

        # Trimming reports
        expand(config['TRIMMED'] + '{sample}_R{read}_fastq.gz_trimming_report.txt', 
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
