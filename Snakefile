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

#TODO: 
# automatize config.yaml creation, 
# replace 2> with tee

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

SAMPLES=['SRR5876159']

rule all:
    input:
        # FASTQC output for RAW reads
        expand(os.path.join(config['FASTQC'], '{sample}_R{read}_fastqc.zip'), 
               sample = SAMPLES, 
               read = ['1', '2']),

        # Trimming 
        expand(os.path.join(config['TRIMMED'],
                            '{sample}_R{read}_val_{read}.fq.gz'), 
               sample = SAMPLES, 
               read = ['1', '2']),

        # Alignment 
        expand(os.path.join(config['ALIGNMENT'], '{sample}_sorted.bam'), 
               sample = SAMPLES), 

        # Marking Duplicates
        expand(os.path.join(config['ALIGNMENT'], '{sample}_sorted_md.bam'),
               sample = SAMPLES),

        # Filtering 
        expand(os.path.join(config['ALIGNMENT'], 
                            '{sample}_sorted_filtered.bam'),
               sample = SAMPLES),
        expand(os.path.join(config['ALIGNMENT'],
                            '{sample}_sorted_filtered.bam.bai'),
               sample = SAMPLES),

        # multiqc report
        "multiqc_report.html"

    message:
        '\n#################### ATAC-seq pipeline #####################\n'
        'Running all necessary rules to produce complete output.\n'
        '############################################################'

#rule peak_calling:
#
rule filter: 
    """Clean up alignments
    Flags to filter out (-F):
      read unmapped (0x4 = 4)
      mate unmapped (0x8 = 8)
      not primary alignment (0x100 = 256)
      read fails platform/vendor quality checks (0x200 = 512)
      read is PCR or optical duplicate (0x400 = 1024)
      supplementary alignment (0x800 = 2048)
  
    Flags to keep (-f):
      read paired (0x1 = 1)
      read mapped in proper pair (0x2 = 2)
    """
    input:
        bam = os.path.join(config['ALIGNMENT'], '{sample}_sorted_md.bam'),
        index = os.path.join(config['ALIGNMENT'], '{sample}_sorted_md.bam.bai')
    output:
        bam = os.path.join(config['ALIGNMENT'], 
                           '{sample}_sorted_filtered.bam'),
        index = os.path.join(config['ALIGNMENT'],
                             '{sample}_sorted_filtered.bam.bai')
    threads: 8
    message: 
        '\n######################### Filtering ########################\n'
        'Filtering alignment followed by indexing\n'
        '{output.bam}\n'
        '{output.index}\n'
        '############################################################'
    shell:
        'samtools view -b -h -f 3 -F 3852 -@ {threads} {input.bam} > {output.bam};'
        ' samtools index {output.bam}'


rule mark_duplicates:
    input: 
        os.path.join(config['ALIGNMENT'], '{sample}_sorted.bam')
    output:
        bam = os.path.join(config['ALIGNMENT'], '{sample}_sorted_md.bam'),
        index = os.path.join(config['ALIGNMENT'], '{sample}_sorted_md.bam.bai'),
        flagstat_metrics = os.path.join(config['ALIGNMENT_QUAL'], 
                                        'flagstat_{sample}.txt'),
        md_metrics = os.path.join(config['ALIGNMENT_QUAL'], 
                                 'DuplicationMetrics_{sample}.txt')
    params:
        tmp_dir = config['ALIGNMENT'],
        picard_path = config['PICARD']
    message: 
        '\n##################### Marking Duplicates ###################\n'
        'Running Picard MarkDuplicates followed by indexing\n'
        '{output.bam}\n'
        '{output.index}\n'
        '############################################################'
    log:
        os.path.join(config['LOGS'], 'MarkDuplicates_{sample}.log'), 
    shell:
        'java -Xmx4g -Xms4g -jar {params.picard_path} MarkDuplicates' 
        ' I={input} O={output.bam} ASSUME_SORTED=true'
        ' METRICS_FILE={output.md_metrics}'
        ' TMP_DIR={params.tmp_dir} 2> {log};'
        ' samtools index {output.bam};'
        ' samtools flagstat {output.bam} > {output.flagstat_metrics}' 

rule alignment:
    input:
        ref = config['REFERENCE'],
        forw = os.path.join(config['TRIMMED'], '{sample}_R1_val_1.fq.gz'),
        rev = os.path.join(config['TRIMMED'], '{sample}_R2_val_2.fq.gz'),
        amb = config['REFERENCE'] + '.amb',
        ann = config['REFERENCE'] + '.ann',
        bwt = config['REFERENCE'] + '.bwt',
        pac = config['REFERENCE'] + '.pac',
        sa = config['REFERENCE'] + '.sa'
    output:
        os.path.join(config['ALIGNMENT'], '{sample}_sorted.bam')
    params: 
        sort_tmp = os.path.join(config['ALIGNMENT'], '{sample}_sort_tmp')
        #rg = omitting this for now, could be useful if multiple lanes per
        #sample, would be added with -R \'{params.rg}\' to bwa mem
    log:
        bwa = os.path.join(config['LOGS'], 'bwa_{sample}.log'), 
        samtools = os.path.join(config['LOGS'], 'samtools_sort_{sample}.log')
    threads: 4
    message: 
        '\n######################### Mapping ##########################\n'
        'Running bwa_mem follow by sorting to produce:\n'
        '{output}\n'
        '############################################################'
    shell:
        'bwa mem -t {threads} {input.ref} {input.forw} {input.rev}'
        ' 2> {log.bwa} |' 
        ' samtools sort -m 1g -@ {threads} -O bam -T {params.sort_tmp}'
        ' -o {output} - 2> {log.samtools}'

rule build_index:
    input:
        config['REFERENCE']
    output:
        amb = config['REFERENCE'] + '.amb',
        ann = config['REFERENCE'] + '.ann',
        bwt = config['REFERENCE'] + '.bwt',
        pac = config['REFERENCE'] + '.pac',
        sa = config['REFERENCE'] + '.sa'
    message: 
        '\n######################### Indexing ##########################\n'
        'BWA index not found, running bwa index command:\n'
        'bwa index -a bwstw {input}\n'
        '############################################################'
    log:
        os.path.join(config["LOGS"], 'indexing.log')
    shell:
        'bwa index -a bwtsw {input} 2> {log}'

rule trimming:
    input:
        forw = os.path.join(config['READS'], '{sample}_R1.fastq.gz'),
        rev = os.path.join(config['READS'], '{sample}_R2.fastq.gz') 
    output:
        os.path.join(config['TRIMMED'], '{sample}_R1_val_1.fq.gz'),
        os.path.join(config['TRIMMED'], '{sample}_R2_val_2.fq.gz'),
        os.path.join(config['TRIMMED'],
                     '{sample}_R1.fastq.gz_trimming_report.txt'),
        os.path.join(config['TRIMMED'],
                     '{sample}_R2.fastq.gz_trimming_report.txt')
    params:
        qc_outdir = config["FASTQC"],
        outdir = config['TRIMMED']
    log:
        os.path.join(config["LOGS"], 'trim_galore_{sample}.log')
    message:
        '\n######################### Trimming #########################\n'
        'Running trim_galore on:\n'
        '{input.forw}\n{input.rev}\n'
        '############################################################'
    shell:
        'trim_galore --fastqc_args "--outdir {params.qc_outdir}" --gzip'
        ' -o {params.outdir} --paired {input.forw} {input.rev} &> {log}'

rule fastqc:
    input:
        forw = os.path.join(config['READS'], '{sample}_R1.fastq.gz'),
        rev = os.path.join(config['READS'], '{sample}_R2.fastq.gz') 
    output:
        os.path.join(config['FASTQC'], '{sample}_R1_fastqc.zip'),
        os.path.join(config['FASTQC'], '{sample}_R2_fastqc.zip')
    params:
        outdir = config["FASTQC"]
    log:
        os.path.join(config["LOGS"], 'fastqc_{sample}.log')
    shell:
        'fastqc {input.forw} {input.rev} -o {params.outdir} &> {log}'

rule multiqc:
    input:
        # FASTQC output for RAW reads
        expand(os.path.join(config['FASTQC'], 
                            '{sample}_R{read}_fastqc.zip'),
               sample = SAMPLES,
               read = ['1', '2']),

        # Trimming reports
        expand(os.path.join(config['TRIMMED'], 
                            '{sample}_R{read}.fastq.gz_trimming_report.txt'), 
               sample = SAMPLES,
               read = ['1', '2']),

        # Mark Duplicates metrcis
        expand(os.path.join(config['ALIGNMENT_QUAL'], 
                            'DuplicationMetrics_{sample}.txt'),
               sample = SAMPLES),

        # Alignment Stats
        expand(os.path.join(config['ALIGNMENT_QUAL'], 
                            'flagstat_{sample}.txt'),
               sample = SAMPLES)

    output:
        "multiqc_report.html"
    message:
        '\n######################### Multiqc ##########################\n'
        'Running multiqc on all intermediate\n'
        'quality checks.\n'
        '############################################################'
    shell:
        "multiqc --force ."

