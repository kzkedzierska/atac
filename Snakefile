"""
Author: Katarzyna Kedzierska
Affiliation: WCHG, Unieversity of Oxford
Aim: A simple Snakemake workflow for processing ATAC-seq data.
Date: February 2020
Run: snakemake -s Snakefile

While preparing the first draft of this pipeline I used the following as a ref,
especially for Snakemake syntax as this was my first snakemake pipeline.
https://github.com/porchard/ATACseq-Snakemake/blob/master/src/Snakefile
"""

# LIMITATIONS:
# assumes that initial input files are named in a certain way:
# {sample_id}_R1.fastq.gz; that {sample_id} should match what's
# in the 1st column of the SAMPLE_SHEET file.

#TODO:
# - automatize config.yaml creation
# - keeping MACS2 for now to be able to compare with the TCGA paper
# - figure out how to handle memory in bwa aligment (it's gready!)
# - normalize this naming convention - the name switch at alignment 
# step is problematic

from glob import glob
import re
import os
# import logging as logger

configfile: "config.yaml"
SAMPLES_DICT = dict()

with open(config['SAMPLE_SHEET'], "r+") as fil:
    next(fil)
    for lin in fil.readlines():
        row = lin.strip("\n").split("\t")
        sample_id = row[0]
        sample_name = row[1]
        if sample_name in SAMPLES_DICT.keys():
            SAMPLES_DICT[sample_name].append(sample_id)
        else:
            SAMPLES_DICT[sample_name] = [sample_id]

SAMPLES = list(SAMPLES_DICT.keys())
SAMPLE_IDS = [sample_id for sample in SAMPLES_DICT.values() for sample_id in sample]

print(SAMPLE_IDS)

rule all:
    input:
        # Filtered alignment
        expand(os.path.join(config['FILTERED'],
                            '{sample_id}.bam{ext}'),
               sample_id = SAMPLE_IDS,
               ext = ['', '.bai']),
        # Called peaks
        expand(os.path.join(config['HMMRATAC'], '{sample_id}{ext}'),
               sample_id = SAMPLE_IDS,
               ext = ['.model', '_peaks.gappedPeak',
                      '_summits.bed', '.bedgraph', '.log']),
        expand(os.path.join(config['MACS2'], '{sample_id}{ext}'),
               sample_id = SAMPLE_IDS,
               ext = ['_peaks.xls', '_peaks.narrowPeak',
                      '_summits.bed', '_treat_pileup.bdg']),
        # multiqc report
        "multiqc_report.html"

    message:
        '\n#################### ATAC-seq pipeline #######################\n'
        'Running all necessary rules to produce complete output.\n'
        '##############################################################'

# logger.info(SAMPLES)

rule hmmratac:
    input:
        bam = os.path.join(config['FILTERED'], '{sample_id}.bam'),
        index = os.path.join(config['FILTERED'], '{sample_id}.bam.bai')
    output:
        model = os.path.join(config['HMMRATAC'], '{sample_id}.model'),
        gappedPeak = os.path.join(config['HMMRATAC'], '{sample_id}_peaks.gappedPeak'),
        summits = os.path.join(config['HMMRATAC'], '{sample_id}_summits.bed'),
        states = os.path.join(config['HMMRATAC'], '{sample_id}.bedgraph'),
        logs = os.path.join(config['HMMRATAC'], '{sample_id}.log')
    log:
        os.path.join(config['LOGS'], 'hmmratac', '{sample_id}.log')
    params:
        genomes = config['GENOMES'],
        blacklisted = config['BLACKLIST'],
        main_dir = config['HMMRATAC'],
        dir = os.path.join(config['HMMRATAC'], '{sample_id}')
    resources:
        mem_mb = 42000
    threads: 4
    message:
        '\n######################### Peak calling #######################\n'
        'Peak calling with hmmratac for {output.model}.\n'
        '##############################################################'
    conda:
        config['CONDA_HMMRATAC']
    shell:
        'if [ ! -d {params.main_dir} ]; then mkdir -p {params.main_dir}; fi; '
        'HMMRATAC -Xms2g -Xmx{resources.mem_mb}m '
        '--bam {input.bam} --index {input.index} '
        '--genome {params.genomes} --blacklist {params.blacklisted} '
        '--output {params.dir} --bedgraph true &> {log}'

rule macs2:
   input:
       os.path.join(config['FILTERED'], '{sample_id}.bam')
   output:
       peaks = os.path.join(config['MACS2'], '{sample_id}_peaks.xls'),
       narrowPeak = os.path.join(config['MACS2'], '{sample_id}_peaks.narrowPeak'),
       summits = os.path.join(config['MACS2'], '{sample_id}_summits.bed'),
       pileup = os.path.join(config['MACS2'], '{sample_id}_treat_pileup.bdg')
   params:
       qval = 0.05,
       outdir = config['MACS2']
   log:
       os.path.join(config['LOGS'], 'macs2', '{sample_id}.log')
   message:
       '\n######################### Peak calling #######################\n'
       'Peak calling with MACS2 for {output.peaks}\n'
       '##############################################################'
   conda:
       config['CONDA_MACS2']
   shell:
       'macs2 callpeak --gsize hs --qvalue {params.qval} --format BAM '
       '--keep-dup all --call-summits --bdg --verbose 3 '
       '--extsize 150 --shift 75 --nomodel --nolambda '
       '--treatment {input} --outdir {params.outdir} --name {wildcards.sample_id} '
       '&> {log}'

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

    Also, filtering out blackilsted regions.
    """
    input:
        bam = os.path.join(config['ALIGNMENT'], '{sample_id}_sorted_md.bam'),
        index = os.path.join(config['ALIGNMENT'], '{sample_id}_sorted_md.bam.bai'),
        blacklist = config['BLACKLIST']
    output:
        bam = os.path.join(config['FILTERED'],
                           '{sample_id}.bam'),
        index = os.path.join(config['FILTERED'],
                             '{sample_id}.bam.bai')
    conda:
        config['CONDA_ALIGNMENT']
    threads: 8
    message:
        '\n########################## Filtering #########################\n'
        'Filtering alignment followed by indexing\n'
        '{output.bam}\n'
        '{output.index}\n'
        '##############################################################'
    shell:
        'samtools view -b -h -f 3 -F 3852 -@ {threads} {input.bam} '
        '$(for i in $(echo $(seq 22) X); do echo chr$i; done | xargs) | '
        'bedtools intersect -v -abam stdin -b {input.blacklist} > {output.bam}; '
        'samtools index {output.bam}'

rule mark_duplicates:
    input:
        os.path.join(config['ALIGNMENT'], '{sample_id}_sorted.bam')
    output:
        bam = os.path.join(config['ALIGNMENT'], '{sample_id}_sorted_md.bam'),
        index = os.path.join(config['ALIGNMENT'], '{sample_id}_sorted_md.bam.bai'),
        flagstat_metrics = os.path.join(config['ALIGNMENT_QUAL'],
                                        'flagstat_{sample_id}.txt'),
        md_metrics = os.path.join(config['ALIGNMENT_QUAL'],
                                 'DuplicationMetrics_{sample_id}.txt')
    params:
        tmp_dir = config['ALIGNMENT']
    message:
        '\n##################### Marking Duplicates ###################\n'
        'Running Picard MarkDuplicates followed by indexing\n'
        '{output.bam}\n'
        '{output.index}\n'
        '############################################################'
    log:
        os.path.join(config['LOGS'], 'MarkDuplicates', '{sample_id}.log')
    resources:
        mem_mb = 4000
    conda:
        config['CONDA_ALIGNMENT']
    shell:
        'picard -Xmx{resources.mem_mb}m MarkDuplicates'
        ' I={input} O={output.bam} ASSUME_SORTED=true'
        ' METRICS_FILE={output.md_metrics}'
        ' TMP_DIR={params.tmp_dir} 2> {log};'
        ' samtools index {output.bam};'
        ' samtools flagstat {output.bam} > {output.flagstat_metrics}'

rule alignment:
    input:
        ref = config['REFERENCE'],
        forw = os.path.join(config['TRIMMED'], '{sample_id}_R1_val_1.fq.gz'),
        rev = os.path.join(config['TRIMMED'], '{sample_id}_R2_val_2.fq.gz'),
        amb = config['REFERENCE'] + '.amb',
        ann = config['REFERENCE'] + '.ann',
        bwt = config['REFERENCE'] + '.bwt',
        pac = config['REFERENCE'] + '.pac',
        sa = config['REFERENCE'] + '.sa'
    output:
        os.path.join(config['ALIGNMENT'], '{sample_id}_sorted.bam')
    params:
        sort_tmp = os.path.join(config['ALIGNMENT'], '{sample_id}_sort_tmp'),
        platform = config['SEQUENCING_PLATFORM'],
        seq_cen = config['SEQUENCING_CENTRE'],
        sample_sheet = config['SAMPLE_SHEET']
    log:
        bwa = os.path.join(config['LOGS'], 'bwa', '{sample_id}.log'),
        samtools = os.path.join(config['LOGS'], 'samtools', 'sort_{sample_id}.log')
    conda:
        config['CONDA_ALIGNMENT']
    threads: 4
    message:
        '\n######################### Mapping ##########################\n'
        'Running bwa_mem follow by sorting to produce:\n'
        '{output}\n'
        '############################################################'
    shell:
        'orgnl_id="$(basename {input.forw} | sed \'s/_[R]*[12]*_val_[R]*[12]*.f[ast]*q.gz//g\')";'
        ' sample_name=$(grep "$orgnl_id" {params.sample_sheet} | cut -f2); '
        ' bwa mem -t {threads} {input.ref} {input.forw} {input.rev}'
        ' -R "@RG\\tID:$orgnl_id\\tSM:$sample_name\\tCN:{params.seq_cen}\\tPL:{params.platform}"'
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
    conda:
        config['CONDA_ALIGNMENT']
    shell:
        'bwa index -a bwtsw {input} 2> {log}'

rule trimming:
    input:
        os.path.join(config['READS'], '{sample_id}_R1.fastq.gz'),
        os.path.join(config['READS'], '{sample_id}_R2.fastq.gz')
    output:
        os.path.join(config['TRIMMED'], '{sample_id}_R1_val_1.fq.gz'),
        os.path.join(config['TRIMMED'], '{sample_id}_R2_val_2.fq.gz'),
        os.path.join(config['TRIMMED'],
                     '{sample_id}_R1.fastq.gz_trimming_report.txt'),
        os.path.join(config['TRIMMED'],
                     '{sample_id}_R2.fastq.gz_trimming_report.txt')
    params:
        qc_outdir = config["FASTQC"],
        outdir = config['TRIMMED']
    conda:
        config['CONDA_QUALITY']
    log:
        os.path.join(config['LOGS'], 'trim_galore', '{sample_id}.log')
    message:
        '\n######################### Trimming #########################\n'
        'Running trim_galore on:\n'
        '{input}\n'
        '############################################################'
    shell:
        'trim_galore --fastqc_args "--outdir {params.qc_outdir}" --gzip '
        '-o {params.outdir} --paired {input} &> {log}'

rule fastqc:
    input:
        os.path.join(config['READS'], '{sample_id}_R1.fastq.gz'),
        os.path.join(config['READS'], '{sample_id}_R2.fastq.gz')
    output:
        os.path.join(config['FASTQC'], '{sample_id}_R1_fastqc.zip'),
        os.path.join(config['FASTQC'], '{sample_id}_R2_fastqc.zip')
    params:
        outdir = config['FASTQC']
    conda:
        config['CONDA_QUALITY']
    log:
        os.path.join(config['LOGS'], 'fastqc', '{sample_id}.log')
    shell:
        'fastqc {input} -o {params.outdir} &> {log}'

rule multiqc:
    input:
        # FASTQC output for RAW reads
        expand(os.path.join(config['FASTQC'],
                            '{sample_id}_R{read}_fastqc.zip'),
               sample_id = SAMPLE_IDS,
               read = ['1', '2']),

        # Trimming reports
        expand(os.path.join(config['TRIMMED'],
                            '{sample_id}_R{read}.fastq.gz_trimming_report.txt'),
               sample_id = SAMPLE_IDS,
               read = ['1', '2']),

        # Mark Duplicates metrcis
        expand(os.path.join(config['ALIGNMENT_QUAL'],
                            'DuplicationMetrics_{sample_id}.txt'),
               sample_id = SAMPLE_IDS),

        # Alignment Stats
        expand(os.path.join(config['ALIGNMENT_QUAL'],
                            'flagstat_{sample_id}.txt'),
               sample_id = SAMPLE_IDS)

    output:
        "multiqc_report.html"
    conda:
        config['CONDA_QUALITY']
    message:
        '\n######################### Multiqc ##########################\n'
        'Running multiqc on all intermediate\n'
        'quality checks.\n'
        '############################################################'
    shell:
        "multiqc --force ."
