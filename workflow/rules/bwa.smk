'''
GATK DNA-Seq Singularity workflow
Copyright (C) 2024, Pedro Garrido Rodríguez

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
'''

# Workflow for read realignment using BWA

# https://gatk.broadinstitute.org/hc/en-us/articles/360035535912

rule bwa_mem:
    # Illumina/454/IonTorrent < 70bp:
        # bwa aln ref.fa reads.fq > reads.sai; bwa samse ref.fa reads.sai reads.fq > aln-se.sam
    # Illumina/454/IonTorrent > 70bp:
        # bwa mem ref.fa read1.fq read2.fq > aln-pe.sam
    # Illumina paired-end < 70bp:
        # bwa aln ref.fa read1.fq > read1.sai; bwa aln ref.fa read2.fq > read2.sai
        # bwa sampe ref.fa read1.sai read2.sai read1.fq read2.fq > aln-pe.sam
    # PacBio subreads or ONT:
        # bwa mem -x pacbio ref.fa reads.fq > aln.sam
        # bwa mem -x ont2d ref.fa reads.fq > aln.sam
    input:
        ref=config['ref_fa'],
        reads=config['input_dir'] + '{sample}'
    output: temp(config['out_dir'] + 'results/BAM/{sample}.sam')
    threads: workflow.cores
    shell: 'bwa mem -t {threads} {input.ref} $(find {input.reads} -type f | sort) > {output}'

rule samtools_sort:
    input: config['out_dir'] + 'results/BAM/{sample}.sam'
    output: temp(config['out_dir'] + 'results/sorted/{sample}.bam')
    shell: 'samtools sort -@ {threads} -o {output} {input}'

rule add_read_groups:
    input: config['out_dir'] + 'results/sorted/{sample}.bam'
    output: temp(config['out_dir'] + 'results/RG/{sample}.bam')
    shell:
        '''
        picard-tools AddOrReplaceReadGroups \
            -I {input} \
            -O {output} \
            -RGLB lib1 \
            -RGPL ILLUMINA \
            -RGPU unit1 \
            -RGSM {wildcards.sample}
        '''

rule mark_duplicates:
    # NOTE: DO NOT RUN mark_duplicates if sequencing is amplicon-based
    # NOTE: if not running this rule, modify base_recalibrator & apply_recalibration input
    input: 'results/RG/{sample}.bam'
    output:
        bam=temp(config['out_dir'] + 'results/duplicates/{sample}.bam'),
        metrics=temp(config['out_dir'] + 'results/duplicates/{sample}.txt')
    conda: 'envs/gatkcondaenv_4.4.0.0.yml'
    shell:
        '''
        picard-tools MarkDuplicates \
            -I {input} \
            -O {output.bam} \
            -M {output.metrics}
        '''

rule base_recalibrator:
    input:
        ref=config['ref_fa'],
        bam=config['out_dir'] + 'results/duplicates/{sample}.bam',
        dbsnp=config['dbSNP']
    output: temp(config['out_dir'] + 'results/BAM/{sample}.recalibrator.table')
    conda: 'envs/gatkcondaenv_4.4.0.0.yml'
    params:
        max_cycle=600
    shell:
        '''
        /opt/gatk-4.4.0.0/gatk BaseRecalibrator \
            -I {input.bam} \
            -R {input.ref} \
            --known-sites {input.dbsnp} \
            -max-cycle {params} \
            -O {output}
        '''

rule apply_recalibration:
    input:
        ref=config['ref_fa'],
        bam=config['out_dir'] + 'results/duplicates/{sample}.bam',
        table=config['out_dir'] + 'results/BAM/{sample}.recalibrator.table'
    output:
        bam=protected(config['out_dir'] + 'results/BAM/{sample}.bam'),
        bai=temp(config['out_dir'] + 'results/BAM/{sample}.bai')
    conda: 'envs/gatkcondaenv_4.4.0.0.yml'
    shell:
        '''
        /opt/gatk-4.4.0.0/gatk ApplyBQSR \
            -R {input.ref} \
            -I {input.bam} \
            --bqsr-recal-file {input.table} \
            -O {output.bam}
        '''

rule index_bam:
    input: config['out_dir'] + 'results/BAM/{sample}.bam'
    output: protected(config['out_dir'] + 'results/BAM/{sample}.bam.bai')
    shell: 'samtools index {input} > {output}'
