'''
GATK DNA-Seq Singularity workflow
Copyright (C) 2024, Pedro Garrido Rodríguez

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
'''

# Workflow for VCF annotation with snpEff

rule snpsift:
    input:
        dbsnp=config['dbSNP'],
        vcf=config['out_dir'] + 'results/filtered/{sample}.vcf.gz',
        tbi=config['out_dir'] + 'results/filtered/{sample}.vcf.gz.tbi'
    output: 'results/rs/{sample}.vcf' # temp('results/rs/{sample}.vcf')
    shell: 'java -Xmx4g -jar /opt/snpEff/SnpSift.jar annotate {input.dbsnp} {input.vcf} > {output}'

rule snpeff:
    input: config['out_dir'] + 'results/rs/{sample}.vcf'
    output: temp(config['out_dir'] + 'results/VCF/{sample}.vcf')
    params:
        genome='hg38'
    threads: 2
    shell:
        # RFE: mark snpEff useless output as temp on snpeff output
        '''
        java -Xmx8g -jar /opt/snpEff/snpEff.jar ann \
            -noStats \
            {params.genome} \
            {input} > {output}
        '''

rule compress_vcf:
    input: config['out_dir'] + 'results/VCF/{sample}.vcf'
    output: protected(config['out_dir'] + 'results/VCF/{sample}.vcf.gz')
    shell: 'bgzip {input}'

rule index_vcf:
    input: config['out_dir'] + 'results/VCF/{sample}.vcf.gz'
    output: protected(config['out_dir'] + 'results/VCF/{sample}.vcf.gz.tbi')
    shell: 'tabix {input}'
