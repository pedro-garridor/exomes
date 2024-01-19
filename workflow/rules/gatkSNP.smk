'''
GATK DNA-Seq Singularity workflow
Copyright (C) 2024, Pedro Garrido Rodríguez

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
'''

# Workflow for germline variant calling with GATK

# https://gatk.broadinstitute.org/hc/en-us/articles/360035535932

rule haplotype_caller:
    input:
        ref=config['ref_fa'],
        bam=config['out_dir'] + 'results/BAM/{sample}.bam',
        bai=config['out_dir'] + 'results/BAM/{sample}.bam.bai',
        cov='resources/SureSelect_Human_AII_Exon_V6_coverage_hg38.bed'
    output:
        vcf=temp(config['out_dir'] + 'results/HC/{sample}.vcf.gz'),
        tbi=temp(config['out_dir'] + 'results/HC/{sample}.vcf.gz.tbi')
    conda: 'envs/gatkcondaenv_4.4.0.0.yml'
    threads: 2
    shell:
        '''
        /opt/gatk-4.4.0.0/gatk HaplotypeCaller \
            -R {input.ref} \
            -I {input.bam} \
            -L {input.cov} \
            -O {output.vcf}
        '''

rule cnn_score_variants:
    # NOTE: 2D mode. For 1D, remove BAM & --tensor-type param
    input:
        vcf=config['out_dir'] + 'results/HC/{sample}.vcf.gz',
        tbi=config['out_dir'] + 'results/HC/{sample}.vcf.gz.tbi',
        ref=config['ref_fa'],
        cov='resources/SureSelect_Human_AII_Exon_V6_coverage_hg38.bed'
    output:
        vcf=temp(config['out_dir'] + 'results/CNN/{sample}.vcf'),
        index=temp(config['out_dir'] + 'results/CNN/{sample}.vcf.idx')
    conda: 'envs/gatkcondaenv_4.4.0.0.yml'
    threads: 2
    shell:
        '''
        /opt/gatk-4.4.0.0/gatk CNNScoreVariants \
            -V {input.vcf} \
            -R {input.ref} \
            -L {input.cov} \
            -O {output.vcf}
        '''

rule filter_variant_tranches:
    input:
        vcf=config['out_dir'] + 'results/CNN/{sample}.vcf',
        index=config['out_dir'] + 'results/CNN/{sample}.vcf.idx',
        dbsnp=config['dbSNP'],
        cov='resources/SureSelect_Human_AII_Exon_V6_coverage_hg38.bed'
    output:
        vcf=temp(config['out_dir'] + 'results/FVT/{sample}.vcf'),
        index=temp(config['out_dir'] + 'results/FVT/{sample}.vcf.idx')
    conda: 'envs/gatkcondaenv_4.4.0.0.yml'
    threads: 2
    shell:
        '''
        /opt/gatk-4.4.0.0/gatk FilterVariantTranches \
            -V {input.vcf} \
            --resource {input.dbsnp} \
            -L {input.cov} \
            --info-key CNN_1D \
            --snp-tranche 99.95 \
            --indel-tranche 99.4 \
            -O {output.vcf}
        '''
