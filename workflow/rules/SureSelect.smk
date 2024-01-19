rule sure_select_compress:
	input: config['out_dir'] + 'results/FVT/{sample}.vcf'
	output: config['out_dir'] + 'results/FVT/{sample}.vcf.gz' # temp('results/FVT/{sample}.vcf.gz')
	shell: 'bgzip {input}'

rule sure_select_index_1:
	input: config['out_dir'] + 'results/FVT/{sample}.vcf.gz'
	output: config['out_dir'] + 'results/FVT/{sample}.vcf.gz.tbi' # temp('results/FVT/{sample}.vcf.gz.tbi')
	shell: 'tabix {input}'

# SureSelect coverage gathered with UCSC Table Browser: https://genome.ucsc.edu/cgi-bin/hgTables
rule sure_select_filter:
	input:
		cov='resources/SureSelect_Human_AII_Exon_V6_coverage_hg38.bed',
		vcf=config['out_dir'] + 'results/FVT/{sample}.vcf.gz',
		index=config['out_dir'] + 'results/FVT/{sample}.vcf.gz.tbi'
	output: config['out_dir'] + 'results/filtered/{sample}.vcf.gz' # temp('results/filtered/{sample}.vcf.gz')
	shell: '''bcftools filter -R {input.cov} -i 'FILTER="PASS"' {input.vcf} -O z -o {output}'''

rule sure_select_index_2:
	input: config['out_dir'] + 'results/filtered/{sample}.vcf.gz'
	output: config['out_dir'] + 'results/filtered/{sample}.vcf.gz.tbi' # temp('results/filtered/{sample}.vcf.gz.tbi')
	shell: 'tabix {input}'
