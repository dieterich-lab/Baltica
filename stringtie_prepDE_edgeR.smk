# -*- coding: utf-8
"""
snakemake pipeline for stringtie to edgeR analysis

"""
configfile: "config.yml"
NAMES = config["samples"].keys()
SAMPLES = config["samples"].values()

rule all:
    input: expand('DEanalysis/{sample}', sample=sample),
    'DEanalysis/sample.tab'

rule create_sample_list:
    output: 'DEanalysis/sample.tab'
    run:

        with open(str(output), 'w') as fou:
            for k, v in zip(NAMES, SAMPLES):
                fou.write('{} {}\n'.format(
                k, v.replace('/Aligned.noS.bam',
                '_StringTieBallgown/stringtieB_outfile.gtf')))

rule prepDE:
    input: 'mappings/{names}.bam'
    output:
        g='DEanalysis/{names}_gene_counts.csv',
        t='DEanalysis/{names}_trans_counts.csv'
    params:
        min_iso = config.get('min_iso', ""),
        min_anchor_len = config.get('min_anchor_len', ""),
        min_junction_cov = config.get('min_junction_cov', ""),
        rf = config.get('rf', "")
    shell:
        '''
        module load python
        python prepDE.py -i {input} --length=75
        '''
