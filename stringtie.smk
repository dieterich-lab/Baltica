# -*- coding: utf-8
"""
Created on 17:07 27/02/2018 
Snakemake file for stringtie.
.. usage:
    snakemake -s stringtie.smk --keep-going \
    --cluster 'sbatch --job-name stringtie_pipeline' \
    --cluster-config  cluster.json --jobs 100
    
"""
__author__ = "Thiago Britto Borges"
__copyright__ = "Copyright 2018, Dieterichlab"
__email__ = "Thiago.BrittoBorges@uni-heidelberg.de"
__license__ = "MIT"

import os
configfile: "config.yml"
NAMES = config["samples"].keys()
SAMPLES = config["samples"].values()

rule all:
    input: expand('stringtie/{names}.gtf', names=NAMES)

rule stringtie:
    input: 'mappings/{names}.bam'
    output: 'stringtie/{names}.gtf'
    threads: 8 
    params:
        min_iso = config.get('min_iso', ""),
        min_anchor_len = config.get('min_anchor_len', ""),
        min_junction_cov = config.get('min_junction_cov', ""),
        rf = config.get('rf', "")
    shell:
        '''
        module load stringtie
        stringtie {input} -v -p {threads} -o {output} \
        {params.rf} {params.min_iso} {params.min_anchor_len} \
        {params.min_junction_cov}
        '''

rule stringtie_merge:
    input: stringtie/{cond}_{rep,\d+}.gtf
    output: 'stringtie/{cond}.gtf'
    shell:
        '''
        module load stringtie
        stringtie --merge -o {output} -v -T 5 -f 0.5 {input}
        '''

rule gffcompare:
    input: rules.stringtie_merge.output
    output: '{cond}.combined.annotated.gtf'
    params:
        gff = config['gff_path']
    shell:
        """
        ~tbrittoborges/bin/gffcompare/gffcompare \ 
        {input} -r {params.gff} -R -Q -M -o {wildcards.cond}.combined
        """

rule gtf2gff3:
    input:
        {cond}.combined.annotated.gtf
    output:
        {cond}.gff
    shell:
        """
        perl ~tbrittoborges/bin/gtf2gff3.pl {input} > {output}
        """
