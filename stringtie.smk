# -*- coding: utf-8
"""
Created on 17:07 27/02/2018 
Snakemake file for stringtie.
for some reason this pipelining is falling at the merge rule (KeyError) if run in
the cluster, but works fine locally
.. usage:
    snakemake -s stringtie.smk --jobs 100

"""
__author__ = "Thiago Britto Borges"
__copyright__ = "Copyright 2018, Dieterichlab"
__email__ = "Thiago.BrittoBorges@uni-heidelberg.de"
__license__ = "MIT"

import os

def create_mapping(config):
    """
    Generate the mapping bettwen samples and replicates
    :param config: snakemake configuration file
    :return dict: mapping bettwen sample and replicates
    """
    names = config["samples"].keys()
    conditions = sorted(set([x.split('_')[0] for x in names]))

    return {c: [x for x in names if x.split('_')[0] == c] for c in conditions}


configfile: "config.yml"
NAMES = list(config["samples"].keys())
mapping = create_mapping(config)


rule all:
    input: expand('stringtie/{names}.gtf', names=NAMES),
           expand('stringtie/{cond}.gff', cond=list(mapping.keys()))

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
    input:
        lambda w: expand('stringtie/{n}.gtf', n=mapping[w.cond])
    output: 'stringtie/{cond}.gtf'
    shell:
        '''
        module load stringtie
        stringtie --merge -o {output} -v -T 5 -f 0.5 {input}
        '''

rule gffcompare:
    input: rules.stringtie_merge.output
    output: 'stringtie/{cond}.combined.{cond}.gtf.tmap'
    params:
        gff = config['gff_path']
    shell:
        '''
        ~tbrittoborges/bin/gffcompare/gffcompare {input} -r {params.gff} -R -Q -M -o {wildcards.cond}.combined
        '''

rule gtf2gff3:
    input: 'stringtie/{cond}.gtf'
    output: 'stringtie/{cond}.gff'
    shell:
        '''
        perl ~tbrittoborges/bin/gtf2gff3.pl {input} > {output}
        '''
