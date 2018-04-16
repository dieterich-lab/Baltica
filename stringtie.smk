# -*- coding: utf-8
"""
Created on 17:07 27/02/2018 
Snakemake file for stringtie.
.. usage:
    snakemake -s stringtie.Snakemake --cluster "sbatch --mem=24000 -n stringtie" --jobs 100
    
"""
__author__ = "Thiago Britto Borges"
__copyright__ = "Copyright 2018, Dieterichlab"
__email__ = "Thiago.BrittoBorges@uni-heidelberg.de"
__license__ = "MIT"

import os
# input file with target files 1 per line
configfile: "config.yml"
NAMES = config["samples"].keys()
SAMPLES = config["samples"].values()

rule all:
    input: expand('StringTieNoRef/{names}.gtf', names=NAMES)


rule symlink:
    input:
        bam = expand('{samples}', samples=SAMPLES),
        bai = expand('{samples}.bai', samples=SAMPLES)
    output:
        bam = expand('mappings/{names}.bam', names=NAMES),
        bai = expand('mappings/{names}.bam.bai', names=NAMES)
    run:
        for bam_in, bai_in, bam_out, bai_out in zip(
            input.bam, input.bai, output.bam, output.bai):
            os.symlink(bam_in, bam_out)
            os.symlink(bai_in, bai_out)

rule stringtie:
    input: 'mappings/{names}.bam'
    output: 'StringTieNoRef/{names}.gtf'
    threads: 8 
    params:
        min_iso = config.get('min_iso', ""),
        min_anchor_len = config.get('min_anchor_len', ""),
        min_junction_cov = config.get('min_junction_cov', ""),
        rf = config.get('rf', "")
    shell:
        '''
        module load stringtie
        stringtie {input}  -v -p {threads} -o {output} \
        {params.rf}    {params.min_iso} {params.min_anchor_len} \
        {params.min_junction_cov}
        '''