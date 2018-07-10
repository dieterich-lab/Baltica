# -*- coding: utf-8
"""
Created on 17:07 29/02/2018 2018
Snakemake pipeline to prepare all needed files for pygenomictracks
splice junction visualisation
.. usage:
    sbatch submit_smk.sh pygenomictrack.smk
"""
__author__ = "Thiago Britto Borges"
__copyright__ = "Copyright 2018, Dieterichlab"
__email__ = "Thiago.BrittoBorges@uni-heidelberg.de"
__license__ = "MIT"

configfile: "config.yml"
NAMES = list(config["samples"].keys())

rule all:
    input:
        expand('tracks/{names}.bw', names=NAMES),
        expand('tracks/{names}.arc', names=NAMES),
        '/tracks/genes.bed'

rule bamCoverage:
      # --scaleFactor can be used for lib correction
      # problems with inputfiles when using smk job in the clusterK
    input: 'mappings/{names}.bam'
    output: 'tracks/{names}.bw'
    shell:
        '''
        bamCoverage -b {input} -o {output} --outFileFormat bigwig \
        --numberOfProcessors {threads} --normalizeUsing BPM
        '''

rule juctionCounts:
    input: 'leafcutter/{names}.junc'
    output: 'tracks/{names}.arc'
    params:
        factor = 1
    shell:
        "awk '{{print $1, $2, $2+1, $1, $3, $3+1, $5*a}}' a={params.factor} OFS='\t' {input} > {output}"


rule gtf2bed:
    input: config['gff_path']
    output: 'tracks/genes.bed'
    shell:
        """
        module load bedops
        convert2bed -i gff < {input} > {output}
        """
