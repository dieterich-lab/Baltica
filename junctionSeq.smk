# -*- coding: utf-8
"""
Created on 17:07 27/02/2018
Snakemake file for stringtie.
.. usage:
    snakemake -s junctionSeq.smk --keep-going --cluster 'sbatch ' --cluster-config  cluster.json --jobs 100

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
    input: expand('junctionseq/{NAMES}/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt.gz', NAMES=NAMES)


for condition, replicate in mapping.items():
    rule:
        input:
            'mappings/{NAMES}.bam'
            'majiq/{condition}/{condition}.ini'.format(condition=condition)
        output:
            expand('majiq/{condition}/{replicate}.majiq',
                replicate=replicate, condition=condition)
        params:
            output='majiq/{condition}'.format(condition=condition),
            annotation=config['gtf_path']
        threads: 20
        shell:
            '''
            module load majiq
            majiq build --conf {input} --nproc {threads} \
            --output {params.output} {params.annotation}
            '''

rule qorts:
    input: 'mappings/{NAMES}.bam'
    output: 'junctionseq/{NAMES}/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt.gz',
    params:
        gtf=config['gtf_path'],
        output='junctionseq/{NAMES}',
        max_read=config['max_read_length'],
    shell:
        """
        module load java qorts
        qorts QC --stranded {params.max_read}   \
        --runFunctions writeKnownSplices,writeNovelSplices,writeSpliceExon  \
        {input} {params.gtf} {params.output}
        """

rule merge:
    input: 'junctionseq/{NAMES}/'
    output: '
    params:
        gtf=config['gtf_path']
        min_count=config['min_count']
    shell:
        """
        qorts mergeNovelSplices  --minCount 6 --stranded  ~/RNPS1/rawCts/ ~/RNPS1/decoder.txt \
        {params.gtf} ~/RNPS1/Cts/
        """

