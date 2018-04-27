# -*- coding: utf-8
"""
Created on 17:07 27/02/2018
Snakemake file for stringtie.

Install qorts
Install JunctionSeq from R:
    source("http://hartleys.github.io/JunctionSeq/install/JS.install.R");
    JS.install();

.. usage:
    snakemake -s junctionSeq.smk --keep-going --cluster \
    'sbatch --job-name junctionseq_pipeline' --cluster-config  \
     cluster.json --jobs 100

"""
__author__ = "Thiago Britto Borges"
__copyright__ = "Copyright 2018, Dieterichlab"
__email__ = "Thiago.BrittoBorges@uni-heidelberg.de"
__license__ = "MIT"

import os
from itertools import combinations

configfile: "config.yml"
NAMES = config["samples"].keys()
SAMPLES = config["samples"].values()
conditions = sorted(set([x.split('_')[0] for x in NAMES]))
comp_names =  ['{}_{}'.format(*x) for x in
               combinations(conditions, 2)]
mapping = {c: [x for x in NAMES if x.startswith(c)] for c in conditions}

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
    input: rules.create_decoder.output
    output: 'junctionseq/mergedOutput/withNovel.forJunctionSeq.gff.gz'
    params:
        gtf=config['gtf_path'],
        min_count=config['min_count']
    shell:
        """
        module load java qorts
        qorts mergeNovelSplices \
        {params.min_count} \
        {params.stranded} \
        junctionseq/rawCts \
        {input} \
        {params.gtf} \
        junctionseq/mergedOutput/
        """

rule junctioseq_analysis:
    input: rules.merge.output
    threads: 10
    shell:
        """
        module load R
        Rscript junctionSeq.R {threads}
        """


