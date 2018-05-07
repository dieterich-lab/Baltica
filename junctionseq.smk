# -*- coding: utf-8
"""
Created on 17:07 27/02/2018
Snakemake file for stringtie.
.. install :
Install qorts
Install JunctionSeq from R:
    source("http://hartleys.github.io/JunctionSeq/install/JS.install.R");
    JS.install();

.. usage:
    sbatch submit_smk.sh junctionseq.smk

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
comp_names = ['{}_{}'.format(*x) for x in
               combinations(conditions, 2)]

mapping = {c: [x for x in NAMES if x.split('_')[0] == c] for c in conditions}

rule all:
    input:
        expand('junctionseq/rawCts/{NAMES}/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt.gz',
            NAMES=NAMES),
        'junctionseq/decoder.tab',
        'junctionseq/mergedOutput/withNovel.forJunctionSeq.gff.gz'


for condition, replicate in mapping.items():
    rule:
        input: 'mappings/{NAMES}.bam'
        output: 'junctionseq/rawCts/{NAMES}/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt.gz'
        params:
            gtf=config['gtf_path'],
            output='junctionseq/rawCts/{NAMES}/',
            max_read=config['max_read_length'],
        shell:
            """
            module load java qorts
            qorts QC --stranded {params.max_read} \
            --runFunctions writeKnownSplices,writeNovelSplices,writeSpliceExon \
            {input} {params.gtf} {params.output}
            """

rule create_decoder:
    output: 'junctionseq/decoder.tab'
    run:
        with open(str(output), 'w') as fou:
            fou.write('sample.ID\tgroup.ID\n')
            for k, v in mapping.items():
                for name in v:
                    fou.write('{}\t{}\n'.format(
                        name, k))

rule merge:
    input:
        rules.create_decoder.output,
        expand('junctionseq/rawCts/{NAMES}/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt.gz', NAMES=NAMES)
    output: 'junctionseq/mergedOutput/withNovel.forJunctionSeq.gff.gz'
    params:
        gtf=config['gtf_path'],
        min_count=config['min_count']
    shell:
        """
        module load java qorts
        qorts mergeNovelSplices \
        --minCount 6 \
        --stranded \
        junctionseq/rawCts \
        {input[0]} \
        {params.gtf} \
        junctionseq/mergedOutput/
        """

rule junctioseq_analysis:
    input: rules.merge.output, rules.create_decoder.output
    threads: 10
    shell:
        """
        module load R
        Rscript junctionSeq.R {threads}
        """
