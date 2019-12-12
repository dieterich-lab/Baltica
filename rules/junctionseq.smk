# -*- coding: utf-8
"""
Created on 17:07 27/02/2018

.. usage:
    sbatch submit_smk.sh 'junctionseq.smk --allow-ambiguity junctioseq_analysis'
"""
__author__ = "Thiago Britto Borges"
__copyright__ = "Copyright 2018, Dieterichlab"
__email__ = "Thiago.BrittoBorges@uni-heidelberg.de"
__license__ = "MIT"

from itertools import groupby

configfile: "config.yml"
name = config["samples"].keys()
gtf_path = config["ref"]
comp_names = config['contrasts'].keys()

cond, rep = glob_wildcards("mappings/{cond}_{rep}.bam")
d = {k: list(v) for k, v in groupby(
    sorted(zip(cond, rep)), key=lambda x: x[0])}
cond = set(cond)

rule all:
    input:
        expand('junctionseq/rawCts/{name}/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt.gz',
            name = name),
        'junctionseq/decoder.tab',
        'junctionseq/mergedOutput/withNovel.forJunctionSeq.gff.gz',
        'junctionseq/analysis/'


rule qc:
    input:
        "mappings/{name}.bam"
    output:
        'junctionseq/rawCts/{name}/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt.gz'
    params:
        gtf = config['ref'],
        output = 'junctionseq/rawCts/{name}/',
        max_read = config['read_len'],
    shell:
        """
        module load java qorts
        qorts QC --stranded --maxReadLength {params.max_read} \
        --runFunctions writeKnownSplices,writeNovelSplices,writeSpliceExon \
        {input} {params.gtf} {params.output}
        """


rule create_decoder:
    output: 'junctionseq/decoder.tab'
    run:
        with open(str(output), 'w') as fou:
            fou.write('sample.ID\tgroup.ID\n')
            for n in name:
                fou.write('{}\t{}\n'.format(n, n.split('_')[0]))


rule merge:
    input:
        rules.create_decoder.output,
        expand('junctionseq/rawCts/{name}/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt.gz',
            name = name)
    output: 'junctionseq/mergedOutput/withNovel.forJunctionSeq.gff.gz'
    params:
        gtf = config['ref']
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
    input:
        rules.merge.output,
        rules.create_decoder.output
    output:
        directory('junctionseq/analysis/')
    threads: 10
    shell:
        """
        module load R
        Rscript scripts/junctionSeq.R {threads}
        """
