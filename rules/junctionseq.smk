# -*- coding: utf-8
"""
Created on 17:07 27/02/2018
Snakemake workflow for JunctionSeq.

If you use JunctionSeq, please cite:

Hartley, Stephen W., and James C. Mullikin. "Detection and visualization of differential splicing in RNA-Seq data with
 JunctionSeq." Nucleic acids research 44.15 (2016): e127-e127.

..Usage:
    snakemake -s junctionseq.smk --configfile config.yml
"""
__author__ = "Thiago Britto Borges"
__copyright__ = "Copyright 2018, Dieterichlab"
__email__ = "Thiago.BrittoBorges@uni-heidelberg.de"
__license__ = "MIT"

from itertools import groupby

name = config["samples"].keys()
gtf_path = config["ref"]
comp_names = config["contrasts"].keys()

cond, rep = glob_wildcards("mappings/{cond}_{rep}.bam")
d = {k: list(v) for k, v in groupby(
    sorted(zip(cond, rep)), key=lambda x: x[0])}
cond = set(cond)

rule all:
    input:
         expand("junctionseq/rawCts/{name}/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt.gz",
                name=name),
         "junctionseq/decoder.tab",
         "junctionseq/mergedOutput/withNovel.forJunctionSeq.gff.gz",
         "junctionseq/analysis/"

include: "symlink.smk"

rule qc:
    input:
         "mappings/{name}.bam"
    output:
          "junctionseq/rawCts/{name}/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt.gz"
    params:
          gtf=config["ref"],
          output="junctionseq/rawCts/{name}/",
          strandness="--stranded" if config.get("stradness") else "",
          max_read=config["read_len"],
    shell:
         "module load java qorts "
         "qorts QC --stranded --maxReadLength {params.max_read} "
         "--runFunctions writeKnownSplices,writeNovelSplices,writeSpliceExon "
         "{input} {params.gtf} {params.output} "

rule create_decoder:
    output: "junctionseq/decoder.tab"
    run:
        with open(str(output), "w") as fou:
            fou.write("sample.ID\tgroup.ID\n")
            for n in name:
                fou.write("{}\t{}\n".format(n, n.split("_")[0]))

rule merge:
    input:
         rules.create_decoder.output,
         expand("junctionseq/rawCts/{name}/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt.gz",
                name=name)
    output: "junctionseq/mergedOutput/withNovel.forJunctionSeq.gff.gz"
    params:
          gtf=config["ref"],
          min_count=config.get("mincount", 6),
          strandness="--stranded" if config.get("stradness") else ""
    shell:
         "module load java qorts "
         "qorts mergeNovelSplices "
         "--minCount {params.min_count} "
         "{params.stranded} "
         "junctionseq/rawCts "
         "{input[0]} "
         "{params.gtf} "
         " junctionseq/mergedOutput/ "

rule junctioseq_analysis:
    input:
         rules.merge.output,
         rules.create_decoder.output
    output:
          directory("junctionseq/analysis/")
    threads: 10
    shell:
         "module load R "
         "Rscript scripts/junctionSeq.R {threads} "
