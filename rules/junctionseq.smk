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
__copyright__ = "Copyright 2020, Dieterichlab"
__email__ = "Thiago.BrittoBorges@uni-heidelberg.de"
__license__ = "MIT"

from itertools import groupby

workdir: config.get("path", ".")
name = config["samples"].keys()
sample = config["samples"].values()
gtf_path = config["ref"]
comp_names = config["contrasts"].keys()
cond, rep = glob_wildcards("mappings/{cond}_{rep}.bam")
d = {k: list(v) for k, v in groupby(
    sorted(zip(cond, rep)), key=lambda x: x[0])}
cond = set(cond)
strandness = {
    'forward': '--stranded --fr_secondStrand',
    'reverse': '--stranded'
}
rule all:
    input:
        "logs/",
        expand(
            "junctionseq/rawCts/{name}/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt.gz",
                name=name),
        expand("junctionseq/{comparison}_decoder.tab",
            comparison=comp_names),
        "junctionseq/mergedOutput/withNovel.forJunctionSeq.gff.gz",
        expand("junctionseq/analysis/{comparison}_sigGenes.results.txt.gz",
            comparison=comp_names)


include: "symlink.smk"

rule qc:
    input:
        "mappings/{name}.bam"
    output:
        "junctionseq/rawCts/{name}/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt.gz"
    params:
        gtf=config["ref"],
        output="junctionseq/rawCts/{name}/",
        strandness= strandness.get(config.get("strandness"),  ""),
        max_read=config["read_len"],
	is_paired_end='--singleEnded' if config.get("is_single_end") == True  else ''
    envmodules:
        "java qorts "
    shell:
         "qorts QC {params.strandness} {params.is_paired_end} --maxReadLength {params.max_read} "
         "--runFunctions writeKnownSplices,writeNovelSplices,writeSpliceExon "
         "{input} {params.gtf} {params.output} "


rule create_decoder:
    output:
        "junctionseq/{comparison}_decoder.tab"
    run:
        ref, alt = wildcards[0].split('-vs-')
        with open(str( output), "w") as fou:
            fou.write("sample.ID\tgroup.ID\n")
            for x in [*d[ref], *d[alt]]:
                fou.write("{}_{}\t{}\n".format(x[0], x[1], x[0]))


rule merge:
    input:
         decoder=expand("junctionseq/{comparison}_decoder.tab", comparison=comp_names ),
         counts=expand("junctionseq/rawCts/{name}/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt.gz",
                name=name)
    output: "junctionseq/mergedOutput/withNovel.forJunctionSeq.gff.gz"
    params:
          gtf=config["ref"],
          min_count=config.get("mincount", 6),
          strandness= strandness.get(config.get("strandness"),  "")
    envmodules:
        "java qorts"
    shell:
         "cat {input.decoder} > junctionseq/decoder.tab; "
         "qorts mergeNovelSplices "
         "--minCount {params.min_count} "
         "{params.strandness} "
         "junctionseq/rawCts "
         "junctionseq/decoder.tab "
         "{params.gtf} "
         "junctionseq/mergedOutput/ "


rule junctioseq_analysis:
    input:
        reads= "junctionseq/mergedOutput/withNovel.forJunctionSeq.gff.gz",
        decoder= "junctionseq/{comparison}_decoder.tab"
    output:
        "junctionseq/analysis/{comparison}_sigGenes.results.txt.gz"
    threads: 10
    envmodules:
	"R/3.5.1 junctionseq"
    script:
        "../scripts/junctionSeq.R"
