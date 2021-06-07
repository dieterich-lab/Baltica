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

from itertools import groupby, chain
try:
    import baltica
    baltica_installed = True
except ImportError:
    baltica_installed = False

workdir: config.get("path", ".")
name = config["samples"].keys()
sample = config["samples"].values()
gtf_path = config["ref"]
comp_names = config["contrasts"].keys()
cond, rep = glob_wildcards("mappings/{cond}_{rep}.bam")
cond = set(cond)
container: "docker://tbrittoborges/junctionseq:1.16.0"

strandness = {
    'forward': '--stranded_fr_secondstrand',
    'reverse': '--stranded'
}
def dir_source(script, ex):
    return script if baltica_installed else srcdir(f"{ex} ../scripts/{script}")


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

localrules: cat_decoder, create_decoder
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
    envmodules: "java qorts "
    shadow: "shallow"
    shell:
         "qorts QC {params.strandness} {params.is_paired_end} --maxReadLength {params.max_read} "
         "--runFunctions writeKnownSplices,writeNovelSplices,writeSpliceExon "
         "{input} {params.gtf} {params.output} "


rule create_decoder:
    output:
        "junctionseq/{comparison}_decoder.tab"
    run:
        ref, alt = wildcards[0].split('-vs-')
        with open(str(output), "w") as fou:
            fou.write("sample.ID\tgroup.ID\n")
            for n in name:
                if n.startswith(ref) or n.startswith(alt):
                    fou.write("{0}\t{1}\n".format(n, n.split('_')[0]))

rule cat_decoder:
    input: decoder=expand("junctionseq/{comparison}_decoder.tab", comparison=comp_names )
    output: "junctionseq/decoder.tab"
    shadow: "shallow"
    shell:
        "awk 'FNR>1 || NR==1 ' {input.decoder} | awk '!x[$0]++' > {output} "

rule merge:
    input:
         decoder='junctionseq/decoder.tab',
         counts=expand("junctionseq/rawCts/{name}/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt.gz",
                name=name)
    output: "junctionseq/mergedOutput/withNovel.forJunctionSeq.gff.gz"
    params:
          gtf=config["ref"],
          min_count=config.get("mincount", 6),
          strandness= strandness.get(config.get("strandness"),  "")
    envmodules: "java qorts"
    shadow: "shallow"
    shell:
         "qorts mergeNovelSplices "
         "--minCount {params.min_count} "
         "{params.strandness} "
         "junctionseq/rawCts "
         "{input.decoder} "
         "{params.gtf} "
         "junctionseq/mergedOutput/ "


rule junctioseq_analysis:
    input:
        reads= "junctionseq/mergedOutput/withNovel.forJunctionSeq.gff.gz",
        decoder= "junctionseq/{comparison}_decoder.tab"
    output:
        "junctionseq/analysis/{comparison}_sigGenes.results.txt.gz"
    threads: 10
    params: script = dir_source("junctionSeq.R", "Rscript")
    envmodules: "R/3.6.3_deb10 junctionseq/1.16.0_deb10"
    shadow: "shallow"
    shell:
        "{params.script} {input.decoder} {output} {threads}"
