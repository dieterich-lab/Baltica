# -*- coding: utf-8
"""
Created on 17:07 27/07/2018
Snakemake de novo transcriptomics and transcript abundance estimation
- Find the the fastq files used by DCC read alignments workflow
- compute de novo tx with StringTie [doi:10.1038/nprot.2016.095]
- use de novo annotation to compute transcript abundance with salmon [https://doi.org/10.1038/nmeth.4197]

If you use this workflow, please cite

Patro, R., et al. "Salmon provides fast and bias-aware quantification of transcript expression. Nat Meth. 2017; 14 (4): 417–9."
Pertea, Mihaela, et al. "Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown." Nature protocols 11.9 (2016): 1650.
"""
__author__ = "Thiago Britto Borges"
__copyright__ = "Copyright 2019, Dieterichlab"
__email__ = "Thiago.BrittoBorges@uni-heidelberg.de"
__license__ = "MIT"

from itertools import groupby
import re

def extract_samples_replicates(samples, _pattern=re.compile('^(.+)_(.+)$')):
    """
    Extract pairs of condition and replicate name from sample files

    :param str _pattern: pattern to . Default uses {condition}_{replicate} template
    :param list samples:
    :return:
    :rtype: list
    """
    return list(zip(*[re.match(_pattern, x).groups() for x in samples]))


strand = {
  'reverse': '--fr',
  'foward': '--rf'}

workdir: config.get("path", ".")
cond, rep = extract_samples_replicates(config["samples"].keys())
name = config["samples"].keys()
raw_name = config["samples"].values()
sample_path = config["sample_path"]

d = {k: list(v) for k, v in groupby(
  sorted(zip(cond, rep)), key=lambda x: x[0])}

cond = set(cond)

include: "symlink.smk"

rule all:
  input:
       expand("mappings/{name}.bam", name=name),
       expand("denovo_tx/merged_bam/{group}.bam", group=cond),
       expand("denovo_tx/denovo_tx/{group}.gtf", group=cond),
       "denovo_tx/merged/merged.combined.gtf"


rule merge_bam:
  input: lambda wc: ["mappings/{}_{}.bam".format(*x) for x in d[wc.group]]
  output: bam="denovo_tx/merged_bam/{group}.bam",
        bai="denovo_tx/merged_bam/{group}.bam.bai"
  conda: srcdir("../envs/scallop.yml")
  threads: 10
  wildcard_constraints: group="|".join(cond)
  envmodules: "samtools"
  shell: "samtools merge {output.bam} {input} --threads {threads};" \
       "samtools index {output.bam} {output.bai} "


rule denovo_transcriptomics:
  input: "denovo_tx/merged_bam/{group}.bam"
  output: "denovo_tx/denovo_tx/{group}.gtf"
  conda: srcdir("../envs/scallop.yml")
  params: strandness=strand.get(config.get("strandness", ""), ""),
        min_junct_coverage=3,
        min_isoform_proportion=.1
  wildcard_constraints: group="|".join(cond)
  log: "logs/denovo_tx_{group}.log"
  envmodules: "stringtie"
  shell: "stringtie {input} -o {output} " \
       "-p {threads} " \
       " {params.strandness} " \
       "-j {params.min_junct_coverage} " \
       "-f {params.min_isoform_proportion} " \
       "2> {log} "


rule merge_gtf:
  input: expand("denovo_tx/denovo_tx/{cond}.gtf", cond=cond)
  output: "denovo_tx/merged/merged.combined.gtf"
  conda: srcdir("../envs/scallop.yml")
  log: "logs/gffcompare.log"
  params: gtf=config["ref"]
  envmodules: "rnaseqtools"
  shell: "gffcompare {input} -r {params.gtf} " \
      "-R -V -o denovo_tx/merged/merged 2>{log}"
