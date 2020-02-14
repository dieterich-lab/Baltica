# -*- coding: utf-8
"""
Created on 17:07 27/07/2018
Snakemake file for denovo transcriptomics followed by quantification
with Scallop and Salmon pipeline.
"""
__author__ = "Thiago Britto Borges"
__copyright__ = "Copyright 2019, Dieterichlab"
__email__ = "Thiago.BrittoBorges@uni-heidelberg.de"
__license__ = "MIT"

from itertools import groupby
from os.path import join


def rename_bam_to_fastq(x):
  # reads star log for input files
  f = x.replace("Aligned.noS.bam", "Log.out")
  # f.replace(f.name, 'Log.out'
  with open(f) as fin:
    for line in fin:
      if line.startswith("##### Command Line:"):
        line = next(fin)
        break

  for args in line.split("--"):
    if args.startswith("readFilesIn"):
      break
  return [x.split("mapping")[0] + arg for arg in args.split()[1:3]]

strand = {
  'reverse': '--fr',
  'foward': '--rf'}

cond, rep = glob_wildcards("mappings/{cond}_{rep}.bam")

name = config["samples"].keys()
raw_name = config["samples"].values()
sample_path = config["sample_path"]
FILES = {k: rename_bam_to_fastq(join(sample_path, v))
         for k, v in zip(name, raw_name)}


d = {k: list(v) for k, v in groupby(
    sorted(zip(cond, rep)), key=lambda x: x[0])}

cond = set(cond)

rule all:
  input:
    expand("mappings/{name}.bam", name=name),
    expand("denovoTx/merged_bam/{group}.bam", group=cond),
    expand("denovoTx/denovoTx/{group}.gtf", group=cond),
    expand("salmon/{sample}/quant.sf", sample=name),
    "denovoTx/merged/merged.fa",
    "denovoTx/salmon/salmon_index/",


rule merge_bam:
  input:
    lambda wc:
      ["mappings/{}_{}.bam".format(*x) for x in d[wc.group]]
  output:
    bam = "denovoTx/merged_bam/{group}.bam",
    bai = "denovoTx/merged_bam/{group}.bam.bai"
  conda:
       srcdir("../envs/scallop.yml")
  threads:
    10
  wildcard_constraints:
    group = "|".join(cond)
  envmodules: "samtools"
  shell:
    "samtools merge {output.bam} {input} --threads {threads};"
    "samtools index {output.bam} {output.bai} "


rule denovo_transcriptomics:
  input:
    "denovoTx/merged_bam/{group}.bam"
  output:
    "denovoTx/denovoTx/{group}.gtf"
  conda:
       srcdir("../envs/scallop.yml")
  params:
    strandness = strand[config.get("strandness", "")],
    min_junct_coverage = 3,
    min_isoform_proportion = .1
  wildcard_constraints: group = "|".join(cond)
  log: "logs/denovoTx_{group}.log"
  envmodules: "stringtie"
  shell:
    "stringtie {input} -o {output} "
    "-p {threads} "
    " {params.strandness} "
    "-j {params.min_junct_coverage} "
    "-f {params.min_isoform_proportion} "
    "2> {log} "


# Merge all denovo annotations
# to use a single index for Salmon
rule merge_gtf:
  input: expand("denovoTx/denovoTx/{cond}.gtf", cond=cond)
  output: "denovoTx/merged/merged.combined.gtf"
  conda: srcdir("../envs/scallop.yml")
  log: "logs/gffcompare.log"
  params: gtf = config["ref"]
  envmodules: "stringtie"
  shell:
    "gffcompare {input} -r {params.gtf} "
    "-R -V -o denovoTx/merged/merged 2>{log}"


rule extract_sequences:
  input: rules.merge_gtf.output
  output: "denovoTx/merged/merged.fa"
  conda: srcdir("../envs/scallop.yml")
  envmodules: "salmon"
  log: "logs/gffread.log"
  params: fasta = config["ref_fa"]
  shell: "gffread {input} -g {params.fasta} -w {output} 2>{log}"


rule salmon_index:
  input: rules.extract_sequences.output
  output: directory("denovoTx/salmon/salmon_index/")
  conda: srcdir("../envs/scallop.yml")
  envmodules: "salmon"
  threads: 10
  log: "logs/salmon_index.log"
  shell: "salmon index -t {input} -i {output} -p {threads} 2> {log}"


rule salmon_quant:
  input:
    r1 = lambda wildcards: FILES[wildcards.sample][0],
    r2 = lambda wildcards: FILES[wildcards.sample][1],
    index = "denovoTx/salmon/salmon_index/"
  output: "salmon/{sample}/quant.sf"
  conda: srcdir("../envs/scallop.yml")
  envmodules: "salmon"
  threads: 10
  log: "logs/{sample}_salmons_quant.log"
  shell:
    "salmon quant -p {threads} -i {input.index} "
    "--libType A "
    "-1 <(gunzip -c {input.r1}) "
    "-2 <(gunzip -c {input.r2}) "
    "-o salmon/{wildcards.sample}/ &> {log}"
