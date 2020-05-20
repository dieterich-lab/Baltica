# -*- coding: utf-8
"""
Created on 10:07 15/07/2019

"""
__author__ = "Thiago Britto Borges"
__copyright__ = "Copyright 2019, Dieterichlab"
__email__ = "Thiago.BrittoBorges@uni-heidelberg.de"
__license__ = "MIT"


name, = glob_wildcards("mappings/{name}.bam")

configfile: "config.yml"
config["whippet"] = "~/miniconda3/envs/whippet/share/julia/site/v0.6/Whippet/src/"

def rename_bam_to_fastq(x):
  # reads star log for input files
  f = x.replace("Aligned.noS.bam", "Log.out")
  # f.replace(f.name, "Log.out"
  with open(f) as fin:
    for line in fin:^
      if line.startswith("##### Command Line:"):
        line = next(fin)
        break

  for args in line.split("--"):
    if args.startswith("readFilesIn"):
      break
  return [x.split("mapping")[0] + arg for arg in args.split()[1:3]]

rule all:
    input:
         pass

rule merge_bam:
    input: expand("mappings/{name}.bam", name=name)
    output: tmpdir("merged.bam")
    shell: "samtools merge {output} {input}"

rule sort_bam:
    input: rules.merge_bam.output
    output: tmpdir("merged.sort.bam")
    shell: "samtools sort {input} {output}"

rule rmdup_bam:
    input: rules.sort_bam.output
    output: "whippet/merged.sort.bam"
    shell: "samtools rmdup -S {input} {output}"

rule index_bam:
    input: rules.rmdup_bam.output
    output: "whippet/merged.sort.bam.bai"
    shell: "samtools index filename.sort.rmdup.bam"

rule whippet_index:
    input:
         bam = "whippet/merged.sort.bam",
         bai = "whippet/merged.sort.bam.bai"
    output: "whippet/index.jls"
    params:
        fa = config["fa"],
        gtf = config["gtf"]
    shell:
        "{config.whippet}/whippet-index.jl --fasta {params.fa} --gtf {params.gtf} --bam {input} -o {output}"

rule whippet_quat:
    input:
        idx = rules.whippet_index.output,
        R1 = "",
    #          r1 = lambda wildcards: FILES[wildcards.sample][0],
    # r2 = lambda wildcards: FILES[wildcards.sample][1],
    #     #
        R2 = ""
    output: "whippet/output_{wc.name}"
    shell:
        "{config.whippet}/whippet-quant.jl {input.R1} {input.R2} -o {output} -x {input.idx}"