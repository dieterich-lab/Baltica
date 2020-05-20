# -*- coding: utf-8
"""
Created on 12:00 07/01/2020
Snakemake file for quality control of RNA-Sequencing dataset focused on junction quality.

Wang, Liguo, Shengqin Wang, and Wei Li. "RSeQC: quality control of RNA-seq experiments."
Bioinformatics 28.16 (2012): 2184-2185 and FastQC https://github.com/s-andrews/FastQC
"""

__author__ = "Thiago Britto Borges"
__copyright__ = "Copyright 2020, Thiago Britto Borges"
__email__ = "tbrittoborges@uni-heidelberg.de"
__license__ = "MIT"

name = config["samples"].keys()
sample = config["samples"].values()
workdir: config.get("path", ".")

include: "symlink.smk"

rule all:
    input:
         expand("qc/fastqc/{name}_fastqc.zip", name=name),
         "ref.bed12",
         expand("qc/rseqc/{name}.inner_distance_plot.pdf", name=name),
         expand("qc/rseqc/{name}.GC_plot.pdf", name=name),
         expand("qc/rseqc/{name}.DupRate_plot.pdf", name=name),
         expand("qc/rseqc/{name}.splice_junction.pdf", name=name),
         expand("qc/rseqc/{name}.splice_events.pdf", name=name),
         expand("qc/rseqc/{name}.junctionSaturation_plot.pdf", name=name),
         expand("qc/rseqc/{name}.infer_experiment.txt", name=name),
         "qc/multiqc/multiqc_report.html"


# based on https://stackoverflow.com/a/50882104
rule fastqc:
    input: "mappings/{name}.bam"
    output: html="qc/fastqc/{name}_fastqc.html",
          zip="qc/fastqc/{name}_fastqc.zip"
    threads: 10
    params: prefix="qc/fastqc/"
    envmodules: 'fastqc'
    shell: "fastqc -t {threads} {input} -o {params.prefix}"

rule ref_annotation_gtf_to_bed:
    input: config['ref']
    output: "ref.bed12"
    conda: "envs/ucsc.yaml"
    envmodules: "ucsc"
    shell:
        """
        gtfToGenePred {input} temp.genepred
        genePredToBed temp.genepred {output}
        rm temp.genepred
        """

rule rseqc_gene_body_coverage:
    input: bed="ref.bed12"
    output: "qc/rseqc/geneBodyCoverage.curves.pdf"
    envmodules: "rseqc"
    params: prefix="qc/rseqc/"
    shell: "geneBody_coverage.py -r {input.bed} -i mappings/ -o {params.prefix}"

rule rseqc_inner_distance:
    input: bed="ref.bed12",
         bam="mappings/{name}.bam"
    output: "qc/rseqc/{name}.inner_distance_plot.pdf"
    params: prefix="qc/rseqc/{name}"
    envmodules: "rseqc"
    shell: "inner_distance.py -i {input.bam} -o {params.prefix} -r {input.bed}"

rule rseqc_read_gc:
    input: "mappings/{name}.bam"
    output: "qc/rseqc/{name}.GC_plot.pdf"
    params: prefix="qc/rseqc/{name}"
    envmodules: "rseqc"
    shell: "read_GC.py -i {input} -o {params.prefix}"

### fix above output
rule rseqc_read_duplication:
    input: "mappings/{name}.bam"
    output: "qc/rseqc/{name}.DupRate_plot.pdf"
    params: prefix="qc/rseqc/{name}"
    envmodules: "rseqc"
    shell: "read_duplication.py -i {input} -o {params.prefix}"

rule rseqc_junction_annotation:
    input: bed="ref.bed12",
         bam="mappings/{name}.bam"
    output: "qc/rseqc/{name}.splice_junction.pdf",
          "qc/rseqc/{name}.splice_events.pdf"
    params: prefix="qc/rseqc/{name}"
    envmodules: "rseqc"
    shell: "junction_annotation.py -i {input.bam} -r {input.bed} -o {params.prefix}"

rule rseqc_junction_saturation:
    input: bed="ref.bed12",
         bam="mappings/{name}.bam"
    output: "qc/rseqc/{name}.junctionSaturation_plot.pdf"
    params: prefix="qc/rseqc/{name}"
    envmodules: "rseqc"
    shell: "junction_saturation.py -i {input.bam} -r {input.bed} -o {params.prefix}"

rule rseqc_infer_experiment:
    input: bed="ref.bed12",
         bam="mappings/{name}.bam"
    output: "qc/rseqc/{name}.infer_experiment.txt"
    envmodules: "rseqc"
    shell: "infer_experiment.py -r {input.bed} -i {input.bam} > {output}"

rule rseqc_bam_stat:
    input: "mappings/{name}.bam",
    output: "qc/rseqc/{name}.bam_stat.txt"
    shell: "bam_stat.py -i {input} > {output}"

rule read_distribution:
    input: bam="mappings/{name}.bam",
         bed="ref.bed12"
    envmodules: "rseqc"
    output: "qc/rseqc/{name}.read_distribution.txt"
    shell: "read_distribution.py -i {input.bam} -r {input.bed} &> {output}"

rule multiqc:
    input: rules.all.input[:-1]
    envmodules: "multiqc"
    output: "qc/multiqc/multiqc_report.html"
    shell: "multiqc --full -d --dirs-depth 1 qc/ -o qc/multiqc"
