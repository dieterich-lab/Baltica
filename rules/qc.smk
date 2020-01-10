# -*- coding: utf-8
"""
Created on 12:00 07/01/2020
Snakemake file for quality control of RNA-Sequencing dataset focused on junction quality.
"""

__author__ = "Thiago Britto Borges"
__copyright__ = "Copyright 2020, Thiago Britto Borges"
__email__ = "tbrittoborges@uni-heidelberg.de"
__license__ = "MIT"

configfile: "config_differentiation.yaml"
sample = config["samples"].keys()
shell.prefix(config["qc_env"])


rule all:
    input:
        expand("qc/fastqc/{sample}_fastqc.zip", sample=sample),
        "ref.bed12",
        "qc/rseqc/geneBodyCoverage.curves.pdf",
        expand("qc/rseqc/{sample}.inner_distance_plot.pdf", sample=sample),
        expand("qc/rseqc/{sample}.GC_plot.pdf", sample=sample),
        expand("qc/rseqc/{sample}.DupRate_plot.pdf", sample=sample),
        expand("qc/rseqc/{sample}.splice_junction.pdf", sample=sample),
        expand("qc/rseqc/{sample}.splice_events.pdf", sample=sample),
        expand("qc/rseqc/{sample}.junctionSaturation_plot.pdf", sample=sample),
        expand("qc/rseqc/{sample}.infer_experiment.txt", sample=sample),
        expand("qc/rseqc/{sample}.bam_stat.txt", sample=sample),
        expand("qc/rseqc/{sample}.read_distribution.txt", sample=sample),
        "qc/multiqc/multiqc_report.html"


# from https://stackoverflow.com/a/50882104
rule fastqc:
    input:
        "mappings/{sample}.bam"
    output:
        html="qc/fastqc/{sample}_fastqc.html",
        zip="qc/fastqc/{sample}_fastqc.zip"
    threads: 10
    params:
        prefix = "qc/fastqc/"
    shell:
        """
        fastqc -t {threads} {input} -o {params.prefix}
        """


# conda install ucsc-gtftogenepred
# conda install ucsc-genepredtobed
rule ref_annotation_gtf_to_bed:
    input:
        config['ref']
    output:
        "ref.bed12"
    shell:
        """
        wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
        wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToBed
        chmod 755 gtfToGenePred genePredToBed
        ./gtfToGenePred {input} temp.genepred
        ./genePredToBed temp.genepred {output}
        rm temp.genepred gtfToGenePred genePredToBed
        """

rule rseqc_gene_body_coverage:
    input:
        bed = "ref.bed12"
    output:
        "qc/rseqc/geneBodyCoverage.curves.pdf"
    params:
        prefix = "qc/rseqc/"
    shell:
        """
        geneBody_coverage.py -r {input.bed} -i mappings/ -o {params.prefix}
        """

rule rseqc_inner_distance:
    input:
        bed = "ref.bed12",
        bam = "mappings/{sample}.bam"
    output: "qc/rseqc/{sample}.inner_distance_plot.pdf"
    params:
        prefix = "qc/rseqc/{sample}"
    shell:
        """
        inner_distance.py -i {input.bam} -o {params.prefix} -r {input.bed}
        """

rule rseqc_read_gc:
    input: "mappings/{sample}.bam"
    output: "qc/rseqc/{sample}.GC_plot.pdf"
    params:
        prefix = "qc/rseqc/{sample}"
    shell:
        """
        read_GC.py -i {input} -o {params.prefix}
        """

### fix above output
rule rseqc_read_duplication:
    input: "mappings/{sample}.bam"
    output: "qc/rseqc/{sample}.DupRate_plot.pdf"
    params:
        prefix = "qc/rseqc/{sample}"
    shell:
        """
        read_duplication.py -i {input} -o {params.prefix}
        """


rule rseqc_junction_annotation:
    input:
        bed = "ref.bed12",
        bam = "mappings/{sample}.bam"
    output:
        "qc/rseqc/{sample}.splice_junction.pdf",
        "qc/rseqc/{sample}.splice_events.pdf"
    params:
        prefix = "qc/rseqc/{sample}"
    shell:
        """
        junction_annotation.py -i {input.bam} -r {input.bed} -o {params.prefix}
        """

rule rseqc_junction_saturation:
    input:
        bed = "ref.bed12",
        bam = "mappings/{sample}.bam"
    output: "qc/rseqc/{sample}.junctionSaturation_plot.pdf"
    params:
        prefix = "qc/rseqc/{sample}"
    shell:
        """
        junction_saturation.py -i {input.bam} -r {input.bed} -o {params.prefix}
        """

rule rseqc.infer_experiment:
    input:
        bed = "ref.bed12",
        bam = "mappings/{sample}.bam"
    output: "qc/rseqc/{sample}.infer_experiment.txt"
    shell:
        """
        infer_experiment.py -r {input.bed} -i {input.bam} > {output}
        """

rule rseqc_bam_stat:
    input: "mappings/{sample}.bam",
    output: "qc/rseqc/{sample}.bam_stat.txt"
    shell:
        """
        bam_stat.py -i {input} > {output}
        """


rule read_distribution:
    input:
        bam = "mappings/{sample}.bam",
        bed = "ref.bed12"
    output:
        "qc/rseqc/{sample}.read_distribution.txt"
    shell:
        """
        read_distribution.py -i {input.bam} -r {input.bed} &> {output}
        """

rule multiqc:
    input:
        rules.all.input[:-1]
    output:
        "qc/multiqc/multiqc_report.html"
    shell:
        """
        multiqc --fullnames -d --dirs-depth 1 qc/ -o qc/multiqc
        """