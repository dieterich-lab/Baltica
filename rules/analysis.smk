# -*- coding: utf-8
"""
Created on 20/05/2020
Snakemake workflow for analysis within Baltica.

usage:
    snakemake -s majiq.smk --configfile config.yaml -j 10
"""
workdir: config.get("path", ".")

name = config["samples"].keys()
contrasts = config["contrasts"]


rule all:
    input:
        "majiq/majiq_junctions.csv",
        "leafcutter/leafcutter_junctions.csv",
        "junctionseq/junctionseq_junctions.csv",
        "results/SJ_annotated.csv",
        "results/SJ_annotated_assigned.csv"

rule parse_majiq:
    input:
        expand("majiq/voila/{contrast}_voila.tsv", contrast=config["contrasts"].keys())
    output:
        "majiq/majiq_junctions.csv"
    envmodules:
        "R/3.6.0"
    params:
        cutoff = 0.90
    script:
        '../analysis/parse_majiq_output.R'


rule parse_leafcutter:
    input:
        expand("leafcutter/{contrast}/{contrast}_cluster_significance.txt", contrast=contrasts.keys())
    output:
        "leafcutter/leafcutter_junctions.csv"
    envmodules:
        "R/3.6.0"
    params:
        cutoff = 0.05
    script:
        "../analysis/parse_leafcutter_output.R"


rule parse_junctionseq:
    input:
        expand("junctionseq/analysis/{contrast}_sigGenes.results.txt.gz", contrast=contrasts.keys())
    output:
        "junctionseq/junctionseq_junctions.csv"
    envmodules:
        "R/3.6.0"
    params:
        cutoff = 0.05
    script:
        "../analysis/parse_junctionseq_output.R"


rule annotate:
    input:
        expand("{method}/{method}_junctions.csv", method=['majiq', 'leafcutter', 'junctionseq']),
        "stringtie/merged/merged.combined.gtf"
    envmodules:
        "R/3.6.0"
    output:
        "results/SJ_annotated.csv"
    script:
        "../analysis/annotate_SJ.R"


rule assign_AS_type:
    input:
        "results/SJ_annotated.csv",
        "stringtie/merged/merged.combined.gtf"
    envmodules:
        "R/3.6.0"
    output:
        "results/SJ_annotated_assigned.csv"
    script:
        "../analysis/assign_AS_type.R"
