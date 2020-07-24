# -*- coding: utf-8
"""
Created on 20/05/2020
Snakemake workflow for analysis within Baltica.

usage:
    snakemake -s majiq.smk --configfile config.yaml -j 10
"""

try:
    import baltica
    baltica_installed = True
except ImportError:
    baltica_installed = False

workdir: config.get("path", ".")

name = config["samples"].keys()
contrasts = config["contrasts"]

def dir_source(script, ex):
    return script if baltica_installed else srcdir(f"{ex} ../scripts/{script}")

rule all:
    input:
        "majiq/majiq_junctions.csv",
        "leafcutter/leafcutter_junctions.csv",
        "junctionseq/junctionseq_junctions.csv",
        "results/SJ_annotated.csv",
        "results/SJ_annotated_assigned.csv",
        "results/SJ_annotated_assigned_simple.xlsx"
        
rule parse_majiq:
    input:
        expand("majiq/voila/{contrast}_voila.tsv", contrast=config["contrasts"].keys())
    output:
        "majiq/majiq_junctions.csv"
    envmodules:
        "R/3.6.0"
    params:
        cutoff = 0.90,
    shell:
        """
        path=$(which parse_majiq_output.R)
        $path --cutoff {params.cutoff}
        """


rule parse_leafcutter:
    input:
        expand("leafcutter/{contrast}/{contrast}_cluster_significance.txt", contrast=contrasts.keys())
    output:
        "leafcutter/leafcutter_junctions.csv"
    envmodules:
        "R/3.6.0"
    params:
        cutoff = 0.05
    shell:
        """
        path=$(which parse_leafcutter_output.R)
        $path --cutoff {params.cutoff}
        """


rule parse_junctionseq:
    input:
        expand("junctionseq/analysis/{contrast}_sigGenes.results.txt.gz", contrast=contrasts.keys())
    output:
        "junctionseq/junctionseq_junctions.csv"
    envmodules:
        "R/3.6.0"
    params:
        cutoff = 0.05
    shell:
        """
        path=$(which parse_junctionseq_output.R)
        $path --cutoff {params.cutoff}
        """


rule annotate:
    input:
        expand("{method}/{method}_junctions.csv", method=['majiq', 'leafcutter', 'junctionseq']),
        ref="stringtie/merged/merged.combined.gtf"
    envmodules:
        "R/3.6.0"
    output:
        "results/SJ_annotated.csv"
    shell:
        """
        path=$(which annotate_SJ.R)
        $path
        """


rule assign_AS_type:
    input:
        "results/SJ_annotated.csv",
        ref="stringtie/merged/merged.combined.gtf"
    envmodules:
        "R/3.6.0"
    output:
        "results/SJ_annotated_assigned.csv"
    shell:
        """
        path=$(which assign_AS_type.R)
        $path
        """
rule simplify:
    input:
        "results/SJ_annotated_assigned.csv",
    envmodules:
        "R/3.6.0"
    output:
        "results/SJ_annotated_assigned_simple.xlsx"
    shell:
        """
        path=$(which simplify.R)
        $path
        """
