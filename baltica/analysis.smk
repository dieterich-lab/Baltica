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


import sys
import pathlib

exec_path = pathlib.Path(sys.executable).parent

rule all:
    input:
        "results/SJ_annotated.csv",
        "results/SJ_annotated_assigned.csv",
        "results/SJ_annotated_assigned_simple.xlsx"
        
rule parse_majiq:
    input:
        expand("majiq/voila/{contrast}_voila.tsv", contrast=config["contrasts"].keys())
    output:
        "majiq/majiq_junctions.csv"
    envmodules:
        "R/4.0.5_deb10"
    params:
        cutoff = 0.05,
    script: str(exec_path / "parse_majiq_output.R")


rule parse_leafcutter:
    input:
        expand("leafcutter/{contrast}/{contrast}_cluster_significance.txt", contrast=contrasts.keys())
    output:
        "leafcutter/leafcutter_junctions.csv"
    envmodules:
        "R/4.0.5_deb10"
    params:
        cutoff = 0.05
    script: str(exec_path / "parse_leafcutter_output.R")


rule parse_junctionseq:
    input:
        expand("junctionseq/analysis/{contrast}_sigGenes.results.txt.gz", contrast=contrasts.keys())
    output:
        "junctionseq/junctionseq_junctions.csv"
    envmodules:
        "R/4.0.5_deb10"
    params:
        cutoff = 0.05
    script: str(exec_path / "parse_junctionseq_output.R")


rule parse_rmats:
    input:
        expand("rmats/{contrast}/{st}.MATS.JC.txt", 
            contrast=contrasts.keys(),
            st=['A3SS', 'A5SS', 'RI', 'MXE', "SE"])
    output:
        "rmats/rmats_junctions.csv"
    envmodules:
        "R/4.0.5_deb10"
    params:
        cutoff = 0.05
    script: str(exec_path / "parse_rmats_output.R")



rule annotate:
    input:
        expand("{method}/{method}_junctions.csv", method=['majiq', 'leafcutter', 'junctionseq', 'rmats']),
        ref="stringtie/merged/merged.combined.gtf"
    params:
        ref=config.get("ref")
    envmodules:
        "R/4.0.5_deb10"
    output:
        "results/SJ_annotated.csv"
    script: str(exec_path / "annotate_SJ.R")


rule assign_AS_type:
    input:
        "results/SJ_annotated.csv",
        ref="stringtie/merged/merged.combined.gtf"
    envmodules:
        "R/4.0.5_deb10"
    output:
        "results/SJ_annotated_assigned.csv"
    script: str(exec_path / "assign_AS_type.R")


rule simplify:
    input:
        "results/SJ_annotated_assigned.csv",
    envmodules:
        "R/4.0.5_deb10"
    output:
        "results/SJ_annotated_assigned_simple.xlsx"
    script: str(exec_path / "simplify.R")