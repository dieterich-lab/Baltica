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
        "junctionseq/junctionseq_junctionse.csv",
    # "{method}/{method}_junctions_annotated.csv"


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
    script: "../analysis/parse_leafcutter_output.R"

rule parse_junctionseq:
    input:
        expand("junctionseq/analysis/{contrast}_sigGenes.results.txt.gz", contrast=contrasts.keys())
    output:
        "junctionseq/junctionseq_junctionse.csv"
    envmodules:
        "R/3.6.0"
    params:
        cutoff = 0.05
    script: "../analysis/parse_junctionseq_output.R"
#
# rule annotate:
#     input:
#         table = "{method}/{method}_junctions.csv",
#         annotation = "stringtie/merged/merged.combined.gtf"
#     output:
#         "{method}/{method}_junctions_annotated.csv"
#     shell:
#         "Rscript analysis/annotate_SJ.R -i {input.table} -a {input.annotation} -o {output}"Rscript analysis/annotate_SJ.R -i majiq/majiq_junctions.csv -a denovo_tx/merged/merged.combined.gtf -o majiq/majiq_junctions_annotated.csv
#
#
# # rule write_simplified_result:
# #     input:
# #     output:
# #     shell: "Rscript simplify.R"
