# -*- coding: utf-8
"""
Created on 17:07 29/02/2018 2018
Snakemake workflow for LeafCutter.

If you use LeafCutter, please cite:

Li, Yang I., et al. "Annotation-free quantification of RNA splicing using LeafCutter." Nature genetics 50.1 (2018): 151-158.

..Usage:
    snakemake -s leafcutter.smk --configfile config.yml
"""
__author__ = "Thiago Britto Borges"
__copyright__ = "Copyright 2020, Dieterichlab"
__email__ = "Thiago.BrittoBorges@uni-heidelberg.de"
__license__ = "MIT"

from pathlib import Path

try:
    import baltica
    baltica_installed = True
except ImportError:
    baltica_installed = False

def dir_source(script, ex):
    return script if baltica_installed else srcdir(f"{ex} ../scripts/{script}")


def basename(path, suffix=None):
    if suffix:
        return str(Path(path).with_suffix(suffix).name)
    return str(Path(path).name)


workdir: config.get("path", ".")
name = config["samples"].keys()
sample = config["samples"].values()
gtf_path = config["ref"]
conditions = [x.split("_")[0] for x in name]
sample_path = config["sample_path"]
comp_names = config["contrasts"].keys()

localrules: all, leafcutter_concatenate, symlink

include: "symlink.smk"

if "leafcutter_env_prefix" in config:
    shell.prefix(config["leafcutter_env_prefix"])

rule all:
    input:
        'logs/',
        expand("mappings/{name}.bam", name=name),
        expand("leafcutter/{comp_names}/{comp_names}_cluster_significance.txt", comp_names=comp_names),
        expand("leafcutter/{name}.junc", name=name)

# step 1.1
rule leafcutter_bam2junc:
    input: "mappings/{name}.bam"
    output:
          bed="leafcutter/{name}.bed",
          junc="leafcutter/{name}.junc"
    params:
          filter_cs_path=dir_source("filter_cs.py", "python"),
          sam2bed_path=dir_source("sam2bed.pl", "perl"),
          bed2junc_path=dir_source("bed2junc.pl", "perl"),
          use_strand="--use-RNA-strand" if config.get("strandness") else ""
    envmodules: "samtools"
    conda: "../envs/leafcutter.yml"
    shadow: "shallow"
    shell:
         """
         samtools view {input} \
         | {params.filter_cs_path} \
         | {params.sam2bed_path} {params.use_strand} - {output.bed}
          {params.bed2junc_path} {output.bed} {output.junc}
         """

rule leafcutter_concatenate:
    input: expand(rules.leafcutter_bam2junc.output, name=name)
    output:
          junc="leafcutter/{comp_names}/juncfiles.txt",
          test="leafcutter/{comp_names}/diff_introns.txt"
    run:
        comp = output.junc.split("/")[1]
        cond_a, cond_b = comp.split("-vs-")
        with open(output.junc, 'w') as file_out:
            work_path = Path(config.get("path", "."))
            for n in name:
                if n.startswith(cond_a) or n.startswith(cond_b):
                    file_out.write(str(work_path / "leafcutter/{}.junc\n".format(n)))

        with open(output.test, 'w') as file_out:
            for n in name:
                if n.startswith(cond_a):
                    file_out.write("{} {}\n".format(n, cond_a))

                elif n.startswith(cond_b):
                    file_out.write("{} {}\n".format(n, cond_b))
# step 2
# m= min numb reads per cluster, l = intron length
rule leafcutter_intron_clustering:
    input: rules.leafcutter_concatenate.output.junc
    params:
          m=config.get('leafcutter_min_cluster_reads', 30),
          l=config.get('leafcutter_max_intron_length', 500000),
          prefix="leafcutter/{comp_names}/{comp_names}",
          n="{comp_names}",
          strand="--strand True" if config.get("strandness") is not None else "",
          script_path=dir_source("leafcutter_cluster.py", "python")
    output: "leafcutter/{comp_names}/{comp_names}_perind_numers.counts.gz"
    conda: "../envs/leafcutter.yml"
    shadow: "shallow"
    shell:
         """
         {params.script_path} \
         -j {input} -m {params.m} -o {params.prefix} -l {params.l} {params.strand}
         rm *{params.n}.sorted.gz
         """

# step 3.1
rule leafcutter_gtf_to_exon:
    input: gtf_path
    output:
          a="leafcutter/" + basename(gtf_path, suffix=".gz"),
          b="leafcutter/exons.gtf.gz"
    params: gtf_to_exon=dir_source("gtf_to_exons.R", "Rscript")
    conda: "../envs/leafcutter.yml"
    envmodules:
        "R/3.6.0"
    shadow: "shallow"
    shell:
         """
         gzip -c {input} > {output.a}
         {params.gtf_to_exon} {output.a} {output.b}
         """

# step 3.2
rule leafcutter_differential_splicing:
    input:
         a=rules.leafcutter_gtf_to_exon.output.b,
         b="leafcutter/{comp_names}/{comp_names}_perind_numers.counts.gz",
         c="leafcutter/{comp_names}/diff_introns.txt"
    output: "leafcutter/{comp_names}/{comp_names}_cluster_significance.txt"
    params:
          min_samples_per_group=config.get("leafcutter_min_samples_per_group", 3),
          min_samples_per_intron=config.get("leafcutter_min_samples_per_intron", 5),
          min_coverage=config.get("leafcutter_min_coverage", 20),
          prefix="leafcutter/{comp_names}/{comp_names}",
          leafcutter_ds_path=dir_source("leafcutter_ds_pair.R", "Rscript")
    threads: 10
    conda: "../envs/leafcutter.yml"
    envmodules:
        "R/3.6.0 leafcutter"
    shadow: "shallow"
    shell:
         """
         {params.leafcutter_ds_path} --exon_file={input.a} \
         --min_coverage {params.min_coverage}  \
         {input.b} {input.c} --num_threads {threads} --output_prefix={params.prefix} \
         -i {params.min_samples_per_intron} -g {params.min_samples_per_group}
         """
