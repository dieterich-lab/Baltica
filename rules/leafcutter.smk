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
__copyright__ = "Copyright 2018, Dieterichlab"
__email__ = "Thiago.BrittoBorges@uni-heidelberg.de"
__license__ = "MIT"

from pathlib import Path

def basename(path, suffix=None):
    if suffix:
        return str(Path(path).with_suffix(suffix).name)
    return str(Path(path).name)

NAMES = config["samples"].keys()
SAMPLES = config["samples"].values()
gtf_path = config["ref"]
conditions = [x.split("_")[0] for x in NAMES]
sample_path = config["sample_path"]
comp_names = config["contrasts"].keys()

localrules: all, concatenate, symlink

rule all:
    input:
         expand("leafcutter/{comp_names}/{comp_names}_cluster_significance.txt",
                comp_names=comp_names),
         expand("leafcutter/{NAMES}.junc", NAMES=NAMES)

include: "symlink.smk"

# step 1.1
rule bam2junc:
    input: "mappings/{NAMES}.bam"
    output:
          bed="leafcutter/{NAMES}.bed",
          junc="leafcutter/{NAMES}.junc"
    params:
          filter_cs_path=srcdir("scripts/filter_cs.py"),
          sam2bed_path=srcdir("scripts/sam2bed.pl"),
          bed2junc_path=srcdir("scripts/bed2junc.pl"),
          use_strand="--use-RNA-strand" if config.get("strandness") else ""
    shell:
         """
         samtools view {input} \
         | python2 {params.filter_cs_path} \
         | perl {params.sam2bed_path} {params.use_strand} - {output.bed}
         perl {params.bed2junc_path} {output.bed} {output.junc}
         """

rule concatenate:
    input:
         expand(rules.bam2junc.output, NAMES=NAMES)
    output:
          junc="leafcutter/{comp_names}/juncfiles.txt",
          test="leafcutter/{comp_names}/diff_introns.txt"
    run:
        comp = output.junc.split("/")[1]
        cond_a, cond_b = comp.split("-vs-")
        with open(output.junc, "w") as fout:
            for name in NAMES:
                if name.startswith(cond_a) or name.startswith(cond_b):
                    fout.write("leafcutter/{}.junc\n".format(name))

        with open(output.test, "w") as fout:
            for name in NAMES:
                if name.startswith(cond_a):
                    fout.write("{} {}\n".format(
                        name, cond_a))

                elif name.startswith(cond_b):
                    fout.write("{} {}\n".format(
                        name, cond_b))

# step 2
# m= min numb reads per cluster, l = intron length
rule intron_clustering:
    input: rules.concatenate.output.junc
    params:
          m=50,
          l=500000,
          prefix="leafcutter/{comp_names}/{comp_names}",
          n="{comp_names}",
          strand="--strand" if config.get("strandness") else ""
    output:
          "leafcutter/{comp_names}/{comp_names}_perind_numers.counts.gz"
    shell:
         """
         python2 scripts/leafcutter_cluster.py \
         -j {input} -m {params.m} -o {params.prefix} -l {params.l} {params.strand}
         rm *{params.n}.sorted.gz
         """

# step 3.1
rule gtf_to_exon:
    input:
         gtf_path
    output:
          a="leafcutter/" + basename(gtf_path, suffix=".gz"),
          b="leafcutter/exons.gtf.gz"
    params:
        gtf_to_exon = srcdir("scripts/gtf_to_exons.R")
    shell:
         """
         gzip -c {input} > {output.a}
         Rscript {params.gtf_to_exon} {output.a} {output.b}
         """

# step 3.2
rule differential_splicing:
    input:
         a=rules.gtf_to_exon.output.b,
         b="leafcutter/{comp_names}/{comp_names}_perind_numers.counts.gz",
         c="leafcutter/{comp_names}/diff_introns.txt"
    output:
          "leafcutter/{comp_names}/{comp_names}_cluster_significance.txt"
    params:
          min_samples_per_group=config["min_samples_per_group"],
          min_samples_per_intron=config["min_samples_per_intron"],
          prefix="leafcutter/{comp_names}/{comp_names}",
          leafcutter_ds_path=srcdir("scripts/leafcutter_ds.R")
    threads: 1-
    shell:
         """
         Rscript {params.leafcutter_ds_path} --exon_file={input.a} \
         {input.b} {input.c} --num_threads {threads} --output_prefix={params.prefix} \
         -i {params.min_samples_per_intron} -g {params.min_samples_per_group}
         """
