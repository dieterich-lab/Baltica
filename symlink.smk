# -*- coding: utf-8
"""
Created on 17:07 29/02/2018 2018
Snakemake file for majiq.
.. usage:
    snakemake -s symlink.Snakemake


"""
__author__ = "Thiago Britto Borges"
__copyright__ = "Copyright 2018, Dieterichlab"
__email__ = "Thiago.BrittoBorges@uni-heidelberg.de"
__license__ = "MIT"


configfile: "config.yml"
NAMES = config["samples"].keys()
SAMPLES = config["samples"].values()

localrules: symlink

rule symlink:
    input:
        bam=expand("{SAMPLES}", SAMPLES=SAMPLES),
        bai=expand("{SAMPLES}.bai", SAMPLES=SAMPLES)
    output:
        bam=expand('mappings/{NAMES}.bam', NAMES=NAMES),
        bai=expand('mappings/{NAMES}.bam.bai', NAMES=NAMES)
    run:
        for bam_in, bai_in, bam_out, bai_out in zip(
            input.bam, input.bai, output.bam, output.bai):
            os.symlink(bam_in, bam_out)
            os.symlink(bai_in, bai_out)
