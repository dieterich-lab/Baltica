# -*- coding: utf-8
"""
Created on 17:07 29/02/2018 2018
Snakemake workflow for Majiq.

If you use Majiq, please cite

Vaquero-Garcia, Jorge, et al. "A new view of transcriptome complexity and
regulation through the lens of local splicing variations." elife 5 (2016):
e11752.

.. usage:

"""
from itertools import combinations

__author__ = "Thiago Britto Borges"
__copyright__ = "Copyright 2019, Dieterichlab"
__email__ = "Thiago.BrittoBorges@uni-heidelberg.de"
__license__ = "MIT"

configfile: "../config.yml"
NAMES = config["samples"].keys()
SAMPLES = config["samples"].values()
conditions = sorted(set([x.split("_")[0] for x in NAMES]))
contrasts = ["{}_{}".format(*x) for x in
             combinations(conditions, 2)]

mapping = {c: [x for x in NAMES if x.startswith(c)]
           for c in conditions}


def comparison(wildcards, index, mapping):
    condition = wildcards.comp_names.split("_")[index]
    return expand("majiq/{name}.majiq", name=mapping[condition])



rule all:
    input:
        "build.ini",
        expand("{names}.majiq", names=NAMES),
        directory(
            expand("{contrast}/voila_deltapsi/",
            contrast=contrasts))


rule create_ini:
    input: expand("../mappings/{names}.bam", names=NAMES)
    output: "build.ini"
    run:
        lines = [
            "[info]",
            "readlen={}".format(config["read_len"]),
            "samdir={}".format("../mappings/"),
            "genome={}".format(config["assembly"]),
            "genome_path={}".format(config["ref_fa"]),
            "strandness={}".format(config["strandness"]),
            "[experiments]"]

        lines.extend(
            ["{}={}".format(k, ",".join(v)) for k, v in mapping.items()])

        with open(str(output), "w") as ini:
            ini.writelines("\n".join(lines))


rule gtf_to_gff:
    output:
        "ref.gff"
    params:
        ref = config["ref"],
        exe = "perl /home/tbrittoborges/bin/gtf2gff3.pl"

    shell:
        "{params.exe} {params.ref} > {output}"


rule build:
    input:
        ini="build.ini",
        ref="ref.gff"
    output:
        expand("{names}.majiq", names=NAMES)
    params:
        exe=". /home/tbrittoborges/bin/miniconda3/etc/profile.d/conda.sh; "
        " conda activate majiq_env; "
    threads: 1
    shell:
        " {params.exe} "
        " majiq build --conf {input.ini} --nproc {threads} "
        " --output majiq/k {input.ref}"

rule deltapsi:
    input: cont=lambda wildcards: comparison(wildcards, 0),
           treat=lambda wildcards: comparison(wildcards, 1)
    output: "{comp_names}/{comp_names}.deltapsi.voila"
    threads: 20
    params:
        output="{comp_names}/",
        names=lambda wildcards: wildcards.comp_names.replace("_", " ")
    shell:
        "{params.exe} "
        "majiq deltapsi -grp1 {input.cont} -grp2 {input.treat} "
        "--nproc {threads} --output {params.output} --names {params.names} "


rule voila_deltapsi:
    input: "majiq/{comp_names}/{comp_names}.deltapsi.voila"
    output: directory("majiq/{comp_names}/voila_deltapsi/")
    shell:
        " {params.exe} "
        " voila deltapsi -o {output} "
        " -s majiq/splicegraph.sql {input}"


rule deltapsi:
    input: cont=lambda wildcards: comparison(wildcards, 0),
           treat=lambda wildcards: comparison(wildcards, 1)
    output: "{comp_names}/{comp_names}.deltapsi.voila"
    threads: 20
    params:
        output="{comp_names}/",
        names=lambda wildcards: wildcards.comp_names.replace("_", " ")
    shell:
        " {params.exe} "
        "majiq deltapsi -grp1 {input.cont} -grp2 {input.treat} "
        "--nproc {threads} --output {params.output} --names {params.names} "


rule voila_deltapsi:
    input: "majiq/{comp_names}/{comp_names}.deltapsi.voila"
    output: directory("majiq/{comp_names}/voila_deltapsi/")
    shell:
        " {params.exe} "
        " voila deltapsi -o {output} "
        " -s majiq/splicegraph.sql {input}"





