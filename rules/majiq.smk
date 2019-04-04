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

__author__ = "Thiago Britto Borges"
__copyright__ = "Copyright 2019, Dieterichlab"
__email__ = "Thiago.BrittoBorges@uni-heidelberg.de"
__license__ = "MIT"

from os.path import join
import re

def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
    return [int(text) if text.isdigit() else text.lower()
            for text in _nsre.split(s)]


def comparison(wc, index, mapping):
    condition = wc.contrast.split("-")[index]
    return expand("majiq/{name}.majiq", name=mapping[condition])


configfile: "config.yml"
name = config["samples"].keys()
raw_name = config["samples"].values()
sample_path = config["sample_path"]
contrasts = config['contrasts']

conditions = sorted(
    set([x.split("_")[0] for x in name]),
    key=natural_sort_key)
mapping = {c: [x for x in name if x[: x.index('_')] == c ]
           for c in conditions}

localrules: symlink, create_ini

rule all:
  input:
    "majiq/build.ini",
    expand("mappings/{name}.bam", name=name),
    expand("majiq/{name}.majiq", name=name),
    expand("majiq/{contrast}/",
        contrast=contrasts.keys()),
    expand("majiq/{contrast}/voila.tsv",
        contrast=contrasts.keys())


rule symlink:
  input:
    expand(join(sample_path, "{raw_name}"),
        raw_name=raw_name)
  output:
    expand("mappings/{name}.bam", name=name)
  run:
    for i, o in zip(input, output):
      os.symlink(i, o)
      os.symlink(i + '.bai', o + '.bai')


rule create_ini:
  input:
    expand("mappings/{name}.bam", name=name)
  output:
    "majiq/build.ini"
  run:
    lines = [
        "[info]",
        "readlen={}".format(config["read_len"]),
        "samdir={}".format("mappings/"),
        "genome={}".format(config["assembly"]),
        "genome_path={}".format(config["ref_fa"]),
        "strandness={}".format(config["strandness"]),
        "[experiments]"]

    lines.extend(
        ["{}={}".format(k, ",".join(v))
        for k, v in mapping.items()])

    with open(str(output), "w") as ini:
        ini.writelines("\n".join(lines))


rule gtf_to_gff:
    output:
        "majiq/ref.gff"
    params:
        ref = config["ref"],
        exe = "perl scripts/gtf2gff3.pl"

    shell:
        "{params.exe} {params.ref} > {output}"


rule build:
    input:
        ini="majiq/build.ini",
        ref="majiq/ref.gff"
    output:
        expand("majiq/{name}.majiq", name=name)
    conda:
        "../envs/majiq.yml"
    threads: len(conditions)
    shell:
        " majiq build --conf {input.ini} --nproc {threads} "
        " --output majiq/ {input.ref}"


rule deltapsi:
    input:
        a=lambda wc: comparison(wc, 0, mapping),
        b=lambda wc: comparison(wc, -1, mapping)
    output:
        directory("majiq/{contrast}/")
    conda:
        "../envs/majiq.yml"
    threads:
        10
    params:
        name=lambda wc: wc.contrast.replace('-vs-', ' '),
    shell:
        "majiq deltapsi -grp1 {input.a} -grp2 {input.b} "
        "--nproc {threads} --output {output} "
        "--names {params.name} --default-prior "
