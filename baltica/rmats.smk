# -*- coding: utf-8
"""
Created on 17:07 27/07/2018
Snakemake workflow for rMATS

If you use Majiq, please cite

Shen, Shihao, et al. "rMATS: robust and flexible detection of differential 
alternative splicing from replicate RNA-Seq data." Proceedings of the 
National Academy of Sciences 111.51 (2014): E5593-E5601.
"""
__author__ = "Thiago Britto Borges"
__copyright__ = "Copyright 2020, Dieterichlab"
__email__ = "Thiago.BrittoBorges@uni-heidelberg.de"
__license__ = "MIT"

from collections import defaultdict
import tempfile


container: "docker://tbrittoborges/rmats:latest"


workdir: config.get("path", ".")


contrasts = config["contrasts"]
keys = config["samples"].keys()
keys = [tuple(x.split("_")) for x in keys]

temp_dir = tempfile.TemporaryDirectory()
strand = {"reverse": "fr-secondstrand", "forward": "fr-firststrand"}

d = defaultdict(list)
for x in keys:
    d[x[0]].append(x[1])


include: "symlink.smk"


localrules:
    create_rmats_input,


rule all:
    input:
        expand("rmats/{group}.txt", group=d.keys()),
        expand("rmats/{contrast}/", contrast=contrasts),


rule rmats_create_input:
    input:
        lambda wc: expand("mappings/{{group}}_{rep}.bam", rep=d.get(wc.group)),
    log:
        "logs/rmats/create_input_{group}.log",
    output:
        "rmats/{group}.txt",
    run:
        with open(str(output), "w") as fou:
            fou.write(",".join(input))


rule rmats_run:
    input:
        alt="rmats/{alt}.txt",
        ref="rmats/{ref}.txt",
    output:
        directory("rmats/{alt}-vs-{ref}/"),
    shadow:
        "shallow"
    log:
        "logs/rmats/run_{alt}-vs-{ref}.log",
    threads: 10
    envmodules:
        "rmats-turbo/4.1.1",
    params:
        gtf=config["ref"],
        is_paired="-t single" if config.get("is_single_end") else "",
        lib="--libType " + strand.get(config["strandness"], "fr-unstranded"),
        read_len=config["read_len"],
        allow_clipping="--allow-clipping",
        tmp=os.path.join(temp_dir.name, "{alt}_vs_{ref}/"),
    shell:
        "rmats.py "
        "--b1 {input.alt} "
        "--b2 {input.ref} "
        "--gtf {params.gtf} "
        "--variable-read-length "
        "--readLength {params.read_len} "
        "--nthread {threads} "
        "--novelSS "
        "{params.lib} "
        "{params.allow_clipping} "
        "--od {output} "
        "--tmp {params.tmp}"
