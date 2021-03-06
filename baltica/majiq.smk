# -*- coding: utf-8
"""
Created on 17:07 29/02/2018 2018
Snakemake workflow for Majiq.

If you use Majiq, please cite

Vaquero-Garcia, Jorge, et al. "A new view of transcriptome complexity and
regulation through the lens of local splicing variations." elife 5 (2016):
e11752.

.. usage:
    snakemake -s majiq.smk --configfile config.yaml -j 10
"""
__author__ = "Thiago Britto Borges"
__copyright__ = "Copyright 2020, Dieterichlab"
__email__ = "Thiago.BrittoBorges@uni-heidelberg.de"
__license__ = "MIT"
try:
    import baltica
    baltica_installed = True
except ImportError:
    baltica_installed = False

def dir_source(script, ex):
    return script if baltica_installed else srcdir(f"{ex} ../scripts/{script}")

def comparison(wc, index, mapping):
    condition = wc.contrast.split("-vs-")[index]
    return expand("majiq/{name}.majiq", name=mapping[condition])


def natural_sort_key(s, _nsre=re.compile("([0-9]+)")):
    return [int(text) if text.isdigit() else text.lower()
            for text in _nsre.split(s)]


workdir: config.get("path", ".")
name = config["samples"].keys()
sample = config["samples"].values()
sample_path = config["sample_path"]
contrasts = config["contrasts"]
conditions = sorted(
    set([x.split("_")[0] for x in name]),
    key=natural_sort_key)

mapping = {c: [x for x in name if x[: x.index("_")] == c]
           for c in conditions}
localrules: symlink, create_ini

if 'majiq_env_prefix' in config:
    shell.prefix(config["majiq_env_prefix"])

include: "symlink.smk"

rule all:
    input:
         'logs/',
         expand("mappings/{name}.bam", name=name),
         "majiq/build.ini",
         expand("majiq/{name}.majiq", name=name),
         expand("majiq/{contrast}/{contrast}.deltapsi.voila", contrast=contrasts.keys()),
         expand("majiq/voila/{contrast}_voila.tsv", contrast=contrasts.keys())


rule create_ini:
    input: expand("mappings/{name}.bam", name=name)
    output: "majiq/build.ini"
    run:
        lines = [
            "[info]",
            "readlen={}".format(config["read_len"]),
            "bamdirs={}".format("mappings/"),
            "genome={}".format(config["assembly"]),
            "genome_path={}".format(config["ref_fa"]),
            "strandness={}".format(config.get("strandness", 'None')),                      "[experiments]"]

        lines.extend(
            ["{}={}".format(k, ",".join(v))
             for k, v in mapping.items()])

        with open(str(output), "w") as ini:
            ini.writelines("\n".join(lines))


rule gtf_to_gff:
    output: "majiq/ref.gff"
    params: ref=config["ref"],
            gtf2gff3_path=dir_source("gtf2gff3.pl", 'perl')
    shadow: "shallow"
    shell: "{params.gtf2gff3_path} {params.ref} > {output}"


rule build:
    input: ini="majiq/build.ini",
         ref="majiq/ref.gff",
    output: expand("majiq/{name}.majiq", name=name),
          "majiq/splicegraph.sql"
    conda: "../envs/majiq.yml"
    threads: len(conditions)
    envmodules: "majiq/2.2"
    shadow: "shallow"
    shell: " majiq build --conf {input.ini} --nproc {threads} --output majiq/ {input.ref}"


rule deltapsi:
    input: a=lambda wc: comparison(wc, 0, mapping),
         b=lambda wc: comparison(wc, -1, mapping)
    output: "majiq/{contrast}/{contrast}.deltapsi.voila"
    conda: "../envs/majiq.yml"
    threads: 10
    params: name=lambda wc: wc.contrast.replace("-vs-", " "),
          name2=lambda wc: wc.contrast.replace("-vs-", "_"),
          cont=lambda wc: wc.contrast,
          majiq_minreads=config.get('minreads', 3)
    envmodules: "majiq/2.2"
    shadow: "shallow"
    shell: "majiq deltapsi -grp1 {input.a} -grp2 {input.b} " 
           "--nproc {threads} --output majiq/{params.cont} "
           "--names {params.name} "
           "--minreads {params.majiq_minreads} ;"
           "mv majiq/{params.cont}/{params.name2}.deltapsi.voila "
           "{output} "


rule voila:
    input: "majiq/splicegraph.sql",
         "majiq/{contrast}/{contrast}.deltapsi.voila"
    output: "majiq/voila/{contrast}_voila.tsv"
    conda: "../envs/majiq.yml"
    envmodules: "majiq/2.2"
    params: 
        threshold=config.get('majiq_threshold', 0.2),
        non_changing_threshold=config.get('majiq_non_changing_threshold', 0.05)
    shadow: "shallow"
    shell: "voila tsv " 
           "--threshold {params.threshold} "
           "--non-changing-threshold {params.non_changing_threshold} "
           "{input} -f {output}"