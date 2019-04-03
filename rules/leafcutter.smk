# -*- coding: utf-8
"""
Created on 17:07 29/02/2018 2018
Snakemake file for leafcutter.
.. installation:

.. usage:

.. notes:

"""
__author__ = "Thiago Britto Borges"
__copyright__ = "Copyright 2018, Dieterichlab"
__email__ = "Thiago.BrittoBorges@uni-heidelberg.de"
__license__ = "MIT"

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
print(mapping)
localrules: symlink, create_ini


rule all:
    input:
        expand('leafcutter/{comp_names}/ds_plots.pdf', comp_names=comp_names),
        expand('leafcutter/{NAMES}.junc', NAMES=NAMES)


rule symlink:
  input: expand(join(sample_path, "{raw_name}"), raw_name=raw_name)
  output: expand("mappings/{name}.bam", name=name)
  run:
    for i, o in zip(input, output):
      os.symlink(i, o)


# step 1.1
rule bam2junc:
    input: 'mappings/{NAMES}.bam'
    output:
        bed='leafcutter/{NAMES}.bed',
        junc='leafcutter/{NAMES}.junc'
    env: "../envs/leafcutter.yml"
    shell:
        """
        samtools view {input} \
        | python2 scripts/filter_cs.py \
        | perl {params.bin_path}/scripts/sam2bed.pl --use-RNA-strand - {output.bed}
        perl {params.bin_path}/scripts/bed2junc.pl {output.bed} {output.junc}
        """

rule concatenate:
    input:
        expand(rules.bam2junc.output, NAMES=NAMES)
    output:
        junc='leafcutter/{comp_names}/juncfiles.txt',
        test='leafcutter/{comp_names}/diff_introns.txt'
    run:
        comp = output.junc.split('/')[1]
        cond_a, cond_b = comp.split('_vs_')
        with open(output.junc, 'w') as fout:
            for name in NAMES:
                if name.startswith(cond_a) or name.startswith(cond_b):
                    fout.write("leafcutter/{}.junc\n".format(name))

        with open(output.test, 'w') as fout:
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
        prefix='leafcutter/{comp_names}/{comp_names}'
    output:
        'leafcutter/{comp_names}/{comp_names}_perind_numers.counts.gz'
    env: "../envs/leafcutter.yml"
    shell:
        """
        python2 scripts/leafcutter_cluster.py \
        -j {input} -m {params.m} -o {params.prefix} -l {params.l}
        """

# step 3.1
rule gtf_to_exon:
    input:
        gtf_path
    output:
        a='leafcutter/' + basename(gtf_path, suffix='.gz'),
        b='leafcutter/exons.gtf.gz'
    env: "../envs/leafcutter.yml"
    shell:
        """
        gzip -c {input} > {output.a}
        Rscript scripts/gtf_to_exons.R {output.a} {output.b}
        """

# step 3.2
rule differential_splicing:
    input:
        a=rules.gtf_to_exon.output.b,
        b='leafcutter/{comp_names}/{comp_names}_perind_numers.counts.gz',
        c='leafcutter/{comp_names}/diff_introns.txt'
    output:
        'leafcutter/{comp_names}/{comp_names}_cluster_significance.txt'
    env: "../envs/leafcutter.yml"
    params:
        min_samples_per_group=config['min_samples_per_group'],
        min_samples_per_intron=config['min_samples_per_intron'],
        prefix='leafcutter/{comp_names}/{comp_names}'
    threads: 4
    shell:
        """
        module load R
        Rscript {bin_path}/scripts/leafcutter_ds.R --exon_file={input.a} \
        {input.b} {input.c} --num_threads {threads} --output_prefix={params.prefix} \
        {params.min_samples_per_intron} {params.min_samples_per_group}
        """

rule plot:
    input:
        exons=rules.gtf_to_exon.output.b,
        test=rules.concatenate.output.test,
        counts=rules.intron_clustering.output,
        signif=rules.differential_splicing.output,
    output:
        'leafcutter/{comp_names}/ds_plots.pdf'
    env: "../envs/leafcutter.yml"
    params:
        fdr=config['fdr']
    shell:
        """
        Rscript {bin_path}/scripts/ds_plots.R -e {input.exons} \
        {input.counts} {input.test} {input.signif} --output={output} \
        {params.fdr}"""
