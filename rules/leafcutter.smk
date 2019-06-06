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

from pathlib import Path

def basename(path, suffix=None):
    if suffix:
        return str(Path(path).with_suffix(suffix).name)
    return str(Path(path).name)

configfile: "config.yml"
NAMES = config["samples"].keys()
SAMPLES = config["samples"].values()
gtf_path = config["ref"]
conditions = [x.split('_')[0] for x in NAMES]
comp_names = config['contrasts'].keys()

localrules: all, concatenate


rule all:
    input:
        expand('leafcutter/{comp_names}/{comp_names}_cluster_significance.txt',
         comp_names=comp_names),
        expand('leafcutter/{NAMES}.junc', NAMES=NAMES)

rule clean:
    shell:
        'rm -rf leafcutter'

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

# step 1.1
rule bam2junc:
    input: 'mappings/{NAMES}.bam'
    output:
        bed='leafcutter/{NAMES}.bed',
        junc='leafcutter/{NAMES}.junc'
    shell:
        """
        samtools view {input} \
        | python2 scripts/filter_cs.py \
        | perl scripts/sam2bed.pl --use-RNA-strand - {output.bed}
        perl scripts/bed2junc.pl {output.bed} {output.junc}
        """

rule concatenate:
    input:
        expand(rules.bam2junc.output, NAMES=NAMES)
    output:
        junc='leafcutter/{comp_names}/juncfiles.txt',
        test='leafcutter/{comp_names}/diff_introns.txt'
    run:
        comp = output.junc.split('/')[1]
        cond_a, cond_b = comp.split('-vs-')
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
        prefix='leafcutter/{comp_names}/{comp_names}',
        n='{comp_names}'
    output:
        'leafcutter/{comp_names}/{comp_names}_perind_numers.counts.gz'
    shell:
        """
        python2 scripts/leafcutter_cluster.py  \
        -j {input} -m {params.m} -o {params.prefix} -l {params.l} --strand True
        rm *{params.n}.sorted.gz
        """

# step 3.1
rule gtf_to_exon:
    input:
        gtf_path
    output:
        a='leafcutter/' + basename(gtf_path, suffix='.gz'),
        b='leafcutter/exons.gtf.gz'
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
    params:
        min_samples_per_group=config['min_samples_per_group'],
        min_samples_per_intron=config['min_samples_per_intron'],
        prefix='leafcutter/{comp_names}/{comp_names}'
    threads: 4
    shell:
        """
        Rscript scripts/leafcutter_ds.R --exon_file={input.a} \
        {input.b} {input.c} --num_threads {threads} --output_prefix={params.prefix} \
        -i {params.min_samples_per_intron} -g {params.min_samples_per_group}
        """