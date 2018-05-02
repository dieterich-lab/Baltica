# -*- coding: utf-8
"""
Created on 17:07 29/02/2018 2018
Snakemake file for leafcutter.
.. install:
    if (!require("devtools")) install.packages("devtools", repos='http://cran.us.r-project.org')
devtools::install_github("davidaknowles/leafcutter/leafcutter")

.. usage:
    sbatch submit_smk.sh leafcutter.smk

Notes:
    They recommend Olego as mapper, but STAR is fine, 
    for example:
    STAR --genomeDir hg19index/
         --twopassMode Basic
         --outSAMstrandField intronMotif
         --readFilesCommand zcat
         --outSAMtype BAM Unsorted
"""
__author__ = "Thiago Britto Borges"
__copyright__ = "Copyright 2018, Dieterichlab"
__email__ = "Thiago.BrittoBorges@uni-heidelberg.de"
__license__ = "MIT"

from itertools import combinations
from pathlib import Path

def basename(path, suffix=None):
    if suffix:
        return str(Path(path).with_suffix(suffix).name)
    return str(Path(path).name)

configfile: "config.yml"
NAMES = config["samples"].keys()
SAMPLES = config["samples"].values()
bin_path = '/home/tbrittoborges/bin/leafcutter'
gtf_path = config["gtf_path"]
conditions = [x.split('_')[0] for x in NAMES]
comp_names = ['{}_vs_{}'.format(*x)
    for x in combinations(sorted(set(conditions)), 2)]

localrules: all, concatenate


rule all:
    input:
        expand('leafcutter/{comp_names}/ds_plots.pdf', comp_names=comp_names),
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
    params: bin_path=bin_path
    shell:
        """
        module load samtools
        samtools view {input} \
        | python2 {params.bin_path}/scripts/filter_cs.py \
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
    shell:
        """
        python2 {bin_path}/clustering/leafcutter_cluster.py \
        -j {input} -m {params.m} -o {params.prefix} -l {params.l}
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
        module load R
        Rscript {bin_path}/scripts/gtf_to_exons.R {output.a} {output.b}
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
    params:
        fdr=config['fdr']
    shell:
        """
        module load R
        Rscript {bin_path}/scripts/ds_plots.R -e {input.exons} \
        {input.counts} {input.test} {input.signif} --output={output} \
        {params.fdr}"""
