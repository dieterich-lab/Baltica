# -*- coding: utf-8
"""
Created on 17:07 29/02/2018 2018
Snakemake file for majiq.

python3 -m venv majiq_env
source pysster_venv/bin/activate

.. usage:
    snakemake -s majiq.Snakemake --cluster "sbatch --mem=24000 --job-name majiq_pipeline" --jobs 100
"""
__author__ = "Thiago Britto Borges"
__copyright__ = "Copyright 2018, Dieterichlab"
__email__ = "Thiago.BrittoBorges@uni-heidelberg.de"
__license__ = "MIT"




from itertools import combinations, chain
from pathlib import Path

configfile: "config.yml"
NAMES = config["samples"].keys()
SAMPLES = config["samples"].values()
conditions = sorted(set([x.split('_')[0] for x in NAMES]))
comp_names =  ['{}_{}'.format(*x) for x in
               combinations(conditions, 2)]
mapping = {c: [x for x in NAMES if x.startswith(c)] for c in conditions}

def basename(path, suffix=None):
    if suffix:
        return str(Path(path).with_suffix(suffix).name)
    return str(Path(path).name)


def comparison(wildcards, index):
    condition = wildcards.comp_names.split('_')[index]
    return expand('majiq/{name}.majiq.hdf5', name=mapping[condition])


rule all:
    input:
        'majiq/build.ini',
        expand('majiq/{names}.majiq.hdf5', names=NAMES),
        expand('majiq/{comp_names}/{comp_names}.deltapsi_deltapsi.tsv',
            comp_names=comp_names)

rule clean:
    shell:
        "rm -rf majiq"

rule symlink:
    input:
        bam=expand("{SAMPLES}", SAMPLES=SAMPLES),
        bai=expand("{SAMPLES}.bai", SAMPLES=SAMPLES)
    output:
        bam=expand('mappings/{names}.bam', names=NAMES),
        bai=expand('mappings/{names}.bam.bai', names=NAMES)
    run:
        for bam_in, bai_in, bam_out, bai_out in zip(
            input.bam, input.bai, output.bam, output.bai):
            os.symlink(bam_in, bam_out)
            os.symlink(bai_in, bai_out)


rule create_ini:
    input: expand('mappings/{names}.bam', names=NAMES)
    output: 'majiq/build.ini'
    run:
        lines = [
            '[info]',
            'readlen={}'.format(config['read_len']),
            'samdir={}'.format('mappings/'),
            'genome={}'.format(config['assembly']),
            'genome_path={}'.format(config['genome_path']),
            '[experiments]']

        lines.extend(
            ['{}={}'.format(k, ','.join(v)) for k, v in mapping.items()])

        with open(str(output), 'w') as ini:
            ini.writelines('\n'.join(lines))

rule build:
    input:
        'majiq/build.ini'
    output:
        expand('majiq/{names}.majiq.hdf5', names=NAMES)
    params:
        output='majiq/',
        annotation=config['gff_path']
    threads: 20
    shell:
        '''
        module load majiq/1.0.6
        majiq build -conf {input} --nthreads {threads} --nogc \
        --output {params.output} {params.annotation}
        '''


rule deltapsi:
    input: cont=lambda wildcards: comparison(wildcards, 0),
           treat=lambda wildcards: comparison(wildcards, 1)
    output: 'majiq/{comp_names}/{comp_names}.deltapsi.voila'
    threads: 20
    params:
        output='majiq/{comp_names}/',
        names=lambda wildcards: wildcards.comp_names.replace('_', ' '),  # RNPS1 Luc
    shell:
        '''
        module load majiq/1.0.6
        majiq deltapsi -grp1 {input.cont} -grp2 {input.treat} \
        --nthreads {threads} --output {params.output} --names {params.names}
        '''


rule voila_deltapsi:
    input: rules.deltapsi.output
    output: 'majiq/{comp_names}/{comp_names}.deltapsi_deltapsi.tsv'
    params:
        output='majiq/{comp_names}/',
    shell:
        '''
        module load majiq/1.0.6
        voila deltapsi --threshold 0.1 -o {params.output} --splice-graph majiq/splicegraph.hdf5 \
         {input}
        '''
