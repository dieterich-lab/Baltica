#!/usr/bin/env python
# -*- coding: utf-8
"""
Created on 17:07 24/01/2018 2018 
Snakemake file for majiq.
.. usage:
    find . -name majiq.slurm.sh -exec sbatch  {} \;
"""
from numpy.f2py import rules

__author__ = "Thiago Britto Borges"
__copyright__ = "Copyright 2017, Dieterichlab"
__email__ = "Thiago.BrittoBorges@uni-heidelberg.de"
__license__ = "MIT"

from itertools import combinations
import os
from os.path import join as path_join

basedir = "/scratch/rbm20/mnt/ARCHIVEDISKS/archive_disk_2_contents/projects/Michael_Gotthardt_Rbm20/"
# wget https://majiq.biociphers.org/download/gtf2gff3.pl
# perl gtf2gff3.pl /biodb/genomes/mus_musculus/GRCm38_79/GRCm38.79.gtf > GRCm38.79.gff
output_dir = "/scratch/tbrittoborges/majiq_Michael_Gotthardt_Rbm20/"
gff3_annotation = output_dir + 'GRCm38.79.gff'
genome =  "GRCm38"
genome_path='/biodb/genomes/mus_musculus/GRCm38_79/GRCm38_79.fa'
 # for f in $basedir/*_STARmapping/Log.final.out; do echo $f; grep 'Average mapped length' $f; done
readlen = 200
sub_template="""#!/bin/bash
#SBATCH -c 20
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=64G
#SBATCH --job-name="majiq_workflow"
#SBATCH --mail-user=thiago.brittoborges@uni-heidelberg.de
#SBATCH --mail-type=ALL
#SBATCH --output=/scratch/tbrittoborges/majiq_Michael_Gotthardt_Rbm20/logs%j.out

module load majiq
module load samtools

parallel  samtools index ::: {output}/*.bam

cd {output_dir}

majiq build \
-conf {output}/majiq.ini \
--nthreads 20 \
--nogc \
--output {output} \
{gff3_annotation}

majiq deltapsi \
-grp1 {output}/{name1}*.majiq.hdf5 \
-grp2 {output}/{name2}*.majiq.hdf5 \
--output {output}/ \
--names {name1} {name2} \
--default_prior

voila deltapsi \
{output}/{name1}_{name2}.deltapsi.voila \
-o {output} \
--splice-graph {output}/splicegraph.hdf5
"""


data={"male_TiKO": [
        "male_TiKO_004_STARmapping/Aligned.noS.correct.bam",
        "male_TiKO_003_STARmapping/Aligned.noS.correct.bam",
        "male_TiKO_002_STARmapping/Aligned.noS.correct.bam"],
    "male_WT": [
        "male_WT_010_STARmapping/Aligned.noS.correct.bam",
        "male_WT_011_STARmapping/Aligned.noS.correct.bam",
        "male_WT_009_STARmapping/Aligned.noS.correct.bam"],
    "female_WT":[
        "female_WT_016_STARmapping/Aligned.noS.correct.bam",
        "female_WT_014_STARmapping/Aligned.noS.correct.bam",
        "female_WT_015_STARmapping/Aligned.noS.correct.bam"],
    "female_TiKO":[
        "female_TiKO_008_STARmapping/Aligned.noS.correct.bam",
        "female_TiKO_006_STARmapping/Aligned.noS.correct.bam",
        "female_TiKO_005_STARmapping/Aligned.noS.correct.bam"]}

sample =  ["_vs_".join(sample) for sample in combinations(data.keys(), 2)]

rule all:
    input: expand('{sample}/majiq.ini', sample=sample)

rule dir:
    output: expand("{sample}/", sample=sample)
    shell: "mkdir -p {sample}"

rule create_ini:
    input: rules.dir.output
    output: expand('{sample}/majiq.ini', sample=sample),
    run:
        for comparison in input:
            sample_a, sample_b = comparison.replace('/', '').split('_vs_')

            samples_names = [sample_a] * len(data[sample_a]) + [sample_b] * len(data[sample_b])
            samples_names = ["{}_{}".format(x, i) for i, x in enumerate(samples_names)]
            samples_paths = data[sample_a] + data[sample_b]


            for filepath, name in zip(samples_paths, samples_names):
                os.symlink(
                    path_join(basedir, filepath),
                    path_join(output_dir, comparison, '{}.bam'.format(name)))


            lines = [
                '[info]',
                'readlen={}'.format(readlen),
                'samdir={}'.format(output_dir + comparison),
                'genome={}'.format(genome),
                'genome_path={}'.format(genome_path),
                '[experiments]',
                '{}={}'.format(sample_a, ','.join(
                    [x for x in samples_names if x.startswith(sample_a)])),
                '{}={}'.format(sample_b, ','.join(
                    [x for x in samples_names if x.startswith(sample_b)]))]

            with open('{0}/majiq.ini'.format(comparison), 'w') as ini:
              ini.writelines('\n'.join(lines))

            with open('{0}/majiq.slurm.sh'.format(comparison), 'w') as sub:
                sub.write(sub_template.format(
                                    output_dir=output_dir,
                                    output=comparison.replace('/', ''),
                                    gff3_annotation=gff3_annotation,
                                    name1=sample_a,
                                    name2=sample_b))