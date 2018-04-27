#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=4G
#SBATCH --job-name="smk"
#SBATCH --mail-user=thiago.brittoborges@uni-heidelberg.de
#SBATCH --partition=long

snakemake -s $1 --keep-going --cluster \
    'sbatch ' --cluster-config  \
    cluster.json --jobs 100