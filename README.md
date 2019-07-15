# Baltica

One stop solution for differential splicing analysis.

## Quickstart (LeafCutter example)
[Instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) for SnakeMake installlation

1) Clone this repo: `git clone dieterich-lab/Baltica/`

2) [Install LeafCutter](#install-leafcutter)

3) Within Baltica folder fill the `config.yaml` 

4) to run in a slurm cluster, use `sbatch submit_snakemake_onslurm.sh rules/leafcutter.smk` 
5) to run locally `snakemake -s rules/leafcutter.smk -j 4` (using 4 CPUs)
6) otherwise check the other options in the SnakeMake [manual](https://snakemake.readthedocs.io/en/latest/executable.html) 



## Install Leafcutter
```bash
conda env create -n leafcutter python=2.7
conda activate leafcutter
conda install -c bioconda samtools r-base=3.3.3

Rscript -e "if (!require("devtools")) install.packages("devtools", repos='http://cran.us.r-project.org'); devtools::install_github("davidaknowles/leafcutter/leafcutter")"
```

## Install Majiq
```bash
conda create --name majiq python=3.5 pysam numpy cython
conda activate majiq
export HTSLIB_INCLUDE_DIR=$(path_to_env)/envs/MAJIQ/lib/python3.6/site-packages/pysam/include/htslib/
export HTSLIB_LIBRARY_DIR=$(path_to_env)/envs/MAJIQ/lib/python3.6/site-packages/pysam/include/htslib/htslib/

pip install git+https://bitbucket.org/biociphers/majiq_stable.git#egg=majiq
```

## Install JunctionSeq
```bash
conda create --name junctionseq qorts r-biocmanager -c bioconda -c conda-forge --yes  
conda activate junctionseq
Rscript -e "BiocManager::install('JunctionSeq', version = '3.8')"
```

## Instal Whippet

```bash
conda create -n whippet -c conda-forge julia=0.6 --yes 
conda activate whippet
julia -e 'Pkg.add("Whippet"); using Whippet'
# /homes/tbrittoborges/miniconda3/envs/whippet/share/julia/site/v0.6/Whippet/src/
# /home/tbrittoborges/miniconda3/envs/whippet/bin/julia
```
