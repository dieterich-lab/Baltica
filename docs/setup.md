# Getting started

Baltica comprise of a collection of workflows and analysis scripts. Workflows are powered by [Snakemake](https://snakemake.readthedocs.io/en/stable/) [^1]. Analysis are done with the Rlang. Bellow we document how to obtain install the methods on which Baltica depends.    

## Install miniconda

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```  

Or get miniconda [here](https://docs.conda.io/en/latest/miniconda.html) 

Follow the instructions to finish your installation, which by default is at `$HOME/miniconda3`

Make sure you initialize conda with `conda init`
You can test whether your installation was successful or not by running `conda --version` and you may need to restart 
your shell instance. 

## Clone Baltica

```bash
git clone git@github.com:dieterich-lab/baltica.git
```

Baltica can be used with the [modules system](https://modules.readthedocs.io/en/latest/index.html) or conda environments. Here we describe the installation with conda. 

## Install Snakemake 
```bash
conda install -c bioconda snakemake==5.2 --yes
```

Some of dependencies can be directly created with conda, go to `baltica/envs` directory and install with:  
    - `conda env create -f stringtie.yml --yes`  
    - `conda env create -f qc.yml --yes`  

In general, R packages do not play nicely with conda, but we still use it because it's flexibility and the ability to 
create isolated software environments.
Cite Stringtie[^5].
<!-- TODO cite rseqc fastqc and multiqcgit clone https://github.com/comwes/mkpdfs-design-sample -->


## Install Majiq 

!!! warning
    Majiq requires a Academic or Commercial license for use. Users are required to obtain their license. [Academic download](https://majiq.biociphers.org/app_download/).

Majiq[^2] can installation can be problematic, but the recipe bellow works for us:

```bash
conda create --name majiq_env python=3.6 pysam numpy cython --yes -c bioconda
conda activate majiq_env
conda install --yes waitress==1.1.0 h5py>=2.8.0 Flask==1.0.2 Flask-WTF==0.14.2 GitPython>=2.1.11 gunicorn==19.9.0 psutil>=5.4.8 h5py>=2.8.0 scipy>=1.1.0
conda install --yes -c bioconda htslib 

python3 -m pip install --upgrade pip
python3 -m pip install git+https://bitbucket.org/biociphers/majiq_stable.git#egg=majiq
```

## Installation Leafcutter

Users can install Leafcutter[^3] with conda using the following recipe: 

```bash
conda create --name leafcutter python=2.7 --yes
conda activate leafcutter
conda install -c bioconda samtools r-base=3.6 --yes

Rscript -e "install.packages('devtools', repos='http://cran.us.r-project.org', dependencies=TRUE)"
Rscript -e "Sys.setenv(TAR = '/bin/tar'); devtools::install_github('stan-dev/rstantools')"
Rscript -e "Sys.setenv(TAR = '/bin/tar'); devtools::install_github('davidaknowles/leafcutter/leafcutter')"
```
<!-- R CMD INSTALL --no-lock <pkg> -->

!!! Important
    Only use `Sys.setenv(TAR = '/bin/tar')` if you have problems with devtools selecting `gtar` instead of `tar`


!!! Important
    If you are experiencing the following the `ERROR: failed to create lock directory` error when trying to install R packages add the following option to install.package `INSTALL_opts = c('--no-lock')`

## Install Junctionseq

JunctionSeq[^4] should be installed directly from BioConductor

```bash
conda env create -f junctionseq-env.yml --yes
conda activate leafcutter
Rscript -e "BiocManager::install('JunctionSeq')"
```

## Clone or installing Baltica?
Baltica can either installed as python package or cloned from Github for each project.
Users who intend to modify the workflows should clone the framework and keep the change under version.

See (workflows)[workflows.md] for details on the configuration and parameters for each available workflow.


[^1]: If you use Baltica, please also [cite Snakemake](https://bioinformatics.oxfordjournals.org/content/28/19/2520)
[^2]: If you use Majiq results, please [cite it]( https://elifesciences.org/articles/11752)
[^3]: If you use Leafcutter results, please [cite it](https://www.nature.com/articles/s41588-017-0004-9)
[^4]: If you use Junctionseq results, please [cite it](http://nar.oxfordjournals.org/content/early/2016/06/07/nar.gkw501.full)
[^5]: If you use the Baltica's analysis module, please also [cite Stringtie](http://www.nature.com/nbt/journal/vaop/ncurrent/full/nbt.3122.html)
