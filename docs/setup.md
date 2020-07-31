## Getting started

Baltica contains a collection of workflows and analysis scripts. Workflows are powered by [Snakemake](https://snakemake.readthedocs.io/en/stable/) [^1]. Analysis is done with the R using Bioconductor packages. 
We developed and tested the workflows with the Debian Linux distribution (v8.11 Jesse) and conda (4.8.3).
We use the module system to test the workflows, but conda usage is similar. 
Below, we document how to install the Baltica dependencies. 

## Install miniconda

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Or get miniconda [here](https://docs.conda.io/en/latest/miniconda.html). 

Follow the instructions to finish then installation, which by default is at `$HOME/miniconda3`.

Make sure you initialize conda with `conda init`.
You can test whether your installation was successful or not by running `conda --version` and you may need to restart your shell instance. 

## Clone and install Baltica

```bash
git clone https://github.com/dieterich-lab/baltica
cd baltica
python setup.py install
```


## Install Snakemake 
```bash
conda install -c bioconda snakemake">=5.2" --yes
```

!!! danger
    Snakemake requires python versions between > 3.4 and < 3.6.

Some of dependencies can be directly created with conda, go to `baltica/envs` directory and install with:  
```
conda env create -f stringtie.yml
```

and

```
conda env create -f qc.yml
```

In general, R packages do not play nicely with conda, but we still use it because it's flexibility and the ability to 
create isolated software environments.

Stringtie citation [^5].

## Install Majiq 

!!! warning
    Majiq requires an Academic or Commercial license for use. Users are required to obtain their license. [Academic download](https://majiq.biociphers.org/app_download/).

Majiq[^2] can installation can be problematic, but the recipe below works for us:

```bash
conda create --name majiq_env python=3.6 htslib "Cython==0.29.14" git pip --yes -c bioconda
conda activate majiq_env

ENV_PATH=$HOME/miniconda3/envs/majiq_env
export HTSLIB_LIBRARY_DIR=$ENV_PATH/lib
export HTSLIB_INCLUDE_DIR=$ENV_PATH/include
pip install git+https://bitbucket.org/biociphers/majiq_academic.git#egg=majiq 
```

!!! danger
    Confirm that `pip` used is the one from the `majiq_env` with `which pip`  

Please inform us if you have issues with this recipe.

## Installation Leafcutter

Users can install Leafcutter[^3] with conda using the following recipe: 

```bash
conda create --name leafcutter python=2.7
conda activate leafcutter
conda install -c bioconda samtools r-base=3.6

TAR='/bin/tar'
Rscript -e "install.packages('devtools', repos='http://cran.us.r-project.org')"
Rscript -e "devtools::install_github('stan-dev/rstantools')"
```

!!! warning
    Only use `TAR='/bin/tar'` or `set TAR '/bin/tar'` (fishshell) if you have problems with devtools selecting `gtar` instead of `tar`.

!!! warning
    If you are experiencing the following the `ERROR: failed to create lock directory` error when trying to install R packages, add the following option to install.package `INSTALL_opts = c('--no-lock')`.

## Install Junctionseq

JunctionSeq[^4] should be installed directly from BioConductor:

```bash
conda env create -f envs/junctionseq.ym
conda activate leafcutter
Rscript -e "BiocManager::install('JunctionSeq')"
```

## Clone or installing Baltica?
Baltica can either installed as a python package or cloned from Github for each project.
The installed version of Baltica is more convenient to be used with:  
`baltica qc config.yml` (as long the dependencies are available).
Users who intend to modify the workflows should clone the framework and keep the change under version. 
See (workflows)[workflows.md] for details on each available workflow configuration and parameters.

## Baltica command-line arguments

Use the command below to list the command line arguments and their options: 
```
baltica --help
```

## Executing a Baltica workflow

* with the [modules system](https://modules.readthedocs.io/en/latest/index.html):
    `baltica qc config.yml --use-envmodule`
* with conda enviroments:
    `baltica qc config.yml --use-conda`
* using an external conda enviroment, like the one we used for Majiq installation:
    set `majiq_env_prefix = conda activate majiq_env;` in the configuration file

There are alternatives to provide the software dependencies to Snakemake workflows, so feel free to contact them if you need an option.

[^1]: If you use Baltica, please also [cite Snakemake](https://bioinformatics.oxfordjournals.org/content/28/19/2520)
[^2]: If you use Majiq results, please [cite it]( https://elifesciences.org/articles/11752)
[^3]: If you use Leafcutter results, please [cite it](https://www.nature.com/articles/s41588-017-0004-9)
[^4]: If you use Junctionseq results, please [cite it](http://nar.oxfordjournals.org/content/early/2016/06/07/nar.gkw501.full)
[^5]: If you use the Baltica's analysis module, please also [cite Stringtie](http://www.nature.com/nbt/journal/vaop/ncurrent/full/nbt.3122.html)

\bibliography