# Baltica

Reconciling alternative splicing identification

## Features

Based on the results from different tools, Baltica can:

- Document method usage and SnakeMake workflows   
- Process and integrate the results  
- Use data integration methods to detect consequences of AS  

# Installation

0) optional create a Baltica environment:

```bash
conda create --name baltica python=3.6 snakemake PyYAML
```

1) Clone and install Baltica:  
 ```bash
git clone dieterich-lab/Baltica/  
cd Baltica; pip install .
```
1) [Install LeafCutter](#install-leafcutter)
1) [Install Majiq](#install-majiq)
1) [Install JunctionSeq](#install-junctionseq)
1) [Install Whippet](#instal-whippet)

The installation of every supported tools is optional. For detail on the SnakeMake installation, please see
 [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

# Configuration

For project specific configuration, please see [Instructions](#Instruction_for_the_configuration_file)

1) Executing 
1) Analysis
1) Contributing

# Example Usage


# Software requirement installation

I recommend restarting the shell after each successful installation.

## Using conda

1) Install miniconda
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
``` 
Follow the instructions in the screen to finish your installation, which by default is at /home/username/miniconda3

Make sure you initialize conda with `conda init`
You can test whether your installation was successful or not by running `conda --version`
and you may start a new shell instance. 

2) Install SnakeMake
As of today, SnakeMake requires python 3.6

Install with 

```bash
conda install python=3.6 --yes
conda install snakemake -c bioconda --yes
```
and test with `snakemake --version`

In general, R packages do not play nicely with conda, but we still use it because it's flexibility and the ability to create isolated software environments.

3) Install LeafCutter

We can create conda enviroments with LeafCutter dependencies: 
Leafcutter is distributed with an [Apache License 2.0](https://github.com/davidaknowles/leafcutter/blob/master/LICENSE)

```bash
conda create --name leafcutter python=2.7 --yes
conda activate leafcutter
conda install -c bioconda samtools r-base=3.6 --yes
```
Or simply `conda env create -f leafcutter-env.yml`

And next, install it from github:
```
Rscript -e "install.packages('devtools', repos='http://cran.us.r-project.org', dependencies=TRUE, INSTALL_opts = c('--no-lock'))"
Rscript -e "Sys.setenv(TAR = '/bin/tar'); devtools::install_github('stan-dev/rstantools')"
Rscript -e "Sys.setenv(TAR = '/bin/tar'); devtools::install_github('davidaknowles/leafcutter/leafcutter')"

```

5) Install JunctionSeq


```bash
conda env create -f junctionseq-env.yml
Rscript -e "BiocManager::install('JunctionSeq',  INSTALL_opts = c('--no-lock'))"
```

4) Install Majiq
Majiq installation requires an license that can be obtained here:
[Academic download](https://majiq.biociphers.org/app_download/).   
Majiq have multiple licenses that should reviewed at the developers website.

```bash
conda create --name majiq_env python=3.6 pysam numpy cython --yes -c bioconda
conda activate majiq_env
conda install --yes waitress==1.1.0 h5py>=2.8.0 Flask==1.0.2 Flask-WTF==0.14.2 GitPython>=2.1.11 gunicorn==19.9.0 psutil>=5.4.8 h5py>=2.8.0 scipy>=1.1.0
conda install --yes -c bioconda htslib 

python3 -m pip install --upgrade pip
python3 -m pip install git+https://bitbucket.org/biociphers/majiq_stable.git#egg=majiq
```

6) Install Stringtie
```bash
conda env create -f stringtie-env.yml
```

## Using [Environment Modules](https://modules.readthedocs.io/en/latest/index.html) for requirements 

1) LeafCutter

```bash
module load R/3.6 samtools python/2.7

export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8
export R_LIBS_USER=$(realpath ~/R)

Rscript -e "install.packages('devtools', repos='http://cran.us.r-project.org', dependencies=TRUE, INSTALL_opts = c('--no-lock'))"
Rscript -e "Sys.setenv(TAR = '/bin/tar'); devtools::install_github('stan-dev/rstantools')"
Rscript -e "Sys.setenv(TAR = '/bin/tar'); devtools::install_github('davidaknowles/leafcutter/leafcutter')"
```

This assumes perl is available. Now you can export this as a module named leafcutter

2) JunctionSeq

```bash
module load R/3.6 qorts

Rscript -e "install.packages('devtools', repos='http://cran.us.r-project.org', dependencies=TRUE, INSTALL_opts = c('--no-lock'))"
Rscript -e "BiocManager::install('JunctionSeq', version = '3.8')"
```
And export this 

3) Majiq
```bash
module load python3
wget "https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2"
tar -xvf htslib-1.9.tar.bz2
cd htslib-1.9/
./configure && make && make prefix=~/htslib-1.9 install
cd
export HTSLIB_LIBRARY_DIR=~/htslib-1.9/lib
export HTSLIB_INCLUDE_DIR=~/htslib-1.9/include
python3 -m pip install -U wheel setuptools numpy GitPython Flask==1.0.2 waitress==1.1.0 scipy>=1.1.0 psutil>=5.4.8 h5py>=2.8.0 gunicorn==19.9.0 Flask-WTF==0.14.2 Werkzeug==0.16.0
python3 -m pip install git+https://bitbucket.org/biociphers/majiq_stable.git#egg=majiq
```


## Instruction for the configuration file
TODO

### Citation
TODO

## Acknowledgements

We acknowledge Tobias Jakobi for helping with infrastructure support. We thank the Dieterich lab at the Heidelberg 
University Hospital for helpful discussions. 

## Contribution

Feedback or issue reports are welcome. Please submit via [the GitHub issue tracker](https://github.com/dieterich-lab/Baltica/issues)
Please provide the following information while submitting issue reports:
- Software and operation system versions
- Command line call
- If possible, also provide sample data so we can try to reproduce the error