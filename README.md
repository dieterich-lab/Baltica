# Baltica

One stop solution for differential splicing analysis.


## Features

- Document example usage of the different tools 

## Usage 


## Quickstart (LeafCutter example)
[Instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) for SnakeMake installlation

1) Clone and install Baltica:  
 ```bash
 git clone dieterich-lab/Baltica/  
 cd Baltica; pip install .
```
2) [Install LeafCutter](#install-leafcutter)
3) Within Baltica folder fill the template `config.yml` [Instructions](#Instruction_for_the_configuration_file)
4) run baltica with ` baltica ` 
## Installation guide



### Install Leafcutter
```bash
conda create --name leafcutter python=2.7 --yes
conda activate leafcutter
conda install -c bioconda samtools r-base=3.5 --yes

Rscript -e "install.packages('devtools', repos='http://cran.us.r-project.org')"
Rscript -e "devtools::install_github('davidaknowles/leafcutter/leafcutter')"
```

## Install Majiq

Majiq installation requires an license that can be obtained here:
[Acamdeic download](https://majiq.biociphers.org/app_download/).   
Majiq have multiple licenses that should reviewed at the developers website.

```bash
conda create --name majiq python=3.5 pysam numpy cython --yes
conda activate majiq
env_path=$(dirname $(dirname $(which python)))
export HTSLIB_INCLUDE_DIR=$env_path/envs/MAJIQ/lib/python3.6/site-packages/pysam/include/htslib/
export HTSLIB_LIBRARY_DIR=$env_path/envs/MAJIQ/lib/python3.6/site-packages/pysam/include/htslib/htslib/

pip install git+https://bitbucket.org/biociphers/majiq_stable.git#egg=majiq
```

## Install JunctionSeq
Leafcutter is distributed with an [Apache License 2.0](https://github.com/davidaknowles/leafcutter/blob/master/LICENSE)

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
```

## Instruction for the configuration file


### Citation
TODO

### Acknowledgements

### Contribute 
Contributions as feed-back and issues should be submitted at this project GitHub issue tracker TODO


