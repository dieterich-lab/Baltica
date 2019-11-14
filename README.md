# Baltica

One stop solution for differential splicing analysis.


## Features

- Document example usage of the different tools 

## Usage 

## Documentation

1) Installation  
1) Configuration
1) Executing 
1) Analysis
1) Contributing

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
Note: if you are having problems with devtools trying using `gtar` instead of `tar`
```bash
Rscript -e "Sys.setenv(TAR = '/bin/tar'); devtools::install_github('davidaknowles/leafcutter/leafcutter')"
```


### Install Majiq

Majiq installation requires an license that can be obtained here:
[Acamdeic download](https://majiq.biociphers.org/app_download/).   
Majiq have multiple licenses that should reviewed at the developers website.

```bash
conda create --name majiq python=3.6 pysam numpy cython --yes
conda activate majiq
conda install --yes h5py>=2.8.0 Flask==1.0.2 Flask-WTF==0.14.2 GitPython>=2.1.11 gunicorn==19.9.0 psutil>=5.4.8 h5py>=2.8.0 scipy>=1.1.0 
env_path=$(dirname $(dirname $(which python)))
python_ver=$(python --version 2>&1 | awk '{print substr($2,1,3)}')
export HTSLIB_INCLUDE_DIR=$env_path/lib/python$python_ver/site-packages/pysam/include/htslib/
export HTSLIB_LIBRARY_DIR=$env_path/lib/python$python_ver/site-packages/pysam/include/htslib/htslib/

pip install git+https://bitbucket.org/biociphers/majiq_stable.git#egg=majiq
```

For older versions of pip (python<3.6), you may have to use `--process-dependency-links` to install any package missing 
from the conda installation. 

### Install JunctionSeq
Leafcutter is distributed with an [Apache License 2.0](https://github.com/davidaknowles/leafcutter/blob/master/LICENSE)

```bash
conda create --name junctionseq qorts r-biocmanager -c bioconda -c conda-forge --yes  
conda activate junctionseq
Rscript -e "BiocManager::install('JunctionSeq', version = '3.8')"
```

### Instal Whippet

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


