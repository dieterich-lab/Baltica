# Baltica

Reconciling alternative splicing identification

# Quick start

```bash
pip install Baltica
```

`Baltica leafcutter ` 

Baltica depends on python 3 (>=3.4) and snakemake, please see [Installation](#installation) for details

## Features

Based on the results from different tools, Baltica can:

- Document example usage, pros and cons  
- Process and join the results  
- Produce reports  
- Use data integration methods to detect consequences of AS  

# Introduction 

The increasing volume and quality in RNA-Sequencing data has reveled the importance of correct splicing to human health.
Many computational methods to identify differential splicing from RNA-Seq experiments have been developed in the past 10
 years. These methods can be classified in three categories:
- Methods that compute exon usage  
- Methods that estimate transcript abundance  
- Methods that identify events from exon-exon junction reads   

To understand differential splicing, due to genetic variation or mis-regulation, is critical to understand in detail
 the molecular mechanism of disease, thus important for new . We take a pragmatic approach to this problem working with
 the state of the art tools for tools that identify events from exon-exon junction reads. This tools have a the 
 advantage of identifying un-annotated exon-exon junctions, which pivotal for . Nonetheless, each tool has it's own pros and cons [TODO LINK]
 and we suggest users to analyse the results from 2 or more tools.

# Installation

0) optional create a Baltica environment:

```bash
conda create --name Baltica python=3.6 snakemake PyYAML
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
[Academic download](https://majiq.biociphers.org/app_download/).   
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

### Install Whippet

```bash
conda create -n whippet -c conda-forge julia=0.6 --yes 
conda activate whippet
julia -e 'Pkg.add("Whippet"); using Whippet'
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