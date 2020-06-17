# Setting up baltica

- [Install miniconda](*https://docs.conda.io/en/latest/miniconda.html)
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
``` 
Follow the instructions to finish your installation, which by default is at $HOME/miniconda3

Make sure you initialize conda with `conda init`
You can test whether your installation was successful or not by running `conda --version` and you may need to restart 
your shell instance. 

- Clone baltica
```bash
git clone git@github.com:dieterich-lab/baltica.git
```

Baltica can be used with the [modules system](https://modules.readthedocs.io/en/latest/index.html) or conda
environments. Here we describe the installation with conda. 

- Install snakemake 
	```bash
    conda install -c bioconda snakemake==5.2 --yes
	```

Many of dependencies can created with environments with conda, go to `baltica/envs` directory and install with:
    - `conda env create -f stringtie.yml --yes` 
    - `conda env create -f qc.yml --yes`


In general, R packages do not play nicely with conda, but we still use it because it's flexibility and the ability to create isolated software environments.


## Majiq installation

IMPORTANT: Majiq requires licenses for Academic or Commercial use. Users are required to obtain their license. 
[Academic download](https://majiq.biociphers.org/app_download/).   

Majiq can installation can be problematic, but the recipe bellow works for us:

```bash
conda create --name majiq_env python=3.6 pysam numpy cython --yes -c bioconda
conda activate majiq_env
conda install --yes waitress==1.1.0 h5py>=2.8.0 Flask==1.0.2 Flask-WTF==0.14.2 GitPython>=2.1.11 gunicorn==19.9.0 psutil>=5.4.8 h5py>=2.8.0 scipy>=1.1.0
conda install --yes -c bioconda htslib 

python3 -m pip install --upgrade pip
python3 -m pip install git+https://bitbucket.org/biociphers/majiq_stable.git#egg=majiq
```

## LeafCutter installation

Users can install Leafcutter with conda using the following recipe: 

```bash
conda create --name leafcutter python=2.7 --yes
conda activate leafcutter
conda install -c bioconda samtools r-base=3.6 --yes

Rscript -e "install.packages('devtools', repos='http://cran.us.r-project.org', dependencies=TRUE, INSTALL_opts = c('--no-lock'))"
Rscript -e "Sys.setenv(TAR = '/bin/tar'); devtools::install_github('stan-dev/rstantools')"
Rscript -e "Sys.setenv(TAR = '/bin/tar'); devtools::install_github('davidaknowles/leafcutter/leafcutter')"
```

*Note*: if you are having problems with devtools trying using `gtar` instead of `tar` use the following:

```{bash eval=FALSE, echo=TRUE}
Rscript -e "Sys.setenv(TAR = '/bin/tar'); devtools::install_github('davidaknowles/leafcutter/leafcutter')"
```
## JunctionSeq installation

JunctionSeq should be installed directly from BioConductor

```bash
conda env create -f junctionseq-env.yml --yes
conda activate leafcutter
Rscript -e "BiocManager::install('JunctionSeq',  INSTALL_opts = c('--no-lock'))"
```
