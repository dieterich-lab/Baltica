# Baltica

One stop solution for differential splicing analysis.


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