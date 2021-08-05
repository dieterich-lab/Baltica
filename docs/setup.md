# Getting started

## Quick example:

If Baltica dependencies, baltica configuration and cluster configuration are available, use:
```bash
baltica <workflow> <config> --use-singularity
```
- **workflow**:  
    - all: run end-to-end wokflows
    - qc: run quality control 
    - stringtie: run *de novo* and guided transcriptome assembly
    - rmats
    - junctionseq
    - majiq
    - leafcutter 
    - analysis: run scripts for integration, annotation and reporting
- **config**: [project configuration file](proj-config.md)

Use `--profile <cluster>` with a [Snakemake cluster profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) 
or set the number of avaiable cores with `--cores`. For more Snakemake parameters, [check their documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html).

!!! warning
    Baltica is under active development. Please [contact us](https://github.com/dieterich-lab/Baltica/issues) if you have any issues with this documentation.

## Software environment:
Baltica framework is based on:  
- A python command-line interface  
- [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflows  
- Docker containers used with [Singularity](https://sylabs.io/singularity/)  
- R scripts for processing, integrating, annotating, assigning biological features, and reporting  
- a Rmarkdown report  

We have developed it on the following computer environments:
- Linux version 4.19.0-16-amd64 Debian 4.19.181-1 (2021-03-19)  
- gcc version 8.3.0  
- Python version 3.7.7  
- Singularity version 3.7.3  
- Snakemake version 6.4.1  
- Git version 2.20.1  

These versions should not matter because the workflows are run within Docker containers, as long **Snakemake version > 6** and a **recent Singularity version**.

Baltica depends on python3, Singulary, and Snakemake:  
- [How to install Singularity](https://sylabs.io/guides/3.0/user-guide/installation.html)  
- [How to install Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)  

## Installation
```bash
git clone https://github.com/dieterich-lab/baltica
cd baltica
pip install .
```
Will install Baltica and its python dependencies. You may want to create a [virtual environment](https://realpython.com/python-virtual-environments-a-primer/) before installing Baltica.  
All other requirements are resolved with singularity containers. Baltica store its singularity containers at `$HOME/.baltica/singularity/`.

!!! note
    We plan to submit Baltica to the Python Package Index.

!!! warning
    majiq requires an Academic or Commercial license for use. Users are required to [obtain their licenses.](https://majiq.biociphers.org/app_download/).

## Executing Baltica
Use baltica `cli` for current help documentation:
```bash
baltica --help
```

Baltica executor takes a single optional argument `--verbose`, to detail its execution. Every other option is passed to Snakemake.

## Test dataset
Baltica ships with a test dataset, located at the `data/` directory. There is a configuration file for the test dataset. **Users are required to update this configuration file**. Please see the [Baltica project configuration](proj-config.md) for further details.

## Cluster profile
Snakemake supports distributed workflow execution in many different high-performance computer clusters, as detailed [here](https://snakemake.readthedocs.io/en/stable/executing/cluster.html?highlight=profile#cluster-execution). We recommend using [cluster profiles](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) and using it like: 

```bash
baltica <workflow> <config> --use-singularity --profile <cluster> 
```

## (Advanced) Baltica workflows directly from Snakemake
Baltica workflows can be used directly with Snakemake without installation. However, there is limited support for it.
