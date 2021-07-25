# Getting started

## Quick example:

If Baltica dependencies, baltica configuration and cluster configuration are available, use:
```bash
baltica <workflow> <config> --use-singularity --profile <cluster>
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
- **cluster**: [Snakemake cluster profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) 

Or, read below.

!!! warning
    Snakemake is under active development. Please [contact us](https://github.com/dieterich-lab/Baltica/issues) if you have any issues with this documentation.

## Software environment:

Baltica framework is based on:  
- A python command-line interface  
- [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflows  
- Docker containers used with [Singularity](https://sylabs.io/singularity/)  
- R scripts for processing, integrating, annotating, assigning biological features, and reporting  
- a Rmarkdown report  

We have developed it on the following computer environments:
<!--  cat /proc/version -->
- Linux version 4.19.0-16-amd64 Debian 4.19.181-1 (2021-03-19)
- gcc version 8.3.0 
- Python version 3.7.7
- Singularity version 3.7.3
- Snakemake version 6.4.1
- Git version 2.20.1

These versions should not matter because the workflows are run within Docker containers, as long **Snakemake version > 6** and a **recent Singularity version**.

Baltica depends on python3, Singulary, and Snakemake.
- [How to install Singularity](https://sylabs.io/guides/3.0/user-guide/installation.html)  
- [How to install Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)  

## Installation

```bash
git clone https://github.com/dieterich-lab/baltica
cd Baltica
pip install .
```
Will install Baltica and its python dependencies. You may want to create a [virtual environment](https://realpython.com/python-virtual-environments-a-primer/) before installing Baltica.  
All other requirements are resolved with singularity containers.

!!! note
    We plan to submit Baltica to the Python Package Index.

!!! warning
    majiq requires an Academic or Commercial license for use. Users are required to [obtain their licenses.](https://majiq.biociphers.org/app_download/).

## Executing Baltica
Use baltica `cli` for current help documentation:
```bash
baltica --help
```

Baltica executor takes a single optional argument `--verbose`, to detail it execution. Every other option is passed to Snakemake.

## Cluster profile

Snakemake supports distributed workflow execution in many different high-performance computer clusters, as detailed [here](https://snakemake.readthedocs.io/en/stable/executing/cluster.html?highlight=profile#cluster-execution). We recommend using [cluster profiles](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) and using it like: 

```bash
baltica <workflow> <config> --use-singularity --profile <cluster> 
```

## (Advanced) Baltica workflows directly from Snakemake

Baltica workflows can be used directly with Snakemake without installation. However, there is limited support for it.

## References

[^1]: If you use Baltica, also please [cite Snakemake](https://bioinformatics.oxfordjournals.org/content/28/19/2520)
[^2]: If you use majiq results, please [cite it]( https://elifesciences.org/articles/11752)
[^3]: If you use leafcutter results, please [cite it](https://www.nature.com/articles/s41588-017-0004-9)
[^4]: If you use rmats, please [cite it](https://www.pnas.org/content/111/51/E5593) 
[^5]: If you use junctionseq results, please [cite it](http://nar.oxfordjournals.org/content/early/2016/06/07/nar.gkw501.full)
[^6]: If you use the Baltica's analysis module, please also [cite Stringtie](https://www.nature.com/articles/nbt.3122)


\bibliography
