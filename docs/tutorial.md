# Step by step tutorial with the sample data set

## Setup
Make sure you have Singularity and Snakemake up and running. We have experienced problems with `singularity pull` and docker images when using a temporary directory (TMPDIR) in the shared file system. Setting `TMPDIR=/tmp/` resolves this issue. 

## Installation
Follow [Installation guide](setup.md#installation)

## Configuration
Open `Baltica/data/config.yml` and replace `/beegfs/homes/tbrittoborges/` or `/home/tbrittoborges` to your desired path. The __path__ parameter specifies where the project will be located. There is no needs to change **baltica_path**, as Baltica resolves it.

## Execution
Execute Baltica with `baltica all /Baltica/config.yml --use-singularity`. Add  __--quiet__ to reduce the logging level or __--verbose__ to increase it. Use a cluster profile (__--profile__) to take advantage of the cluster scheduler. See more useful parameters for snakemake below.

You should expect a `results/` directory containing the most relevant files by the end of the run. The report and excel table, `results/baltica_report{project_title}.html` and `results/baltica_table{project_title}.xlsx`, are the most relevant files. Intermediate files are kept in directories for each method (named for the methods), and can be used or deleted.

## Important snakemake parameters
- __--cores__: only required if you are not using a cluster scheduler. Use `baltica ... --cores all` to specify the maximum number of cores available.  
- __--profile__: setup a cluster configuration [profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles). 
- __--dry-run__: only computes the DAG, does not execute rules. But execution order may not reflect the actual order of execution.  
- __--unlock__: unlock the directory if a previous run had problems.  
- __--list-untracked__: list files that are not tracked by workflows, useful if you want to clean up the directory to save disk space after a successful run.
- __--quiet__: less verbose output. 
- __--reason__, __--printshellcmds__, __--verbose__: verbose output. `Simple use baltica ... --verbose` to get maximum debug information (you may want to redirect it to a file).

See more detail at the Snakemake [docs](https://snakemake.readthedocs.io/en/stable/executing/cli.html#useful-command-line-arguments).

