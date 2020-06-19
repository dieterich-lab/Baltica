# Baltica

*Baltica*: integrated differential junction usage (DJU) and consequence analysis with an ensemble of methods

## Features

- SnakeMake workflows for DJU: JunctionSeq, Majiq and LeafCutter
- Process and integrate the results from the methods  
- Summarise AS class of differently spliced junctions
- Benchmark using a set of artificial transcript in RNA-Sequencing experiment
- (WIP) Consequence of differently spliced junctions

## Instructions for using Baltica

This instructions are specific for computer cluster using Slurm. Baltica uses Snakemake, which support many difference 
grid engines systems, for instructions on how to use it on different systems please see 
[Setting up Baltica: long edition](docs/other_file.md) or contact us. 

### Setting up Baltica (once):
- Clone Baltica
	```bash
	git clone git@github.com:dieterich-lab/baltica.git
	```

- (optional) Setup snakemake cluster profile
	```bash
	module load python
	pip install cookiecutter --user
	mkdir -p ~/.config/snakemake
	cd ~/.config/snakemake
	cookiecutter https://github.com/Snakemake-Profiles/slurm.git
	``` 
    Make sure your $PATH is append by the locally installed python packages:
    ```bash
    export PATH="$HOME/.local/bin/:$PATH"
    ```
    
 - (optional) Test with the majiq workflow with a toy dataset. Should take less than a minute:
    ```bash
    snakemake -s rules/majiq.smk --configfile config.yml --cores 10 --use-envmodule
    ```
  
   
### For a specific project (for each project)
- (optional) Clone baltica to keep changes under version control
- set the configuration file
- run baltica
	- (optional) quality control 
	- DJU workflows
	- de novo transcripts 
	- Analysis 

If Baltica or snakemake is not in your $PATH use

The cluster-profile allow us to centralize the cluster configuration for snakemake workflows. 
Snakemake can be run with `snakemake --profile slurm`.
For other cluster configuration templates see [Snakemake profiles](https://github.com/Snakemake-Profiles/)
## Documentation
Please see documentation:
   - [Setting up Baltica: long edition](docs/setup.md)
   - [Workflows](docs/worflows.md)  
   - [Analysis](docs/analysis.md)
   - [Benchmark](docs/benchmark.md)
 
## Contribution

Feedback or issue reports are welcome. Please submit via [the GitHub issue tracker](https://github.com/dieterich-lab/Baltica/issues)
Please provide the following information while submitting issue reports:
- Software and operation system versions
- Command line call
- If possible, also provide sample data so we can try to reproduce the error

## Citation
TBD