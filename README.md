# Baltica
___

 <img align="right" src="https://gist.githubusercontent.com/tbrittoborges/3c86ffbaa62e671771f443c65cb04fdc/raw/7ae0ea4a76e8f5464139ef34164c67de7a297ce8/baltica_logo.png" height="300"> *Baltica*: integrated differential junction usage (DJU) and consequence analysis with an ensemble of methods

## Features
- SnakeMake workflows for DJU: JunctionSeq, Majiq and LeafCutter
- Process and integrate the results from the methods  
- Summarise AS class of differently spliced junctions
- Benchmark using a set of artificial transcript in RNA-Sequencing experiment
- (WIP) Consequence of differently spliced junctions  

## Instructions for using Baltica
This instructions are specific for computer cluster using Slurm. Baltica uses Snakemake, which support many difference 
grid engines systems, for instructions on how to use it on different systems please see 
[Setting up Baltica: long edition](docs/setup.md) or contact us.

### Setting up Baltica (once):
- Clone Baltica
	```bash
	git clone git@github.com:dieterich-lab/baltica.git
	```
- (_optional_) Setup snakemake cluster profile
The cluster-profile allow us to centralize the cluster configuration for snakemake workflows. 
Snakemake can be run with `snakemake --profile slurm`.
For other cluster configuration templates see [Snakemake profiles](https://github.com/Snakemake-Profiles/)
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
 - Install baltica with:
    ```bash 
    cd Baltica/
    python setup.py install
    ```
 - Test with the toy dataset:
    ```bash 
    baltica majiq --configfile baltica/config.yml --use-envmodule
    ```
 - (optional) Use baltica without installing in the toy dataset :
    ```bash
    snakemake -s baltica/majiq.smk --configfile config.yml --cores 10 --use-envmodule
    ```

### For a specific project (for each project)
- (_optional_) Clone baltica to keep changes under version control
- set the configuration file
- run baltica
	- (optional) quality control 
	- DJU workflows
	- de novo transcripts 
	- Analysis 
     
## Documentation
Please see documentation:
   - [Setting up Baltica: long edition](docs/setup.md)
   - [Workflows](docs/workflows.md)  
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
