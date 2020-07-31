# Baltica
___

 <img align="right" src="https://gist.githubusercontent.com/tbrittoborges/3c86ffbaa62e671771f443c65cb04fdc/raw/7ae0ea4a76e8f5464139ef34164c67de7a297ce8/baltica_logo.png" height="300"> *Baltica*: integrated splice junction usage analysis

## Features
- Snakemake workflows for DJU: JunctionSeq, Majiq, and Leafcutter
- Process and integrate the results from the methods  
- Summarise AS class of differently spliced junctions
- Benchmark using a set of artificial transcript in RNA Sequencing experiment
- (WIP) Consequence of differently spliced junctions  

## Instructions for using Baltica
We tested the instruction bellow in a computer cluster running Debian Linux distribution (v8.11 Jesse) and using the Slurm job scheduler. Baltica uses Snakemake, which support many difference 
grid engines systems, for instructions on how to use it on different systems, please see 
[Setting up Baltica: long edition](docs/setup.md) or [contact us](https://github.com/dieterich-lab/Baltica/issues),

## Documentation

The documentation is hosted by GitHub pages at https://dieterich-lab.github.io/Baltica/.

### Setting up Baltica (once):

[Here](https://dieterich-lab.github.io/Baltica/setup.html)
     

## Contribution  

Feedback and contribution are welcome. Please submit via [the GitHub issue tracker](https://github.com/dieterich-lab/Baltica/issues)
Please provide the following information while submitting issue reports:
- Software and operation system versions
- Command-line call
- If possible, also provide sample data so we can try to reproduce the error  


## Citation
Thiago Britto-Borges, Volker Boehm, Niels H. Gehring and Christoph Dieterich (2020)__Baltica: integrated splice junction usage analysis__. 
Manuscript in preparation.