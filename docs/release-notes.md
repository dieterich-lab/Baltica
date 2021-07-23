## Change log

### Master July 23, 2021 (unreleased)
* Add rMATs workflow
* Add scrips for parsing for rMATS and updated analysis to support the method
* Create the benchmark with the ONT Nanopore-seq
* Update benchmaks, included difference comparison for SIRV benchmark
* Splite annotation and AS type assigment functions
* Update baltica table algorithm 
* Add support for singularity container via snakemake, with container recipes `baltica qc config.yaml --use-singularity`
* Add parsing method for gffcompare tracking output
* Update configuration file to expose important parameters from the DJU methods
* Add end-to-end analysis with `baltica all config`  
* Experiment with meta-score (gradient boosted trees) 
* Add baltica report and improved on report summaries 
* Add orthogonal dataset use-case, to integrate third generation sequencing to the baltica table
* Change strand parameter to "fr-firststrand": "reverse", "fr-secondstrand": "forward" or unstranded, fix error in rMATS strand 

### 1.1 <small> September 17, 2020</small>
* Add `is_novel` column, indication introns not into the reference annotation
* Remove unitended columns (X1, ...) from merge

### 1.0 <small>- July 23, 2020</small>

* First public release comprises of DJU methods Leafcutter, Junctionseq and Majiq. Stringtie for *de novo* transcriptomics assembly. FastQC and MultiQC (#1).
