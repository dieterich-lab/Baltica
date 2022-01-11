## Change log

### v1.1.2 <small> Unreleased </small>

* Add support to unstranded RNA-seq data
* Add scripts for benchmark
* Add a new configuration `bind_paths` that allow integrating bam files from different projects

### v1.1.1 <small> July 23, 2021 (released in September 7 2021) </small>

* Add rmats workflow
* Add scrips for parsing for rmats and updated analysis to support the method
* Create the benchmark with the ONT Nanopore-seq
* Update benchmaks, included difference comparison for SIRV benchmark
* Split annotation and AS type assigment functions
* Update baltica table algorithm
* Add support for singularity container via snakemake, with container recipes `baltica qc config.yaml --use-singularity`
* Add parsing method for gffcompare tracking output
* Update configuration file to expose important parameters from the DJU methods
* Add end-to-end analysis with `baltica all config`  
* Experiment with meta-score (gradient boosted trees)
* Add baltica report and improved on report summaries
* Add orthogonal dataset use-case, to integrate third generation sequencing to the baltica table
* Change strand parameter to "fr-firststrand": "reverse", "fr-secondstrand": "forward" or unstranded, fix error in rmats strand

### v1.0 <small> September 17, 2020</small>

* Add `is_novel` column, indication introns not into the reference annotation
* Remove unitended columns (X1, ...) from the report

### v1.0 <small>- July 23, 2020</small>

* First public release comprises of DJU methods Leafcutter, Junctionseq and Majiq. Stringtie for *de novo* transcriptomics assembly. FastQC and MultiQC (#1).
