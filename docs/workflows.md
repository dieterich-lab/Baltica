# Workflows

This document details on the implementation and usage for each workflow in Baltica.

Baltica comprises a collection of Snakemake workflows (snakemake files with extension .smk). Each file determines a series of sub-tasks (rules). Generally, the sub-tasks run in a specific order and only successful if their output exits once the execution finishes. 

1. Quality control   
1. Differential Junction Usage (DJU)   
1. _de novo_ transcriptomics with Stringtie   
1. Analysis  

Which are achieved by successively calling: 

```{bash}
baltica qc config.yml  
baltica junctionseq config.yml  
baltica majiq config.yml  
baltica leafcutter config.yml  
baltica stringtie config.yml  
baltica analysis config.yml  
```

Please the [Tutorial](tutorial.md) for a step-by-step guide on who to runs these analysis on sample data. 

Baltica requires RNA-Seq read alignments as input. At this time, the impact of the read aligment in DJU methods results is unknown, but we suggest [STAR](https://github.com/alexdobin/STAR) for reads alignment. 

!!! important
    The transcriptome annotation is also important, so make sure you use the same annotation for read aligment and in baltica parameter.
   

## Baltica configuration

Baltica configuration files are JSON file used by Snakemake [configuration file](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html) 
in the `json` or `yaml` formats. Please see a minimal working example [here](https://github.com/dieterich-lab/Baltica/blob/master/baltica/config.yml). The following parameters are mandatory:

Parameter name | Description | Note     
-------------- | ----------- | ---- 
`path` |  absolute path to the logs and results  |  
`sample_path` | path to the parent directory to the alignment files |  
`samples` | sample name to alignment files (BAM format) a condition name | [^1]
`comparison` |  pair of groups to be tested |
`ref` |  full path to the reference transcriptome annotation in the (GTF format) | 
`ref_fa` | full path to the reference genome sequence in the FASTA format | [^2]  
`*_env` | Used if the required dependency is __not__ available in the path | [^3]
`strandedness` | choice between `reverse`, `forward` or None | [^4] 
`read_len` | positive integer representing the maximum read length | [^5]

[^1]: We use the following convention for the sample name: `{condition}_{replicate}`, where condition is the experimental group name without spaces or underline character, and replicate is a positive integer
[^2]: Use by Majiq for for GC content correction
[^3]: User should prefer `--use-conda` or `--use-envmodules`; however if it's not possible to load the requirements for, this hack may help 
[^4]: Check RSeQC infer_experiment result 
[^5]: Check FastQC Sequence Length Distribution report

<!-- TODO detail on this usecase -->

!!! note
    Junctionseq and Leafcutter support more complex designs but these are not currently implemented in Baltica

## Quality control workflow

This section describes the quality control steps used to check for problems in RNA-Sequencing data sets, focusing on
  exon-exon junctions reads. Here is an incomplete lite of parameters from RNA-Sequencing experiment that are 
  necessary for splicing analysis:
  
  - Sequencing depth
  - Number of technical replicate
  - Read length
  - Read overhang
  - Splice junction saturation
  - Annotation quality
  
For RNA-Sequencing experiments aiming to detect lowly-expression genes and transcripts, a higher sequencing depth (40-60 
 million reads) is required, in contrast, to experiment that only aim finding the most abundant genes, and so only demand around 10 million reads [see](https://support.illumina.com/bulletins/2017/04/considerations-for-rna-seq-read-length-and-coverage-.html). This parameter is particularly relevant for samples with novel splice junctions (SJ). Read length and paired-end are also critical
 for splice junction identification, and longer reads offer more coverage of the exons boundaries (
 see [Chhangawala et al., 2015](https://doi.org/10.1186/s13059-015-0697-y)). The target nominal read length 
 (without adaptor or barcodes) should be between 75-100 nucleotide, maximize the read overhang size, and, consequently, maximize the quality of spliced suggest a higher number of replicates versus deeper sequencing 

Also, databases such as the CHESS [2](http://ccb.jhu.edu/chess/) can provide additional evidence for splice sites absent in the annotation. The sequencing depth is related to the splice junction saturation metric. 
 This metric is defined by the RSeQC tool, which implements a sampling procedure to identify the percentage of annotated and novel introns are observed in sub-samples of the data.
 
For more details, please consult the documentation:
  - [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/) 
  - [RSeQC](http://rseqc.sourceforge.net/#usage-information)
  - [MultiQC](https://multiqc.info/docs/)

### Software dependencies

The workflow wes tested with the following software version:  


- RSeQC == 2.6.4
- FastQC == 0.11.8 
- MultiQC == 0.8 

## Differential splicing workflow

This section details the implementation of the DJU. In overview, the workflows have 
common steps:
1. Extracting split reads that support the SJ 
1. Defining the splicing events  
1. Modeling the SJ counts

### Leafcutter workflow

Leafcutter uses a series of scrips to extract the split reads from the BAM files. Recently, this step was changed 
 to use (regtools)[https://github.com/griffithlab/regtools] to speed up the process, but our preliminary test 
 show this change affects the workflow results, instead only speeding up the process.
1. Extracting intron from the alignments files: reads with M and N cigar are extracted from the alignments, giving a
minimum read overhang
1. Intron clustering: introns with at least `minclureads` (default: 30) reads and up to `maxintronlen` kb (
 default: 100000) are clustered. The clustering procedure iteratively discards introns supported by less than 
 `mincluratio` reads within a cluster.
1. Differential splicing analysis: Leafcutter uses a Dirichlet-Multinominal model to model the usage (proportion) of a 
giving SJ within a cluster and compare this usage among conditions

!!! info
  If you annotation use non-canonical chromosome names, you may have to change Leafcutter clustering scrip
!!! info
  By default, Leafcutter clustering does not use the read strandedness information. In Baltica, we override these parameters and use the strand information for clustering.

#### Software dependencies

- python == 2.7
- base R == 3.5
- samtools == 1.9

Leafcutter also depends on a series of R packages. 

#### Parameters

- `minclureads` (default: 30) 
- `maxintronlen` (default: 100000)
- `mincluratio` (default: 0.001)
- `checkchrom` (default: TRUE)
- `min_samples_per_group` (default: 2) 
- `min_samples_per_intron` (default: 2) 
- `fdr` (default: 0.05)

#### Output

The relevant output from Leafcutter is `_cluster_significance.txt` and 
`_effect_sizes.txt`, which are computed for each comparison. Baltica parses these files.

Column description:

- `*_cluster_significance.txt`:
1. `cluster`: TODO check identifier on the format `{chromosome}:{intron_start}:{intron_end}`
1. `Status`: is this cluster testable?
1. `loglr`: the log-likelihood ratio between the null model and alternative 
1. `df`: degrees of freedom, equal to the number of introns in the cluster minus one (assuming two groups)
1. `p` unadjusted p-value dor the under the asymptotic Chi-squared distribution

- `*_effect_sizes.txt`:
1. `intron`: intron identifier on the format `chromosome:intron_start:intron_end:cluster_id`
1. `es`: TODO check fitted log effect size
1. `{cond_1}`: fitted junction usage in condition `cond_1`
1. `{cond_2}`: fitted junction usage in condition `cond_2`
1. `DeltapPSI`: difference between usage in the two conditions  

## Majiq workflow

1. MAJIQ Builder - creates the Splice Graph database with exons and SJ from the RNA-Seq experiment
1. PSI analysis - compute PSI and deltaPSI

Majiq also provides a visualization with the `voila view`. 

### Software dependencies

python == 3.6
htslib == 1.9 
Among other python dependencies that are automatically installed.

### Parameters
- Used in Majiq INI file:
   `assembly:`
    - Description: name of the assembly on the UCSC genome browser
    
    `strandness:`
    - Description: RNA-Sequencing library type 
    - default: reverse 
    
    `read_len:`
    - Description: maximum read lenght to be considered
    - Default: 100

- For Majiq build:  

   `--min-experiments`
    - Description:Iinteger or proportion of the minimum number of experiments a LSV event is observed to be considered
  - Default: 0.5
   
   `--min-intronic-cov`
  - Description: Minimum number of coverage a intron needs to be tested for intron retention events
  - Default: 0.01
  
   `--min-denovo`
  - Description: Minimum number of reads at all positions in a LSV to consider a de novo junction 
  - Default: 2
  
    `--minreads`
    - Description: Minimum number of reads at all positions in an LSV to consider an LSV
  - Default: 3
  
  `--minpos`
  - Description: Minimum number of start positions with at least 1 read in a LSV to consider that the LSV "exist in the data".
  - Default: 2

    `--markstacks`
  - Description: p-value threshold (or negative number to disable) to mark stack positions 
  - Default: 1e-07

    `--k`
  - Description: Number of positions to sample per iteration.
  - Default: 50

    `--m`
  - Description: Number of sampling steps using on bootstrap
  - Default: 30
  - baltica:

- For Majiq deltapsi:  
  
    `--binsize`
  - Description: Number of bins for the PSI value distribution 
  - Default: 0.025

    `--prior-minreads`
    - Description: Minimum number of reads at all to included 
  - Default: 20

    `--prior-minnonzero`
  - Minimum number of positions for the best set.
  - Default: 10

For voila tsv:  

   `--threshold`  
    - Description: Discard LSVs if the probability is lower then this value
    - Default: 0.2

   `--non-changing-threshold`
    - Description: None 
  - Default: 0.05

   `--probability-threshold` 
  - Description: None
  - default: off



## JunctionSeq workflow

### Software dependencies

### Parameters

## Integration workflow

#### Software dependencies

#### Parameters


