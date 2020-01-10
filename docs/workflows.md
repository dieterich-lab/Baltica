This document details on the implementation and usage for each workflow in Baltica.

Baltica comprises a series of SnakeMake workflows (snakemake files with extension .smk). Each file comprise a series 
of sub-tasks (rules). Generally, the sub-tasks should be executed in a specific order. 


`QC (optional) > Differential Splicing > Integration`

Baltica requires reads to be aligned to the transcriptome, and for that we tend to use STAR []() and Ensembl 
transcriptome annotation. 

## Configuration file: common parameters

Baltica configuration files are SnakeMake [configuration file](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html) 
in the `json` or `yaml` formats.

- `sample_path`: parent directory to the alignment files
 
- `workdir`: parent directory where you want your results. This parameters is used when a project have multiple batches 
of alignments that are analysed independently  
 
- `samples`: mapping for the alignment files (.bam format) in the format `sample_name: sample_path`
We use the following convention for the sample name: `{condition}_{replicate}`, where condition is condition name, 
without space or special characters, and replicate is a positive integer. 

- `contrast`: which comparisons should be make? Currently JunctionSeq and LeafCutter support more complex designs, 
than the pairwise comparison, but this is currently unsupported by Baltica.

- `ref`: reference transcriptome annotation in the GTF format

- `ref_fa`: (optional) reference genome sequence in the fasta format. Only used for the `de novo transcriptomics` rules
TODO link.

- `*_env`: (optional) variable to activate environments, such as modules, conda environments or pyenvs. Only used if 
the target software is *not* in path.

## Quality control workflow

This session describes the quality control steps used to check for problems in RNA-Sequencing data sets, focusing on
  exon-exon junctions reads. Here is a incomplete lite of parameters from RNA-Sequencing experiment that are 
  important for splicing analysis:
  
  - Sequencing depth:
  - Number of technical replicate:
  - Read length:
  - Read overhang:
  - Splice junction saturation:
  - Annotation quality:

Here I summarise suggestions for the ideal these factors. This are general suggestions. 

In RNA-Sequencing experiments aiming to detect low-expression genes and transcripts, a higher sequencing depth (40-60 
million reads) is required, in contrast to experiment that only aim finding the most abundant genes, and so only 
demand around 10 million reads (1)[https://support.illumina.com/bulletins/2017/04/considerations-for-rna-seq-read-length-and-coverage-.html] . If the sequencing depth 
This is particularly relevant for tissue or organisms that are poorly annotated. In addition, databases such as the
CHESS (2)[http://ccb.jhu.edu/chess/] can provide additional evidence for splice sites absent in the annotation. 
The sequencing depth is related to the splice junction saturation metric. This metric is defined by RSeQC tool, which
implements a sampling procedure to identify the percentage of annotated and novel introns are observed in sub-samples 
of the data.

The target nominal read length (without adaptor or barcodes) should be between 75-100 nucleotide ()[ref], to maximise
the read overhang size and, consequently, maximise the quality of the alignments of spliced reads. 

### Software dependencies

The workflow wes tested with the following software version:  

- RSeQC == 2.6.4
- FastQC == 0.11.8 
- MultiQC == 0.8 

## Differential splicing workflow

This section details on the implementation of each differential splicing workflow 

### LeafCutter workflow

LeafCutter implements . This workflow can be explained by three tasks:

1. Extracting intron from the alignments:
Here   
1. Intron clustering 
1. Differential splicing analysis


#### Software dependencies

#### Parameters

## Majiq workflow

### Software dependencies

### Parameters

- `threshold:` --threshold 0.1
- `assembly:` GRCh38_96
- `strandness:` reverse
- `read_len:` 100


## JunctionSeq worlflow

### Software dependencies

### Parameters

## Integration workflow

#### Software dependencies

#### Parameters
