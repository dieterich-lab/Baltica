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

1. Extracting intron from the alignments files:
Here reads with the cigar string of type nMnNnM, where n are positive integers, are extracted from the alignments.
These reads have only alignment matches in the overhang.
 
1. Intron clustering:
Here introns with at least `minclureads` (default: 30) reads and up to `maxintronlen` kb (default: 100000) are
 clustered. The clustering procedure iteratively remove introns supported by less than `mincluratio` reads within a 
 cluster.

*NOTE*: Earlier versions of LeafCutter had problem with chromossome names different  from human  (chromossome) names in this step.   
*NOTE*: By default, leafcutter clustering does not use strandness information. In Baltica we set the, over-riding 
this default parameter. 


1. Differential splicing analysis


#### Software dependencies

- python==2.7
- base==3.5
- samtools==1.9

LeafCutter also depends on a series of R packages. 

#### Parameters

The following parameters where implemented in Baltica

- `minclureads` (default: 30) 
- `maxintronlen` (default: 100000)
- `mincluratio` (default: 0.001)
- `checkchrom` (default: TRUE)
- `min_samples_per_group` (default: 2) 
- `min_samples_per_intron` (default: 2) 
- `fdr` (default: 0.05)

#### Output

The relevant output from leafcutter are `leafcutter_ds_cluster_significance.txt` and 
`leafcutter_ds_cluster_significance.txt`, which are computed for each comparison.

Column descripton:

- `leafcutter_ds_cluster_significance.txt`:
1. `cluster`: TODO check identifier on the format `{chromosome}:{intron_start}:{intron_end}`
1. `Status`: whether this was tested or not
1. `loglr`: log likelihood ratio between the null model and alternative 
1. `df`: degrees of freedom, equal to the number of introns in the cluster minus one (assuming two groups)
1. `p` unadjusted p-value dor the under the asymptotic Chi-squared distribution

- `leafcutter_ds_effect_sizes.txt`:
1. `intron`: intron identifier on the format `chromosome:intron_start:intron_end:cluster_id`
1. `es`: TODO check fitted log effect size
1. `{cond_1}`: fitted junction usage in condition `cond_1`
1. `{cond_2}`: fitted junction usage in condition `cond_2`
1. `DeltapPSI`: difference between usage in the two conditions  

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
