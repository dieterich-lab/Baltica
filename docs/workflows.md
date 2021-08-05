# Workflow implementation

This chapter details the implementation and usage of each workflow in Baltica.

Baltica comprises a collection of Snakemake workflows (SMK files). Each file determines a series of sub-tasks (rules). The sub-tasks run in a specific order; once the output of every rule is complete, the workflow is considered successful. We implemented the workflows following instructions and parameters suggested by the methods authors unless otherwise noted.

![baltica_overview](img/Baltica_overview.png){ .center width=80% }
__Fig. 1 - Baltica overview__: Baltica is a framework to execute and integrate differential junction usage (DJU) analysis and further investigation enabled by the data integration.
__1 -- Input__: Baltica takes as input RNA-seq alignments, reference annotation, and a configuration file.
__2 -- Quality control__: As the first step of the pipeline, Baltica performs quality control of alignments with `RSeQC` and `FastQC`, which is reported by `MultiQC`.
__3 -- DJU and transcriptome assembly__: 
Next, Baltica computes DJU with `JunctionSeq`, `Majiq`, and `Leafcutter`, and uses `Stringtie2` to detected new transcripts and exons in the dataset.
__4 -- Downstream analysis__: 
Finally, we integrate the results from the DJU method. Optionally, Baltica can include an extra piece of evidence for DJU, such as DJU obtained from third-generation sequencing, to the integrated table. The set of introns is re-annotated using information from _de novo_ transcriptome annotation, and splice types between SJ and exons are assigned. Finally, a Baltica compiles a report with the most relevant information.

## Quality control workflow

Executed with:
```bash
baltica qc <config> --use-singularity
```

The first workflow comprises the quality control of the read alignments. 
This step aims to determine the success of sequencing and alignment.
Baltica includes workflows for RSeQC [@Wang2012] and FastQC [@andrews2012]. MultiQC [@Ewels_2016] summarizes the output from both tools.
In addition, users can use parameters from QC workflow into Baltica, such as maximum read length and library type.

Beyond the quality control, this step may help to identify differences among the RNA libraries.
For example, RSeQC provides the proportion of reads per feature in the input annotation. 
Differences between case vs. control, such as enrichment of reads aligned to introns, may suggest technical artifacts or global changes in splicing.
In addition, RSeQC provides the `junction_saturation.py` method, which quantifies the abundance of known and novel SJ in the RNA-seq alignments, and diagnoses if the alignment coverage detects known and novel splice junctions in sub-samples for alignments. 
Thus, users can use this functionality to identify the saturation of annotated and unannotated SJ, or a higher level of coverage is needed.
In conclusion, the quality control step serves to identify potential problems with the RNA-Seq library alignment and, potentially, direct on further troubleshooting and downstream analysis.

[Software dependencies](https://github.com/dieterich-lab/Baltica/blob/master/envs/qc.yml)

## RMATs workflow

Executed with:
```bash
baltica rmats <config> --use-singularity
```

RMATs [@Shen_2014] workflow is done in two steps:
    - Determine the experimental groups.  
    - Run `rmats.py`.  


Running RMATs [prep and post tasks separately](https://github.com/Xinglab/rmats-turbo/tree/8a2ad659717a1ccd6dbecd593dc1370ba7c30621#running-prep-and-post-separately) and paired statistical test were not implemented in Baltica. 

[Software dependencies](https://github.com/dieterich-lab/Baltica/blob/master/envs/rmats.yml)


## JunctionSeq workflow

Executed with:
```bash
baltica junctionseq <config> --use-singularity
```

JunctionSeq [@Hartley2016] workflows starts by junction read counts extraction done with QoRTs [@Hartley_2015]. 
In Baltica implementation for JunctionSeq workflow, we only consider reads that span multiple exons (splice junction reads, SJ) for annotated and unannotated introns, ignoring exon counts.
JunctionSeq uses disjoint genomic bins to flatten the transcriptome annotation.
To test the hypothesis that features are differently expressed in experimental groups, JunctionSeq fits a generalized linear model, as described in DEXSeq [@Anders2012], but reporting a test statistic at the genomic feature and gene level.
Unlike other DJU methods, JunctionSeq does not group the introns or S in AS events, so it does not compute PSI events but rather log fold change.

Baltica parses the `*_sigGenes.results.txt.gz` (located at `junctionseq/analysis/`) and discard entries that were flagged as not testable.

[Software dependencies](https://github.com/dieterich-lab/Baltica/blob/master/envs/junctionseq.yml) and [docker image recipe](https://github.com/dieterich-lab/Baltica/blob/master/docker/junctionseq/1.16.0/Dockerfile).

## Majiq workflow

Executed with:
```bash
baltica majiq <config> --use-singularity
```

Majiq workflow includes the following steps:
1. Create a configuration file (`majiq/build.ini`)
1. Converts the reference annotation from gtf to gff with `gtf2gff3.pl`
1. __majiq build__ generates the Splice Graph database with exons and SJ from the RNA-Seq experiment and the reference annotation
1. __majiq deltapsi__: computes &#968; and &#916;&#968; and tests if the &#916;&#968; significantly changes between comparisons. Introns are called significant if the probability of &#916;&#968; > `threshold` is higher than `non-changing-threshold`, where `threshold` and `non-changing-threshold` are the `--threshold` and `--non-changing-threshold` parameters, respectively 
1. __voila tsv__: filter and outputs the Majiq result to a tab-separated value file

Majiq visualization methods, such as `voila view`, are not currently implemented in Baltica but can be used independently. 

Baltica parses the `{comparison}_voila.tsv` files - one per comparison, located at `majiq/voila/`. 

[Docker image recipe](https://github.com/dieterich-lab/Baltica/blob/master/docker/majiq/2.2/Dockerfile).

## Leafcutter workflow

Executed with:
```bash
baltica leafcutter <config> --use-singularity
```

Leafcutter uses Regtools[@Cotto_2018] to extract SJ reads from the BAM files. Next, introns with at least `minclureads` reads clustered.
The clustering procedure iteratively discards introns supported by less than 
 `mincluratio` reads within a cluster. 
Finally, Leafcutter fits a Dirichlet-Multinominal model, which determines the SJ usage for each cluster  the usage (proportion) of a 
giving SJ within a cluster and compare this usage among conditions

The relevant output files from Leafcutter have the `_cluster_significance.txt` and `_effect_sizes.txt` suffix, computed for each comparison.

Column description:

`*_cluster_significance.txt`:
1. `cluster`: `{chromosome}:{intron_start}:{intron_end}`
1. `Status`: is this cluster testable?
1. `loglr`: the log-likelihood ratio between the null model and alternative 
1. `df`: degrees of freedom, equal to the number of introns in the cluster minus one (assuming two groups)
1. `p` unadjusted p-value dor the under the asymptotic Chi-squared distribution

`*_effect_sizes.txt`:
1. `intron`: intron identifier on the format `chromosome:intron_start:intron_end:cluster_id`
1. `es`: effect size
1. `{cond_1}`: fitted junction usage in condition `cond_1`
1. `{cond_2}`: fitted junction usage in condition `cond_2`
1. `deltapsi`: difference between usage in the two conditions  

[Software dependencies](https://github.com/dieterich-lab/Baltica/blob/master/envs/leafcutter.yml) and [Docker image recipe](https://github.com/dieterich-lab/Baltica/blob/master/docker/leafcutter/0.2.7/Dockerfile).

## Stringtie workflow

Baltica uses the gene, transcript, and class code assignments from the Stringtie output to the SJ from the DJU method outputs. In addition, exons defined by this annotation are used for assignments of splicing event types.
 
We process _de novo_ transcriptomic workflow with Stringtie [@pertea_2015]. First, we merge the alignment files from biological replicates. Next, we compute _de novo_ annotation with Stringtie, using the parameters  `-c 3 -j 3 -f 0.01`. Finally, the merge the multiple annotation with `gffcompare -r {reference_annotation.gtf} -R -V`. Details on the parameter selection are in the [Integration chapter](integration.md). 

[Docker image recipe](https://github.com/dieterich-lab/Baltica/blob/master/docker/stringtie/2.1.5/Dockerfile)

## References
\bibliography
