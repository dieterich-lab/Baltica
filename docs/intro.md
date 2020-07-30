# Introduction 
## The life of RNA transcripts is complex 
Alternative promoter usage, Alternative splicing (AS), and alternative polyadenylation site (PAS) usage are processes that contribute to the transcriptome complexity by producing different RNA transcript isoforms.
AS is defined by a series of enzymatic reactions by which ribonucleoprotein complexes (spliceosomes) sequentially excise introns from a premature mRNA (pre-mRNA) and ligate the donor exon, at the 5' end of the intron, to the acceptor exon, at the 3' end of the intron. 
Once complete, the series of splicing events produces an mRNA isoform, which encodes the protein-coding sequence. 
Alternative combinations of exons produce different mRNA isoforms; thus, AS enables one pre-mRNA to generated many transcripts and protein isoforms, consequently diversifying the gene function. 
In addition, the combinatorial usage of exons augments the possibilities of differential regulation on the transcript level.

## Elements of splicing regulation

Cis- and trans-acting elements regulate RNA Splicing by orchestrating processes such as the exon-intron boundaries definition.
The cis-acting elements are local features, such as the primary and secondary mRNA structure.
Trans-acting elements are RNA binding proteins (RBPs) that, once deposited on the mRNA, regulate many steps of splicing.
AS and its regulation are essential for many physiological processes, such as development [@Wong2013] and tissue remodeling [@Liu2018].
Moreover, defective pre-mRNA splicing, or mis-splicing, has been extensively linked to human disease, as reviewed by [@Scotti2015], and is often associated with genetic variation.
A recent study has demonstrated that up to 10% of human genetic variants with causal links to neurodevelopmental disorders are predicted to cause mis-splicing [@Jaganathan2019].

## Consequence of splicing and motivation
Changes in pre-mRNA splicing can have drastic consequences for protein function.
Spliced isoforms may lead to truncated or extended protein domains, with a significant impact on protein function.
Also, AS changes have been implicated in changes to the subcellular trafficking of proteins [@Link_2016].
The encoding of a premature stop codon due to AS activates the nonsense-mediated decay pathway, which leads to the degradation of this transcript and, thus, depletion of the encoded protein [@Lewis_2002,@Baserga_1992].
However, experimental evidence to support AS biological consequence is limited in the scientific literature because of the many challenges in detecting and prioritizing AS events.
For example, VastDB, a database of experimentally validated splicing events, contains 1148 events (version 1.8). [@Tapial_2017].  
<!-- \todo{There are examples of software that touch the topic of protein function changes in AS, notably  IsoformSwitchAnalyser, apprisws and I need to mention those somewhere}  -->

## Classes of computational methods for AS identification
The number of reported AS events has increased enormously over the past ten years.
This increase is primarily due to advances in the RNA-Sequencing (RNA-Seq) technology, notably the longer read-length, and the development of computational methods to detect AS from RNA-Seq.
A read that aligns with two or more exons is the evidence that a intron was removed from the pre-mRNA.
These reads are named exon-exon junctions or splice junctions (SJ) and represented by the gaps marked with the N cigar in the read alignment.
The splice graph is a network representation in which nodes and edges represent exon and SJ, respectively, and the different network paths form the different transcripts.
Methods to detect AS from RNA-Seq build a model based on the feature abundances and its difference between two or more experimental conditions to test for AS differences.
These methods fit into three classes: (a) different transcript usage or expression (DTU) [@Trapnell_2010]<sup>,</sup>[@Nowicka2016]<sup>,</sup>[@Froussios2019]; (b) different exon usage [@Anders2012]; and (c) different junction usage (DJU) [@Li2017]<sup>,</sup>[@Hartley2016]<sup>,</sup>@VaqueroGarcia2016]<sup>,</sup>[@Trincado2018]<sup>,</sup>@SterneWeiler2018].
The genomic feature used for statistical modeling determines the class of each method.  
Most DTU methods depend on a known transcriptome as input, and the results they report rely heavily on the quality of this annotation [@Soneson_2016].

## DEXSeq applied to differential exon usage
DEXSeq is a popular method for DEU that is currently maintained.  
To overcome the complexity of the splice graph representation, DEXSeq uses a flat transcript structure representation.
It split exons in genomic bins, and partly overlapping exons with alternative a start or end coordinates are split in distinct bins. 
The bin RNA-Seq read coverage is then modeled with a generalized linear model, which enables the comparison of experimental conditions in context to its covariates.
This approach, however, also requires the exon coordinates as input. 

## DJU methods and their advantages
In contrast, DJU methods are less dependent on the annotation and are more robust to ambiguous sequencing that reads that map to many features. 
The main advantage of DJU versus DTU methods is that PSI, the proportion of reads supporting a path in the splice graph, can be computed directly from the sequencing reads, opposing the transcript level quantification, which depends on the complexity of the splice graph.  
Thus, DJU methods enable the _de novo_ identification and quantification of SJ from the RNA-Seq data independent of transcriptome assembly methods.
This feature enables the identification of many more splicing events that otherwise would be missed.
There are considerable differences in terms of implementation and assumption among DJU methods. Nevertheless, these methods share a set common set of steps that represent the core of a DJU workflow.
First, the SJ are extracted from the alignment files; second, testable AS events are selected; third, SJ read counts for the selected AS events are modeled; finally, the test results are reported. Most of the methods also include a visualization of the results to facilitate interpretation.

## Current challenges of DJU methods
However, these state-of-the-art DJU methods produce results with a low agreement with each other. 
To demonstrate that, we included the Spike-in RNA Variant Control Mixes (SIRVs) in a recent dataset with knockdown and knockout of the CASC3 gene [@Gerbracht_2020]. 
By computing DJU with JunctionSeq, Majiq and, Leafcutter, we report little intersection of significant identified genes and SJ on the SIRV artificial transcriptome, for which the ground truth is known (see the [Benchmark](benchmark.md)). 
The lack of consensus among the different tools represents a barrier to users who want to compare the various methods to select alternatively spliced exon-exon junctions to be experimentally validated.

## The aim of this project
Despite the increase in the number of known AS events, there are still many challenges to further our understanding of the molecular mechanism underlying tissue and disease-specific RNA splicing. 
The lack of consensus among the tools leads to difficulties in selecting which junctions are biologically relevant. This issue challenges the detection and validation of AS events. 
Here, we introduce Baltica, a framework that simplifies the execution and integration of results from DJU methods and summarizes AS's potential biological consequences from changes in the annotation.
Besides presenting Baltica, we also provide a benchmark of the three DJU methods using the Spike-In RNA Variants (SIRVs) as ground-truth for alternative splicing detection.

\bibliography
