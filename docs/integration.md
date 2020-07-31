## Parsing the results of the method 

The first step in the analysis workflow is parsing and processing the DJU methods' output with `scripts/parse_{method}_output.R` scripts as follows:

1. The resulting text output from the DJU methods is parsed and loaded as R data frames
1. The data frames are pivoted in a longer format to have one junction and one comparison per row 
1. Finally, SJ not called significant are discarded 

The default significant cut-off is an adjusted p-value < 0.05 for JunctionSeq and Leafcutter or probability of changing > 0.90 for Majiq.

The output results in `{method}/{method}_junction.csv` are in a `tidy` format and can be used for downstream analysis. 

!!!note
    Although it is generally good to filter results with small effect sizes, discarding results with small deltaPSI can be problematic. First, it's not trivial to assign deltaPSI value to JunctionSeq results, since the tool does not detect splicing events. Second, SJ with small deltaPSI may indicate RNA degradation in a cell with intact degradation machinery. 


## Annotating the results 

We annotate the results with information from the gene and transcript hosting the SJ. For this, we use the _de novo_ transcript annotation at `stringtie/merged/merged.combined.gtf`. It's common the multiple transcripts share an intron, and so a single intron may be annotated with multiple transcripts.

One challenge for the integration of DJU methods' results is the different genomic coordinates system. The coordinates system's differences are due to the method implementation design: methods can be 0-indexed (BED format) versus 1-indexed (GTF format) or use the exonic versus intronic coordinates to represent the SJ genomic position.   
We propose a `filter_hits_by_diff` (below) that discard any hits with more than two bp differences from overlapping features between the reference (query) and the SJ (subject). Using introns obtained in that annotation as a reference, and `filter_hits_by_diff`, Baltica enables the reconciliation of the multiple DJU results.

```R
filter_hits_by_diff <- function(query, subject, max_start = 2, max_end = 2) {
  stopifnot(is(query, "GRanges"))
  stopifnot(is(subject, "GRanges"))
  hits <- findOverlaps(query, subject)
  query <- query[queryHits(hits)]
  subject <- subject[subjectHits(hits)]
  start_dif <- abs(start(query) - start(subject))
  end_dif <- abs(end(query) - end(subject))
  hits <- hits[start_dif <= max_start & end_dif <= max_end]
  hits
}
```

These are the columns assigned after the annotation:

Column name | Description | Note
------------|-------------|------
comparison | pairwise comparison as `{case}_vs_{control}` |
chr | seqname or genomic contig |
start | intron start position for the SJ | 1-index
end |  intron end position | 1-index
strand | RNA strand that encodes that gene |
gene | the gene symbol | 
e2| acceptor exon number, if in + strand otherwise donor exon | 
e1 | donor exon number, if in the + strand otherwise acceptor exon |
tx_id | transcript identifier from the combined annotation |
transcript_name | transcript name | 
class_code | association between reference transcript and novel transcript | [see fig 1 from this paper for details]( https://doi.org/10.12688/f1000research.23297.1)

__Table 6.1: Columns added after annotation __

## Selectin optimal _de novo_ transcriptome parameters
We found that the parameters used to obtain the _de novo_ transcriptome are critical for maximum integration between the GTF and the SJ from DJU methods. __Fig 6.1__ shows a parameter scan where we vary the group, `-j` (minimum junction coverage), `-c` (minimum coverage), and `-f` (minimum isoform proportion) and compute the number of transcripts that match with SJ called significant. As expected, the merged annotation and not the group-specific annotation have the highest rate of annotated introns. The crucial result here is the dependency of the `-f` parameter, which is also associated with an increased number of annotated introns. As we confirmed this behavior in other datasets, we decided to use `-c 3 -j 3 -f 0.01` as default values in Baltica. The higher coverage (`-c` and `-j`) values counter the potential noise of transcripts with low abundance.

![](img/stringtie_parameter_scan_heatmap.png)  

__ Fig. 6.1: Parameters scan to maximize the number of introns annotated __ We have run Stringtie with multiple for group annotations and merged annotation. We use a series of parameters: junction coverage of 1, 2, or 3; coverage of 1, 2, 3, and minimum isoform fraction of 0.1, 0.01, or 0.001.  

## Assigning AS type

Identifying the type of AS is critical to understand a potential molecular mechanism for AS events. Many factors influence splicing. Among then, splicing factors probably the most well-studied and can have a specific function in splicing regulation. [SRSF2](https://www.uniprot.org/uniprot/Q01130) is a relevant example in this context. SRSF2 is splicing factors from the SR family that are known for auto-regulation. In certain conditions, the SRSF2 transcript can activate the nonsense-mediated decay by either including a new exon containing a premature stop codon or an intron in it's 3' UTR. These changes lead to the transcript degradation and overall reduction of the gene expression. SRSF2 increased the probability of exons with weak splice site signals to be included in the transcript. Thus, the reduction of SRSF2 protein level leads to widespread exon skipping. Identifying such patterns is critical to understanding which splicing regulators are driving the observed splicing changes.

In Baltica, we use a geometric approach to define AS in three classes:
- ES, for exon skipping
- A3SS, for alternative 3' splice-site
- A5SS, for alterantive 5' splice-site

Figure 6.2 details how we use the distance between features start and end to determine the AS type.

![](img/Baltica_as_type.png)  

__ Fig. 6.2: AS type assignment in Baltica. __ Baltica uses the genomic coordinates from the SJ and its overlapping exons to assing AS type to SJ and it's overlapping exons. Because many exons may be affected, multiple assignments are output. Donor and acceptor exons are assigned as JS and JE, respectively. 

## Simplify the AS event
Because most of the final users are only interested in the list of genomic ranges, gene names, or event types, we offer a simplified output that removes the excess of redundant information. This step is useful for generating a final report. 

Below the `{r} sessionInfo()` output for packages loaded for the analysis workflow.

```
R version 3.6.0 (2019-04-26)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Debian GNU/Linux 8 (jessie)

Matrix products: default
BLAS:   /beegfs/biosw/R/3.6.0/lib/R/lib/libRblas.so
LAPACK: /beegfs/biosw/R/3.6.0/lib/R/lib/libRlapack.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets
[8] methods   base

other attached packages:
 [1] openxlsx_4.1.5       reshape2_1.4.4       rtracklayer_1.46.0
 [4] GenomicRanges_1.38.0 GenomeInfoDb_1.22.1  IRanges_2.20.2
 [7] S4Vectors_0.24.4     BiocGenerics_0.32.0  optparse_1.6.6
[10] dplyr_0.8.5          readr_1.3.1          stringr_1.4.0
[13] tidyr_1.0.3

loaded via a namespace (and not attached):
 [1] zip_2.0.4                   Rcpp_1.0.4.6
 [3] plyr_1.8.6                  pillar_1.4.3
 [5] compiler_3.6.0              XVector_0.26.0
 [7] bitops_1.0-6                tools_3.6.0
 [9] zlibbioc_1.32.0             lattice_0.20-38
[11] lifecycle_0.2.0             tibble_3.0.1
[13] pkgconfig_2.0.3             rlang_0.4.6
[15] Matrix_1.2-17               DelayedArray_0.12.3
[17] cli_2.0.2                   GenomeInfoDbData_1.2.2
[19] Biostrings_2.54.0           vctrs_0.2.4
[21] hms_0.5.3                   grid_3.6.0
[23] getopt_1.20.3               tidyselect_1.0.0
[25] Biobase_2.46.0              glue_1.4.0
[27] R6_2.4.1                    fansi_0.4.1
[29] BiocParallel_1.20.1         XML_3.99-0.3
[31] purrr_0.3.4                 magrittr_1.5
[33] matrixStats_0.56.0          GenomicAlignments_1.22.1
[35] Rsamtools_2.2.3             ellipsis_0.3.0
[37] SummarizedExperiment_1.16.1 assertthat_0.2.1
[39] stringi_1.4.6               RCurl_1.98-1.2
[41] crayon_1.3.4

```

\bibliography