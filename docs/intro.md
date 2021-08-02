# Introduction

For motivation and review of state of the art, please check:

Thiago Britto-Borges, Volker Boehm, Niels H. Gehring and Christoph Dieterich (2020) Baltica: integrated splice junction usage analysis. Manuscript in preparation.

## Tips on RNA-Seq aiming differential splicing detection
If you aim to resolve mRNA isoforms with relatively low abundance, you should design the RNA-seq experiment accordingly. [The expert suggestion](https://support.illumina.com/bulletins/2017/04/considerations-for-rna-seq-read-length-and-coverage-.html) is to sequence around 40 to 60 million reads pairs. This parameter is particularly relevant for complex RNA libraries, but it can be insufficient to saturate novel SJ, in our experience. Read length and paired-end reads are also critical for SJ  identification, and longer reads offer more coverage of the exons boundaries. Thus, the target read length should be around 100 nucleotides for Illumina RNA-seq to maximize the read overhang length and, consequently, maximize the quality of the alignments. 

In the Baltica manuscript, we propose an approach to integrate DJU results from Illumina to DJU results from third-generation sequencing. 

Also, databases such as the [CHESS](http://ccb.jhu.edu/chess/) can provide additional evidence for splice sites absent in the annotation.

\bibliography
