# Baltica

Baltia enables integrated differential junction usage (DJU) and consequence analysis with an ensemble of methods. 

Alternative splicing (AS) is a tightly regulated co- and post-transcriptional process that contributes to the transcriptome diversity observed in eukaryotes.
We noted that several methods for detecting differential junction usage (DJU) from RNA-Sequencing data exist yet deliver results that show limited agreement.  

Here we present Baltica, a framework that enables the integration of DJU methods results and provides an extra annotation to prioritize splice junctions (SJ) for experiments validation by associating the AS with its biological consequence.  

Baltica provides workflows for quality control and supports three DJU methods, Junctionseq, Majiq, and Leafcutter. The framework uses the ensemble methods to decide which junctions are called significant and sort SJ with AS events that change a protein annotation to select the ones with higher potential for protein function change.  

!!! info
    To get started, use the menu on the left-hand side or search function to navigate over this documentation.

# Citation
TODO

Baltica comprises of workflows for multiple DJU tools and scripts for downstream analysis of the results from these tools. The DJU tools represent several years of work from research scientists, and so if you use the results of any of these tools, please consider citing them. TODO citation page

# License
Baltica is free, open-source software released under an [MIT License](https://github.com/dieterich-lab/Baltica/blob/master/LICENSE).

# Contact
Please contact us via [the GitHub issue tracker](https://github.com/dieterich-lab/Baltica/issues)