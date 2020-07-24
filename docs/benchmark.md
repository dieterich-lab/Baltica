# Baltica benchmark with SIRV transcripts

We implemented a benchmark method for DJU methods using the Spike-In RNA Variants (SIRVs Set-1, cat 025.03[^1]) as ground-truth for alternative splicing identification.
Splicing identifiation methods often use simulated data for benchmark, which does not fully appreciate the complexity of RNA-Sequencing experiment. We use a complementar approach that aims to overcome this limitation.  

The SIRVs spike-in comprise 7 genes, 69 transcript isoforms, 357 exons and 113 introns that are used in different concetrantion in the 3 mixes. For this particular experiment, we designed the SIRV addition to not confound with the biological factors . SJ differently used among mixes in the SIRV contig were considered true positive calls, while other significant call were considered false positives. 

![](img/experiment_design.png)
* __Experimental design__ : ... *


<!-- 
    ROC CURVE 
    PR CURVE 
    INTRON HEATMAP 
    SCORE RANKS -->
## RNA-Seq processing and mapping

Cell lines, RNA extraction, and RNA-Seq were described by [@Gerbracht_2020]. In short, we obtained 15 libraries from Flp-In T-REx 293 cells, extracted the RNA fraction, with TrueSeq Stranded Total RNA kit (Illumina), followed by ribosomal RNA depletion, with RiboGold Plus kit, and reads were sequenced with an Illumina -HiSeq4000- using PE 100bp protocol and yield around 50 million reads per sample. Data was deposited in ArrayExpress (E-MTAB-8461).
Sequenced reads' adapters and low-quality bases were trimmed, and reads mapping to human precursor ribosomal RNA were discarded.
The remaining reads were aligned to the human transcriptome (version 38, EnsEMBL 90) extended with the SIRV annotation.
Regarding the DJU method benchmarking, we are not interested in the biological condition, but the SIRVs AS event, and so this experiment was designed, so the SIRVs Set-1 (cat 025.03) mixes were not confounded to the biological factors.


[^1]: https://www.lexogen.com/sirvs/#sirvsdownload - Lexogen took no part on design of this experiment and we have not relationshop with this company


The SIRVs spike-in comprise 7 genes, 69 transcript isoforms, 357 exons and 113 intron. \todo{XXX} introns are differently used in the 3 mixes. According to our experimental design, differences in transcript abbundance lead to differetial splicing events in SIRV chromosomes but not in the human contigs, that in this case represent false positive calls if callsed significant. \textbf{(A)} the ROC curve and \textbf{(B)} the precision-recall curve. Overall, in this case where the annotation faithfully represents the transcript structure, JunctionSeq outperforms Majiq and Leafcutter. In this case, Leafcutter suffers from a recall cealing issue, because it can't identify the certain introns from the artifical transcripts.
