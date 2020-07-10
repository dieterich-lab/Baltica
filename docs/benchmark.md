## Baltica benchmark with the SIRVome

We implemented a benchmark method for DJU methods using the Spike-In RNA Variants (SIRVs Set-1, cat 025.03[^1]) as ground-truth for alternative splicing identification.
Splicing identifiation methods often use simulated data that not fully includes the complexity of a RNA-Sequencing experiment, and our approach aims to overcome this limitation.  

The SIRVs spike-in comprise 7 genes, 69 transcript isoforms, 357 exons and 113 introns that are used in different concetrantion in the 3 mixes. For this particular experiment, we designed the SIRV addition to not confound with the biological factors (TODO link figure). SJ differently used among mixes in the SIRV contig were considered true positive calls, while other significant call were considered false positives. 

(TODO link figure)
Detail on the result

More details on the finding can be found on our manuscript: TODO

[^1]: https://www.lexogen.com/sirvs/#sirvsdownload - Lexogen took no part on design of this experiment and we have not relationshop with this company


