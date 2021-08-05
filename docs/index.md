# [Baltica: integrated splice junction usage analysis](https://github.com/dieterich-lab/Baltica)

![Baltica logo](https://gist.githubusercontent.com/tbrittoborges/3c86ffbaa62e671771f443c65cb04fdc/raw/7ae0ea4a76e8f5464139ef34164c67de7a297ce8/baltica_logo.png
){ : .right width=220px loading=lazy  }

Baltica is a framework that facilitates the execution and enables the integration of results from multiple differential junction usage (DJU) methods. The core of the framework is Snakemake workflows [@M_lder_2021], a python command-line interface, and R/Bioconductor scripts for analysis [@r_core][@Lawrence_2013][@Lawrence_2009][@Wickham_2019]. The workflows are include methods for RNA-Seq quality control [@wang2012][@andrews2012][@ewels_2016], four DJU methods: RMATs [@Shen_2014] JunctionSeq [@Hartley2016], Majiq [@VaqueroGarcia2016] and Leafcutter [@Li2017]. We use Stringtie2 [@Kovaka_2019] _de novo_ transcriptome assembly to re-annotate the results. Baltica's main goal is to provide an integrative view of the results of these methods. To do so,  Baltica produces an RMarkdown report with the integrated results and links to UCSC GenomeBrowser for further exploration.

## Features
    - Snakemake workflows for DJU: junctionseq, majiq, rmats, and leafcutter
    - Snakemake workflow for de novo transcriptome annotation with stringtie
    - Process, integrate and annotate the results from the methods
    - Summarise AS class of differently spliced junctions
    - DJU method benchmarks
    - Report on the integrative analysis

**To get started**, use the menu on the left-hand side or search function to navigate over this documentation. 
  
## Citation

Thiago Britto-Borges, Volker Boehm, Niels H. Gehring and Christoph Dieterich (2020) __Baltica: integrated splice junction usage analysis__. 
Manuscript in preparation.

Baltica is based on the work of many scientists and developers. Thus, if you use the results of their tools in your analysis, consider citing their work.

## License
Baltica is free, open-source software released under an [MIT License](https://github.com/dieterich-lab/Baltica/blob/master/LICENSE).

## Contact
Please get in touch with us [the GitHub issue tracker](https://github.com/dieterich-lab/Baltica/issues).

## References
\bibliography