# syntax=docker/dockerfile:1

FROM bioconductor/bioconductor_docker:latest

RUN apt-get update \
    && pip3 install ttp \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN install2.r --error \
    tidyverse \
    reshape2 \
    rmarkdown \
    openxlsx \
    reactable \
    knitr \
    ggrepel \
    scales \
    yaml \
    pander \
    rmdformats \
    optparse \
    Rcpp \
    RColorBrewer \ 
    igraph

RUN R -e 'BiocManager::install(c("GenomicRanges", "GenomicFeatures", "ComplexHeatmap"))'
