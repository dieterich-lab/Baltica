# syntax=docker/dockerfile:1
FROM bioconductor/bioconductor_docker:latest

ENV LANG en_US.utf8

RUN R -e "install.packages(c('GenomicRanges', 'tidyverse', 'openxlsx', 'optparse'), repos = 'https://cloud.r-project.org/', dependencies = T, ask=F)"