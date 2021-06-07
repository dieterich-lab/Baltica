# syntax=docker/dockerfile:1
FROM bioconductor/bioconductor_docker:RELEASE_3_10
LABEL maintainer="Thiago Britto-Borbes <thiagobrittoborges@uni-heidelberg.de>"

RUN apt-get update

RUN R -e 'BiocManager::install("JunctionSeq", dependencies = T, ask=F)' 

