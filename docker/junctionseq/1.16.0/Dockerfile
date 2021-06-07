# syntax=docker/dockerfile:1
FROM bioconductor/bioconductor_docker:RELEASE_3_10
LABEL maintainer="Thiago Britto-Borbes <thiagobrittoborges@uni-heidelberg.de>"

RUN apt-get update

RUN R -e 'BiocManager::install("JunctionSeq", repos = "https://cloud.r-project.org/", dependencies = T, ask=F)' 

