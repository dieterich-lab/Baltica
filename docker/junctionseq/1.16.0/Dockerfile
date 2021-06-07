# syntax=docker/dockerfile:1
FROM bioconductor/bioconductor_docker:RELEASE_3_10
LABEL maintainer="Thiago Britto-Borbes <thiagobrittoborges@uni-heidelberg.de>"

RUN apt-get update

RUN R -e 'install.packages("http://hartleys.github.io/QoRTs/QoRTs_LATEST.tar.gz", repos=NULL, type="source")'

RUN R -e 'BiocManager::install("JunctionSeq", dependencies = T, ask=F)' 

