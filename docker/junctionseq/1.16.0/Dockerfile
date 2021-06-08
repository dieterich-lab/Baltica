# syntax=docker/dockerfile:1
FROM bioconductor/bioconductor_docker:RELEASE_3_10
LABEL maintainer="Thiago Britto-Borbes <thiagobrittoborges@uni-heidelberg.de>"

RUN echo 'deb http://deb.debian.org/debian/ sid main' | sudo tee /etc/apt/sources.list.d/docker.list && \
    apt-get update && \
    apt-get install -y --no-install-recommends openjdk-11-jre && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/*

RUN wget http://hartleys.github.io/QoRTs/QoRTs.jar && \
    mv QoRTs.jar /usr/local/bin/ && \
    R -e 'install.packages("http://hartleys.github.io/QoRTs/QoRTs_LATEST.tar.gz", repos=NULL, type="source")' &&\
    echo $'#!/bin/sh\nprefix=`dirname $(readlink $0 || echo $0)`\nexec java -Xmx16000m -XX:-UseGCOverheadLimit -jar QoRTs.jar "$@"' > /usr/local/bin/qorts


RUN R -e 'BiocManager::install("JunctionSeq", dependencies = T, ask=F)' 