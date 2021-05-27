# syntax=docker/dockerfile:1
FROM bioconductor/bioconductor_docker:latest
LABEL maintainer="Thiago Britto-Borbes <thiagobrittoborges@uni-heidelberg.de>"

ENV LANG=C.UTF-8

RUN apt-get update && \
    apt-get install -y \
    samtools \
    cmake \
    build-essential \
    && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/*

RUN git clone https://github.com/griffithlab/regtools && \
    cd regtools && \
    mkdir build && cd build && \
    cmake .. && \
    make && \
    cd ..

RUN R -e "install.packages(c('devtools', 'rstan', 'rstantools'), repos = 'https://cloud.r-project.org/', dependencies = T)"

RUN R -e 'BiocManager::install(c("DirichletMultinomial", "Biobase"), dependencies = T, ask=F)'

RUN R -e "devtools::install_github('davidaknowles/leafcutter/leafcutter', ask=F, dependencies = T)"


