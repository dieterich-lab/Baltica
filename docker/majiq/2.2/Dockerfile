# syntax=docker/dockerfile:1
FROM python:3.7
LABEL maintainer="Thiago Britto-Borbes <thiagobrittoborges@uni-heidelberg.de>"

ENV HTSLIB_VER=1.9
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y git

RUN wget https://github.com/samtools/htslib/releases/download/$HTSLIB_VER/htslib-$HTSLIB_VER.tar.bz2 && \
    tar -vxjf htslib-$HTSLIB_VER.tar.bz2 && rm htslib-$HTSLIB_VER.tar.bz2 && \
    cd htslib-$HTSLIB_VER/ && \
    make && make install && \
    export HTSLIB_LIBRARY_DIR=/usr/local/lib && \
    export HTSLIB_INCLUDE_DIR=/usr/local/include

RUN cd .. && \
    git clone https://bitbucket.org/biociphers/majiq_academic.git  && \
    cd majiq_academic && \
    pip install --no-cache-dir -r requirements.txt && \
    python setup.py install

RUN wget https://majiq.biociphers.org/download/gtf2gff3.pl &&\
    mkdir scripts && mv gtf2gff3.pl scripts/

ENV PATH="scripts/:${PATH}"
