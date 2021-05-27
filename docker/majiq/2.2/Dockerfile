# syntax=docker/dockerfile:1
FROM python:3.7
LABEL maintainer="Thiago Britto-Borbes <thiagobrittoborges@uni-heidelberg.de>"

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y git

RUN wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
    tar -vxjf htslib-1.9.tar.bz2 && rm htslib-1.9.tar.bz2 && \
    cd htslib-1.9/ && \
    make && make install && \
    export HTSLIB_LIBRARY_DIR=/usr/local/lib && \
    export HTSLIB_INCLUDE_DIR=/usr/local/include

RUN cd .. && \
    git clone https://bitbucket.org/biociphers/majiq_academic.git  && \
    cd majiq_academic && \
    pip install --no-cache-dir -r requirements.txt && \
    python setup.py install
