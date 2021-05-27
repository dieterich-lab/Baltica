# syntax=docker/dockerfile:1
FROM debian:latest

ENV LANG en_US.utf8

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    locales \
    samtools \
    libz-dev \
    ca-certificates \
    build-essential \
    make \
    git && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/*


RUN git clone https://github.com/gpertea/stringtie --branch v2.1.5 --single-branch && \
    cd stringtie && \
    make release && \
    cd 

RUN git clone https://github.com/gpertea/gffcompare --branch v0.12.2 --single-branch && \
    cd gffcompare && \
    make release && \
    cd 

ENV PATH="gffcompare/:${PATH}"
ENV PATH="stringtie/:${PATH}"
