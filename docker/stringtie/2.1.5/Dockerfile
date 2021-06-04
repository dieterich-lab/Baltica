# syntax=docker/dockerfile:1
FROM debian:latest

ENV LANG=C

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

ARG STRINGTIE_VERSION=2.1.5
RUN git clone https://github.com/gpertea/stringtie --branch v$STRINGTIE_VERSION --single-branch && \
    cd stringtie && make release && cd 


ARG GFFCOMPARE_VERSION=0.12.2
RUN git clone https://github.com/gpertea/gffcompare --branch v$GFFCOMPARE_VERSION --single-branch && \
    cd gffcompare && make release && cd

RUN mv stringtie/* /usr/local/bin/ && mv /gffcompare/* /usr/local/bin/ 
RUN rm -rf stringtie/ /gffcompare/
LABEL maintainer='tbrittoborges'