# syntax=docker/dockerfile:1
################# BASE IMAGE ######################
FROM python:3.6.12-buster

################## METADATA ######################
LABEL maintainer="Thiago Britto-Borbes <thiagobrittoborges@uni-heidelberg.de>"

################## INSTALLATION ######################
RUN apt-get update && \
    apt-get install -y \
    unzip \
    default-jdk \
    libfindbin-libs-perl \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* && 

RUN pip install "multiqc==v1.10.1" "RSeQC==4.0.0" 

ENV URL="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/"
ENV FASTQC=fastqc_v0.11.9.zip  

RUN wget ${URL}${FASTQC} \
    && unzip ${FASTQC} \
    && rm ${FASTQC} \
    && chmod 755 /FastQC/fastqc

ENV PATH "/FastQC:$PATH"