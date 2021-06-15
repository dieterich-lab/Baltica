
# syntax=docker/dockerfile:1
FROM mambaorg/micromamba:0.13.1

COPY envs/qc.yml /opt/qc.yml

RUN micromamba install -y -n base -f /opt/qc.yml && \
    micromamba clean --all --yes

ENV URL=http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/
RUN wget ${URL}gtfToGenePred && chmod +x gtfToGenePred && mv gtfToGenePred /usr/local/bin/gtfToGenePred
RUN wget ${URL}genePredToBed && chmod +x genePredToBed && mv genePredToBed /usr/local/bin/genePredToBed