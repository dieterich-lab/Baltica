# syntax=docker/dockerfile:1

FROM continuumio/miniconda3:4.10.3

COPY envs/baltica.yml /opt/env.yaml

RUN conda env create -f /opt/env.yaml

RUN conda init bash

RUN echo "conda activate baltica" >> ~/.bashrc
