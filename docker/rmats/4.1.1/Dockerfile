# syntax=docker/dockerfile:1
FROM mambaorg/micromamba:0.13.1

COPY envs/rmats.yml /root/rmats.yml

RUN micromamba install -y -n base -f /root/rmats.yml && \
    micromamba clean --all --yes
