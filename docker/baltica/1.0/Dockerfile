# syntax=docker/dockerfile:1

FROM mambaorg/micromamba:0.13.1

COPY envs/baltica.yml /opt/env.yaml

RUN micromamba install -y -n base -f /opt/env.yaml && \
    micromamba clean --all --yes