# syntax=docker/dockerfile:1
FROM python:3.6.12-buster

LABEL maintainer="Thiago Britto-Borbes <thiagobrittoborges@uni-heidelberg.de>"

RUN apt-get update && \
    apt-get install -y \
    build-essential \
    cmake \
    git \
    gfortran \
    libblas-dev \
    libgsl-dev \
    liblapack-dev && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* && \
    pip install "Cython==0.29.21" "numpy==1.16.6"


RUN git clone https://github.com/Xinglab/rmats-turbo.git && \
    cd rmats-turbo && \
    ./build_rmats --no-paired-model 


WORKDIR rmats-turbo
ENTRYPOINT ["./run_rmats"]
CMD ["--help"]