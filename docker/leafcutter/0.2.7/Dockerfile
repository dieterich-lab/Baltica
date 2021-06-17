# syntax=docker/dockerfile:1
FROM mambaorg/micromamba:0.13.1

COPY envs/leafcutter.yml /opt/env.yaml

RUN micromamba install -y -n base -f /opt/env.yaml && \
    micromamba clean --all --yes

RUN Rscript -e "install.packages(c('ggplot2', 'reshape2', 'Hmisc', 'dplyr', 'shiny', 'intervals', 'RcppEigen', 'BiocManager'), repos = 'https://cloud.r-project.org/')"

RUN Rscript -e "BiocManager::install(c('DirichletMultinomial', 'Biobase'))" 

ENV LEAFCUTTER_VERSION="0.2.7" 
# Need to wait next release to fix the release version
# python3 clustering scripts are unavaiable 

RUN Rscript -e 'devtools::install_github("davidaknowles/leafcutter/leafcutter", tag="v$LEAFCUTTER_VERSION")'

# RUN git clone https://github.com/davidaknowles/leafcutter --branch v$LEAFCUTTER_VERSION --single-branch && \
#     cd leafcutter && chmod 755 scripts/* clustering/* && \
#     mv scripts/* /opt/conda/bin/ && mv clustering/* /opt/conda/bin/ 
RUN wget https://raw.githubusercontent.com/dieterich-lab/Baltica/master/scripts/{leafcutter_cluster_regtools_py3.py,leafcutter_ds_pair.R,gtf_to_exons.R} && \
    chmod 755 leafcutter_cluster_regtools_py3.py leafcutter_ds_pair.R gtf_to_exons.R  && \
    mv leafcutter_cluster_regtools_py3.py leafcutter_ds_pair.R gtf_to_exons.R /opt/conda/bin/ 
