ARG BASE_IMAGE=metabolomics_benchmark-base:latest
FROM $BASE_IMAGE

WORKDIR /code/
COPY . /code/

RUN Rscript ./install_packages.R