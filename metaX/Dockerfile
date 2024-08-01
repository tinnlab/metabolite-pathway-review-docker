ARG BASE_IMAGE=metabolomics_benchmark-base:latest
FROM $BASE_IMAGE

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    libgit2-dev \
    libnlopt-dev && \
    apt-get clean

WORKDIR /code/
COPY . /code/

RUN Rscript ./install_packages.R


