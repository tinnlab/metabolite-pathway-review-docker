ARG BASE_IMAGE=metabolomics_benchmark-base:latest
FROM $BASE_IMAGE

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    libxt-dev \
    ca-certificates-java \
    openjdk-11-jre-headless \
    r-cran-rjava \
    libgit2-dev && \
    apt-get clean

WORKDIR /code/
COPY . /code/
RUN Rscript ./install_packages.R