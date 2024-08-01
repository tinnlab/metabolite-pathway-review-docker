FROM ubuntu:22.04
ARG DEBIAN_FRONTEND=noninteractive

LABEL maintainer="Quang-Huy Nguyen <huynguyen96.dnu@gmail.com>. My Google Scholar: <https://scholar.google.com/citations?user=LokoolIAAAAJ&hl=en>"
USER root  

# Install R and its keys
RUN apt-get update -qq && \
    apt-get install -y --no-install-recommends software-properties-common dirmngr gnupg2 wget && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* 
#remove cached files (/var/lib/apt/lists/*) after installing dependencies to reduce the image size.

RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc \
  | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc \
  && add-apt-repository -y "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" \
  && apt-get update -y \
  && apt-get install -y r-base r-base-core r-recommended r-base-dev r-cran-dplyr r-cran-devtools r-cran-rcpp

# Update system and install prerequisites
RUN apt-get update && \
    apt-get install -y --no-install-recommends  \
    build-essential \
    liblapack-dev  \
    libopenblas-dev \
    gfortran \
    libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    libcurl4-openssl-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libxml2 libxml2-dev libcairo2-dev \
    libpq-dev libssh2-1-dev unixodbc-dev libssl-dev libudunits2-dev libgdal-dev libnetcdf-dev \
    libgeos-dev libproj-dev python3 python3-dev python3-pip && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
# remove cached files (/var/lib/apt/lists/*) after installing dependencies to reduce the image size.

# Install Python packages
RUN pip3 install --no-cache-dir numpy pandas

# Create a folder to save code
RUN mkdir -p /code/
WORKDIR /code/
COPY . /code/

RUN Rscript ./install_packages.R

CMD ["echo", "Welcome to Quang-Huy Nguyen's Metabolomics Review!"]