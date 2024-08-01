ARG BASE_IMAGE=metabolomics_benchmark-base:latest
FROM $BASE_IMAGE

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    libxt-dev \
    libtirpc-dev \
    python3-scipy \
     && \
    apt-get clean

#dont know why but we need setuptools>51 to install rpy2 successfully.
RUN pip3 install --upgrade pip wheel setuptools>51 
RUN pip3 install --no-cache-dir rpy2 matplotlib lxml

WORKDIR /code/
COPY . /code/

RUN Rscript ./install_packages.R

# Install Python packages
RUN /usr/bin/yes | pip3 uninstall --no-cache-dir scipy
RUN pip3 install --no-binary scipy scikit-learn