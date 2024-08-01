ARG BASE_IMAGE=metabolomics_benchmark-base:latest
FROM $BASE_IMAGE

RUN pip3 install --no-cache-dir networkx==2.6.3 xlrd scipy==1.14.0 scikit-learn openpyxl==3.1.5 tqdm==4.66.4 matplotlib

WORKDIR /code/
COPY . /code/