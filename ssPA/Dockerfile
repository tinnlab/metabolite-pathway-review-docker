ARG BASE_IMAGE=metabolomics_benchmark-base:latest
FROM $BASE_IMAGE

RUN pip3 install --no-cache-dir tqdm==4.66.4 sspa==1.0.2

WORKDIR /code/
COPY . /code/

