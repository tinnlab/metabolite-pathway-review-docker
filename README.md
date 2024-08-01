# A Review on Functional Analysis Tools For Metabolites (Docker Images)

This repository contains Dockerfiles for all methods in the review work on functional analysis on metabolomics data

# QUICK START
Users should make sure they have [Docker Desktop](https://www.docker.com/products/docker-desktop/) installed already.

**Step 1:** we need to clone this repository and place it in a folder in our local machine
```
# If you are users with a Linux/MAC computer, open your terminal
# If you are users with a Windows computer, open your Windows Powershell
cd Desktop                              
git clone https://github.com/tinnlab/metabolite-pathway-review-docker
cd metabolite-pathway-review-docker
```

**Step 2:** To build the docker image for a method of your interest, please run the following command:
```
# build base image first
docker-compose build metabolomics_benchmark-base

# build method image
docker-compose build <*method-name*>

# run the method image
docker-compose run <*method-name*> bash
```
