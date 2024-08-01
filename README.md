# Functional Analysis Tools For Metabolites: A Review (Docker Images)

This repository hosts Dockerfiles for all non-web-based methods covered in the review on functional analysis of metabolomics data.

# QUICK START
Users should make sure they have [Docker Desktop](https://www.docker.com/products/docker-desktop/) installed already.

**Step 1:** Clone this repository and save it to a folder on your local machine
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

# build a certain method image
docker-compose build <method-name>

# run the method image
docker-compose run --rm <method-name> bash
```
**NOTE**: 

+ Users are free to replace `<method-name>` with one of the following tool names: `fella`, `lilikoi`, `metabolyzer`, `metax`, `mummichog`, `metaboanalyst`, `papi`, `puma`, `sspa`, `mbpls`, `ogpls`, `imsea`, or `dcimsea`.

+ We **DO NOT** recommend users to build all method images in one shot (`docker-compose build`) because it will drain your memory. Just build one tool at a time, and remove it (`docker rmi <method-name>`) before coming to another tool.

+ For methods that require MATLAB to run (`MB-PLS`, `iMSEA`, and `ogPLS`), Please visit the [Matlab license center](https://www.mathworks.com/licensecenter/licenses) to obtain the `license.lic` file beforehand. The hostid of this license is the same as the hostid of the host machine. If docker is used, the user of this license must be root. After you had that license file, place it in the `license` sub-folder in the respective folder (i.e., `MBPLS/license`, `iMSEA/license`, or `ogPLS/license`).
