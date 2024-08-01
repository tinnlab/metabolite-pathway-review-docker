# Stage 1: Setup the Conda environment with Miniconda3
FROM continuumio/miniconda3:latest as conda-stage

WORKDIR /app

# Create the conda env
COPY environment.yml .
RUN conda env create -f environment.yml

# Stage 2: Use my base image for other dependencies
FROM r-base:latest

RUN mkdir -p /code/
WORKDIR /code/
COPY . /code/

# Copy the Conda env from the first stage to the same location in the final image being built from the base image.
COPY --from=conda-stage /opt/conda/envs/PUMA /opt/conda/envs/PUMA

# Set the PATH environment variable to use the Conda environment
ENV PATH /opt/conda/envs/PUMA/bin:$PATH


