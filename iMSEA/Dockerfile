ARG BASE_IMAGE=metabolomics_benchmark-base:latest
FROM $BASE_IMAGE

ARG MATLAB_RELEASE=r2022b

# Install necessary dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    ca-certificates \
    wget \
    unzip \
    sudo \
    libx11-6 \
    procps && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Run mpm to install MATLAB in the target location and delete the mpm installation afterwards
RUN wget -q https://www.mathworks.com/mpm/glnxa64/mpm && \
    chmod +x mpm && \
    ./mpm install \
        --release=${MATLAB_RELEASE} \
        --destination=/opt/matlab \
        --products MATLAB Optimization_Toolbox Global_Optimization_Toolbox Statistics_and_Machine_Learning_Toolbox Parallel_Computing_Toolbox&& \
    rm -f mpm /tmp/mathworks_root.log && \
    ln -s /opt/matlab/bin/matlab /usr/local/bin/matlab

# Add "matlab" user and grant sudo permission.
RUN adduser --shell /bin/bash --disabled-password --gecos "" matlab && \
    echo "matlab ALL=(ALL) NOPASSWD: ALL" > /etc/sudoers.d/matlab && \
    chmod 0440 /etc/sudoers.d/matlab

# Set environment variables for Matlab
ENV MW_DDUX_FORCE_ENABLE=true MW_CONTEXT_TAGS=MATLAB:DOCKERFILE:V1

COPY . /code/

# Ensure the script is executable
RUN chmod +x run.sh

# Default command to run the script
CMD ["/bin/sh", "/code/run.sh"]