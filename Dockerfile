# Use the latest LTS version of Ubuntu as the base image
FROM ubuntu:22.04

# Set environment variables to prevent prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Update package list and install necessary packages
RUN apt-get update && apt-get install -y \
    gfortran \
    mpich \
    libblas-dev \
    wget \
    git \
    make \
    vim \
    unzip \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Set a working directory inside the container
WORKDIR /workspace

# Create a bin directory
RUN mkdir -p /workspace/bin
RUN echo 'export PATH="$PATH:/workspace/bin"' >> /root/.bashrc

# Download the specified file
RUN wget https://www.ssisc.org/lis/dl/lis-2.1.6.zip

# Unzip the downloaded file
RUN unzip lis-2.1.6.zip && rm lis-2.1.6.zip

# Install LIS library
RUN cd lis-2.1.6  && ./configure --enable-mpi --enable-f90 && make && make install 

# Download Andromeda
RUN git clone https://github.com/WildSmilodon/Andromeda.git

# Compile Andromeda
RUN cd Andromeda/src && make
RUN cd Andromeda/src && make install

# Default command to start a bash shell
CMD ["/bin/bash"]
