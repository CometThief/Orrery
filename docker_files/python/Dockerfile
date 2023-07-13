# Using a Conda-ready Docker image
FROM continuumio/miniconda3:latest

# Create work directory
RUN mkdir -p /app/Orrery
WORKDIR /app/Orrery

# Install system dependencies
RUN apt-get update && apt-get install -y \
    curl \
    tar \
    && rm -rf /var/lib/apt/lists/*

# Copy files into the container
COPY ./environment.yml ./environment.yml
COPY ./installer_files/mgltools_x86_64Linux2_1.5.7.tar ./installer_files/

# Create the Python 3 environment
RUN conda env create -f environment.yml
RUN /opt/conda/envs/Orrery/bin/pip install pybel
RUN /bin/bash -c ". activate Orrery && conda install -c conda-forge qvina"

# Create Python 2.7 environment
RUN conda create -n mgltools python=2.7

# Install MGLTools into mgltools env (via MGLTools' predownloaded installer)
# this needs prettifying
SHELL ["/bin/bash", "-c"]
RUN source activate mgltools
SHELL ["/bin/sh", "-c"]
RUN tar -xvf ./installer_files/mgltools_x86_64Linux2_1.5.7.tar \
    && cd mgltools_x86_64Linux2_1.5.7 \
    && ./install.sh -d /usr/local \
    && ldconfig
    
# Add MGLTools scripts to PATH
ENV PATH="/usr/local/MGLToolsPckgs/AutoDockTools/Utilities24/:${PATH}"
ENV PYTHONPATH="$PYTHONPATH:/usr/local/MGLToolsPckgs"

# Extra installations
RUN /opt/conda/envs/mgltools/bin/pip install numpy

# Default start of correct env when starting a bash shell (useful with docker exec)
RUN echo "source activate Orrery" >> ~/.bashrc
