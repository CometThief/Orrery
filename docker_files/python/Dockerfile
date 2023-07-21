# Using a Conda-ready Docker image
FROM continuumio/miniconda3:latest

# Create work directory
RUN mkdir -p /app/Orrery
WORKDIR /app/Orrery

# Install system dependencies
RUN apt-get update && apt-get install -y \
    curl \
    tar \
    g++ \
    && rm -rf /var/lib/apt/lists/*

# Copy files into the container
COPY ./environment.yml ./environment.yml
#COPY ./installer_files/mgltools_x86_64Linux2_1.5.7.tar ./installer_files/

# Create the Python 3 environment
RUN conda config --add channels defaults \
    && conda config --set remote_max_retries 10 \
    && conda config --set remote_connect_timeout_secs 60 \
    && conda env create -f environment.yml

# Set the PATH and activate Orrery environment
ENV PATH /opt/conda/envs/Orrery/bin:$PATH
RUN echo "source activate Orrery" > ~/.bashrc

# Install additional packages into Orrery environment
RUN pip install pybel \
    && conda install -c conda-forge qvina \
    && pip install torch==1.13.1+cpu -f https://download.pytorch.org/whl/cpu/torch_stable.html \
    && pip install torch-scatter -f https://data.pyg.org/whl/torch-1.13.1+cpu.html \
    && pip install torch-sparse -f https://data.pyg.org/whl/torch-1.13.1+cpu.html \
    && pip install torch-spline-conv torch-geometric==2.0.1 -f https://data.pyg.org/whl/torch-1.13.1+cpu.html \
    && pip install cairosvg svgutils molvs \
    && pip install git+https://github.com/mayrf/pkasolver.git \
    && pip install meeko

# Create Python 2.7 environment
RUN conda create -n mgltools python=2.7

# Install MGLTools into mgltools env (via MGLTools' predownloaded installer)
SHELL ["/bin/bash", "-c"]
RUN source activate mgltools
SHELL ["/bin/sh", "-c"]
RUN wget -q https://ccsb.scripps.edu/download/532/ -O mgltools.tar.gz && \
    tar -xzf mgltools.tar.gz && \
    rm mgltools.tar.gz && \
    cd mgltools_x86_64Linux2_1.5.7 && \
    ./install.sh -d /usr/local && \
    ldconfig

# Add MGLTools scripts to PATH
ENV PATH="/usr/local/MGLToolsPckgs/AutoDockTools/Utilities24/:${PATH}"
ENV PYTHONPATH="$PYTHONPATH:/usr/local/MGLToolsPckgs"

# Extra installations
RUN /opt/conda/envs/mgltools/bin/pip install numpy
