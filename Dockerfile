# Using a Conda-ready Docker image
FROM continuumio/miniconda3:latest

RUN mkdir -p /app/Orrery
WORKDIR /app/Orrery

# Create the environment
COPY ./environment.yml ./environment.yml
RUN conda env create -f environment.yml

# Activate the environment
RUN echo "source activate Orrery" > ~/.bashrc
ENV PATH /opt/conda/envs/Orrery/bin:$PATH

# Extra stuff
RUN pip install pybel