FROM continuumio/miniconda3:latest

RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    git

RUN conda install -c conda-forge openbabel

RUN mkdir -p /app/Orrery
WORKDIR /app/Orrery

RUN conda install -c anaconda flask
RUN pip install flask_restful
