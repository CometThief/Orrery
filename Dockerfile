FROM ubuntu:latest

RUN mkdir -p /app/Orrery
WORKDIR /app/Orrery

RUN apt-get update && \
    apt-get install -y python3 python3-rdkit python3-pip

RUN pip install --upgrade pip \
    pandas \
    requests \
    pymongo  \
    tqdm

