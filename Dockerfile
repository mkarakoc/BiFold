# hash:sha256:29ad7a02ec02f0b20ad111e2c7a43b2ed8924979053cd14b7716ebcc43385e11
FROM registry.codeocean.com/codeocean/miniconda3:4.9.2-cuda11.7.0-cudnn8-ubuntu20.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        python3-urllib3=1.25.8-2ubuntu0.1 \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install -U --no-cache-dir \
    bifold==0.731.32 \
    html-parser==0.2 \
    matplotlib==3.6.0 \
    numba==0.56.2 \
    numpy==1.23.3 \
    scipy==1.9.1
