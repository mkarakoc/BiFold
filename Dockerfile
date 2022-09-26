# https://hub.docker.com/r/hesap/
FROM hesap/aimpy:main202002270048

MAINTAINER Mesut Karakoç <mesudkarakoc@gmail.com>

# root user
USER root

# add jupyter notebooks
ADD examples    /home/main/examples
RUN chown -R main:main /home/main/examples

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        python3-urllib3=1.25.8-2ubuntu0.1 \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install -U --no-cache-dir \
    bifold==0.731.32 \
    matplotlib==3.6.0 \
    numba==0.56.2 \
    scipy==1.9.1

# main user
USER main

# working directory
WORKDIR /home/main

# jupyter notebook and its ip to use from outside of the Docker container
CMD echo '**************************************************'                          && \
    echo 'check whether the following IP exist if not change the IP in the Dockerfile' && \
    echo http://`hostname -i | awk '{print $1}'`:8080/?token=[COPY PASTE TOKEN HERE]   && \
    echo '**************************************************'                          && \
    echo                                                                               && \     
    jupyter notebook --no-browser --ip=172.17.0.2 --port=8080 

