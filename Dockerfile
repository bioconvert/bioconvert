FROM ubuntu:21.04

RUN apt-get update -y && apt-get -y install wget libgl1-mesa-glx libegl1-mesa libxrandr2 libxrandr2 libxss1 libxcursor1 libxcomposite1 libasound2 libxi6 libxtst6 libcurl4-openssl-dev build-essential git zlib1g zlib1g-dev

# gzip is not needed to install bioconvert
# but we need it to install webapp
# in webapp Dockerfile we running under bioconvert id
RUN apt-get install -y unzip \
 && rm -rf /var/lib/apt/lists/*

RUN useradd bioconvert -g users -m
USER bioconvert

ENV PATH /home/bioconvert/miniconda3/bin:$PATH

WORKDIR /home/bioconvert

RUN wget https://repo.continuum.io/miniconda/Miniconda3-py37_4.11.0-Linux-x86_64.sh -O miniconda.sh && \
    chmod +x miniconda.sh && \
    bash miniconda.sh -b -p $HOME/miniconda3 && \
    export PATH=$HOME/miniconda3/bin:$PATH && \
    hash -r && \
    conda update --yes conda && \
    conda config --add channels r && \
    conda config --add channels conda-forge && \
    conda config --add channels bioconda

#ADD requirements.txt requirements.txt
COPY spec-file.txt spec-file.txt

RUN conda install --yes --file spec-file.txt
RUN conda install --yes pybigwig
RUN pip install --upgrade pip && pip install bioconvert==0.6.3

