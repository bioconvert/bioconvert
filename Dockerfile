FROM ubuntu:19.04

RUN apt-get update && \
    apt-get install -y && \
    apt-get install -y wget bzip2 && \
    apt-get install -y libglu1-mesa 

# gzip is not needed to install bioconvert 
# but we need it to install webapp
# in webapp Dockerfile we running under bioconvert id
RUN apt-get install -y unzip \ 
 && rm -rf /var/lib/apt/lists/*

RUN useradd bioconvert -g users -m 
USER bioconvert

ENV PATH /home/bioconvert/miniconda3/bin:$PATH

WORKDIR /home/bioconvert

RUN wget https://repo.continuum.io/miniconda/Miniconda3-4.3.14-Linux-x86_64.sh -O miniconda.sh && \
    chmod +x miniconda.sh && \
    bash miniconda.sh -b -p $HOME/miniconda3 && \
    export PATH=$HOME/miniconda3/bin:$PATH && \
    hash -r && \
    conda update --yes conda && \
    conda config --add channels r && \
    conda config --add channels defaults && \
    conda config --add channels conda-forge && \
    conda config --add channels bioconda

ADD requirements.txt requirements.txt
ADD requirements_tools.txt requirements_tools.txt

RUN conda install --yes --file requirements.txt && \
    conda install --yes --file requirements_tools.txt && \
    pip install pypandoc && \
    pip install --upgrade pip && \
    pip install bioconvert==0.4.3 
