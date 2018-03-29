FROM python:3.6
ENV TRAVIS_PYTHON_VERSION 3.6.4
ENV PATH /root/miniconda3/bin:$PATH

RUN apt-get update && \
    apt-get install -y \
        nano \
        wget \
        libglu1-mesa \
 && rm -rf /var/lib/apt/lists/*

RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    chmod +x miniconda.sh && \
    ./miniconda.sh -b && \
    export PATH=$HOME/miniconda3/bin:$PATH && \
    hash -r && \
    conda update --yes conda && \
    conda config --add channels r && \
    conda config --add channels defaults && \
    conda config --add channels conda-forge && \
    conda config --add channels bioconda
#     && \
#    rm -rf /dev/shm && \
#    ln -s /run/shm /dev/shm && \
#    export DISPLAY=:99.0 && \
#    sh -e /etc/init.d/xvfb start

ADD requirements.txt /requirements.txt
ADD requirements_tools.txt /requirements_tools.txt
ADD requirements_dev.txt /requirements_dev.txt

RUN conda install --yes python=$TRAVIS_PYTHON_VERSION --file requirements.txt && \
    conda install --yes python=$TRAVIS_PYTHON_VERSION --file requirements_tools.txt && \
    conda install --yes python=$TRAVIS_PYTHON_VERSION --file requirements_dev.txt && \
    conda install --yes python=$TRAVIS_PYTHON_VERSION  coverage && \
    pip install pygatb --no-deps && \
    pip install pypandoc && \
    pip install biocode && \
    pip install --upgrade pip

COPY . /code/
WORKDIR /code
RUN cd /code/ && \
    find -type f -name '*.pyc' -delete && \
    find -type f -name '__pycache__' -delete && \
    cd /code/test/ && \
    find -type f -name '__init__.py' -delete

RUN pip install -e .

CMD pytest -v --durations=10  test/ --cov=bioconvert --cov-report term --timeout 300 -n 1
#docker build . -t bioconvert && docker run -it bioconvert