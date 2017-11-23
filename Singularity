#Bootstrap: debootstrap
#OSVersion: zesty
#MirrorURL: http://archive.ubuntu.com/ubuntu/


BootStrap: docker
From: ubuntu:17.04


%labels
  maintainer Bertrand Neron <bneron@pasteur.fr>
  package.name
  package.version 0.0.7
  package.homepage https://pypi.python.org/pypi/bioconvert/0.0.6
  package.source.url
  package.source.mdm5
  package.license GPLv3

%post
  ######### install system #########

  echo "deb http://archive.ubuntu.com/ubuntu/ zesty main restricted" > /etc/apt/sources.list
  echo "deb http://archive.ubuntu.com/ubuntu/ zesty-updates main restricted" >> /etc/apt/sources.list

  echo "deb http://archive.ubuntu.com/ubuntu/ zesty universe" >> /etc/apt/sources.list
  echo "deb-src http://archive.ubuntu.com/ubuntu/ zesty universe" >> /etc/apt/sources.list
  echo "deb http://archive.ubuntu.com/ubuntu/ zesty-updates universe" >> /etc/apt/sources.list
  echo "deb-src http://archive.ubuntu.com/ubuntu/ zesty-updates universe" >> /etc/apt/sources.list

  echo "deb http://archive.ubuntu.com/ubuntu/ zesty multiverse" >> /etc/apt/sources.list
  echo "deb http://archive.ubuntu.com/ubuntu/ zesty-updates multiverse" >> /etc/apt/sources.list

  echo "deb http://archive.ubuntu.com/ubuntu/ zesty-backports main restricted universe multiverse" >> /etc/apt/sources.list

  echo "deb http://security.ubuntu.com/ubuntu/ zesty-security main restricted" >> /etc/apt/sources.list
  echo "deb http://security.ubuntu.com/ubuntu/ zesty-security universe" >> /etc/apt/sources.list
  echo "deb-src http://security.ubuntu.com/ubuntu/ zesty-security universe" >> /etc/apt/sources.list
  echo "deb http://security.ubuntu.com/ubuntu/ zesty-security multiverse" >> /etc/apt/sources.list

  apt-get update -y
  apt-get install -y wget bzip2
  apt-get install -y libgl1-mesa-glx

  # install anaconda
  if [ ! -d /usr/local/anaconda ]; then
      # wget https://repo.continuum.io/miniconda/Miniconda3-4.3.14-Linux-x86_64.sh\
      # for now, we use 4.2.12 to have python3.5 by default so no need to
      # create a new env saving space in the process. The reason for using 3.5
      # is inherent to the packages used at the moment.
      wget https://repo.continuum.io/miniconda/Miniconda3-4.2.12-Linux-x86_64.sh\
           -O ~/anaconda.sh && \
      bash ~/anaconda.sh -b -p /usr/local/anaconda && \
      rm ~/anaconda.sh
  fi

  # set anaconda path
  export PATH=$PATH:/usr/local/anaconda/bin
  conda update conda

  conda config --add channels r
  conda config --add channels defaults
  conda config --add channels conda-forge
  conda config --add channels bioconda

  # The main packages for sequana:
  conda install --file https://raw.githubusercontent.com/biokit/bioconvert/master/requirements.txt
  conda install --file https://raw.githubusercontent.com/biokit/bioconvert/master/requirements_tools.txt
  #conda install -y graphviz==2.38 pygraphviz

  

  ######### install bioconvert #########
  pip install bioconvert==0.0.7

  # add this directory for Institut Pasteur cluster usage
  if [ ! -d /pasteur ]; then mkdir /pasteur; fi

  # Uses agg as backend instead of qt (less dependencies)
  echo "backend:tkagg" > matplotlibrc

  ######## clean image ########
  apt-get autoremove -y
  apt-get clean -y
  conda clean -y --all

%environment
  export PATH=$PATH:/usr/local/anaconda/bin


%runscript
  exec /usr/local/anaconda/bin/bioconvert "$@"


