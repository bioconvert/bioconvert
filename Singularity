Bootstrap: debootstrap
OSVersion: zesty
MirrorURL: http://archive.ubuntu.com/ubuntu/


%labels
  maintainer Bertrand Neron <bneron@pasteur.fr>
  package.name 
  package.version 0.0.1
  package.homepage https://pypi.python.org/pypi/bioconvert/0.0.1
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
  apt-get install -y --no-install-recommends python3
  apt-get install -y python3-pip

  ######### install bioconvert ######### 
  pip3 install bioconvert==0.0.1

  ######## clean image ########
  apt-get autoremove -y
  apt-get clean -y

%runscript
  exec /usr/local/bin/bioconvert "$@"

%tests

