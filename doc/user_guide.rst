User Guide
============

Overview
------------

You can use **bioconvert** from a developer point of view, or as an end-user.
We provide one standalone called::

    bioconvert

You can obtain help by using::

    bioconvert --help


Installation
-------------

pip method
~~~~~~~~~~~~~
For developers, **bioconvert** standalone is installed with the package available on Pypi so you could type::

    pip install bioconvert 

This method installs bioconvert and its Python dependencies (available on Pypi website). Note, however, that **bioconvert** may use (depending on the conversion you want to use) external dependencies not available on Pypi. You will need to install those third-party dependencies yourself. An alternative is to install bioconvert using **conda**. 

conda / bioconda method
~~~~~~~~~~~~~~~~~~~~~~~~~

::

    conda config --add channels r
    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda

.. warning:: it is important to add them in this order, as mentionned on bioconda webpage    (https://bioconda.github.io/).

If you have already set the channels, please check that the order is correct.
With the following command::

    conda config --get channels

You should see::

    --add channels 'r'   # lowest priority
    --add channels 'defaults'
    --add channels 'conda-forge'
    --add channels 'bioconda'   # highest priority

Finally, get the source, install the dependencies using conda, and install
bioconvert as follows::

    git clone https://github.com/biokit/bioconvert
    cd bioconvert
    conda install --file requirements.txt
    conda install --file requirements_tools.txt
    conda install --file requirements_dev.txt
    python setup.py install

    


Singularity
------------

For production, we would recommend to use the singularity container.

This method will download a file (container) with everything pre-compiled and
pre-installed with the latest version available on Pypi.

First, you will need to install singularity. To install the version 2.4 of
singularity on a Linux plaform, just download and execute this :download:`install_singularity.sh` bash script, or just type these commands::

    VERSION=2.4
    wget https://github.com/singularityware/singularity/releases/download/$VERSION/singularity-$VERSION.tar.gz
    tar xvf singularity-$VERSION.tar.gz
    cd singularity-$VERSION
    ./configure --prefix=/usr/local
    make
    sudo make install

.. note:: here we need to be sudo, but you can install singularity localy if needed. 

For other version, or to install singularity on windows or Mac, please check out the singularity website singularity `<http://singularity.lbl.gov/>`_

First, download the container::

    singularity pull --name bioconvert.img shub://biokit/bioconvert:latest







