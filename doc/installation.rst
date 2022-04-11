
.. _installation:

Installation
============

**Bioconvert** is developed in Python so you can use the **pip** method to install it easily.
Note, however, that you will lack lots of functionalities since **Bioconvert** relies on third-party software.

We therefore recommend to also use **bioconda** to install those dependencies.


pip installation
----------------

For users, **Bioconvert** can be installed with the **pip** tool as follows::

    pip install bioconvert

This method installs **Bioconvert** and its Python dependencies. Note, however, that **bioconvert** may use (depending on the conversion you want to use) external dependencies not available on Pypi. You will need to install those third-party dependencies yourself. An alternative is to install bioconvert using **conda** as explained here after


conda / bioconda installation
-----------------------------

You can install **Bioconvert** and its dependencies with **conda** as follows::

    conda install bioconvert

Note that you will need to set up the **bioconda** channel (see below for
details).
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


From source
-----------

For developers, if you want to use the latest source files and based on
**bioconda** to install the dependencies you just need to do::

    git clone https://github.com/bioconvert/bioconvert
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

For other version, or to install singularity on windows or Mac, please check out the singularity website `<http://singularity.lbl.gov/>`_

First, download the container. For the latest version::

    singularity pull --name bioconvert.img shub://bioconvert/bioconvert:latest

or for a specific version::

    singularity pull --name bioconvert.img shub://bioconvert/bioconvert:0_3_0

You can then create an alias::

    alias bioconvert="singularity run bioconvert.simg bioconvert"
