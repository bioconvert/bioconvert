User Guide
============

Overview
------------

You can use **bioconvert** from a developer point of view, or as an end-user.
We provide one standalone called::

    bioconvert

You can obtain help by using::

    bioconvert --help

To convert a format into another you need to provide the name of the conversion.
For instance, to convert a fastq file into a fasta, you need to use the sub
command **fastq2fasta**::

    bioconvert fastq2fasta  input.fastq output.fasta

The rationale behind the subcommand choice is manyfold. First, you may have dedicated help
for a given conversion::

    bioconvert fastq2fasta --help

Second, the extensions of your input and output may be non-standard or different
from the bioconvert choice. So, using the subcommand you can do::

    bioconvert fastq2fasta  input.fq output.fa

where the extensions can actually be whatever you want.


Installation
-------------

pip pr conda methods
~~~~~~~~~~~~~~~~~~~~~~~~

For users, **bioconvert** standalone is installed with the package **bioconvert** available on Pypi so you could type::

    pip install bioconvert

This method installs bioconvert and its Python dependencies (available on Pypi website). Note, however, that **bioconvert** may use (depending on the conversion you want to use) external dependencies not available on Pypi. You will need to install those third-party dependencies yourself. An alternative is to install bioconvert using **conda** using::

    conda install bioconvert

Note that you will need to set up the **bioconda** channel (see below for
details).

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
    
You can then create an alias::

    alias bioconvert="singularity run bioconvert.simg bioconvert"




Parallelization
--------------------


Some converters can use several threads, but if you have hundreds of files and
wish to use several CPUs, or run bioconvert on a cluster, we provide here below
a simple Snakefile (snakemake) that can be run easily as follows.

You can download the following file :download:`Snakefile` 

.. literalinclude:: Snakefile
    :language: python
    :lines: 12-28

and execute it locally as follows (assuming you have 4 CPUs)::

    snakemake -s Snakefile --cores 4

or on a cluster::

    snakemake -s Snakefile --cluster "--mem=1000 -j 10"








