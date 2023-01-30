
.. _installation_details:

Installation
============

**Bioconvert** is developed in Python so you can use the **pip** method to install it easily. We recommend to use a
virtual environment to not interfere with your system. In any case, install **BioConvert** with::

    pip install bioconvert

Note, however, that you will be able to use only about half of the conversions (pure Python). Others depend 
on third-party software.

One solution is to create a dedicated environment using **conda**. In particular, we use `bioconda <https://bioconda.github.io>`_ to install those dependencies.

conda / bioconda /mamba installation
--------------------------------------

One workable and relatively straightfoward installation is based on conda/mamba::


    conda create --name bioconvert python=3.8
    conda activate bioconvert
    conda  install mamba

Then, use mamba to install the missing executable. For example for **samtools**::

    mamba install samtools 

Third-package executables can be installed with your own method. We recommend and provide solutions for **conda**.
Indeed, **BioConvert** is available on the `bioconda <https://bioconda.github.io>`_ channel (see :ref:`conda_channels` section for details).

So, you could create a conda environment and install **bioconvert** directly with all dependencies. This is, however, pretty slow due to the large
number of dependencies::

    conda create --name bioconvert bioconvert

Instead, we recommend to use an intermediate tool called **mamba** that will provide a more robust and faster
installation::

    conda create -c bioconda --name bioconvert mamba
    conda activate bioconvert
    mamba install bioconvert

In Jan 2023, this method worked out of box and created an environment with Python3.10 and bioconvert 0.6.2 with all its
dependencies.

We also provide a frozen version of an environment with the bioconvert github repository. Note, however, that this file
may change with time. This will create a conda environment called bioconvert. See the link 

    wget https://raw.githubusercontent.com/bioconvert/bioconvert/main/environment.yml -O test.yml
    conda create install create -f test.yml

Docker
------

A Dockerfile (version 0.6.1 of **BioConvert**) is available on dockerhub::

    docker pull bioconvert/bioconvert:0.6.1

Which can be used as follows::

    docker run bioconvert -d /home/user:/home/user bioconvert /home/user/test_file.fastq /home/user/test_file.fasta

Since **bioconvert** is on bioconda, it is also available on quay.io. For instance, version 0.6.2 is reachable here::

    docker pull quay.io/biocontainers/bioconvert:0.6.2--pyhdfd78af_0

Singularity/Apptainer
----------------------

We provide Singularity/Apptainer images of **BioConvert** within the https://damona.readthedocs.io project.

The version 0.6.2 of **BioConvert** is available for downloads.

Using damona::

    pip install damona

    # create and activate an environment
    damona env --create test_bioconvert
    damona activate test_bioconvert
    damona install bioconvert
    bioconvert

You can also install the singularity image yourself by downloading it::

    wget https://zenodo.org/record/7034822/files/bioconvert_0.6.1.img
    singularity exec bioconvert_0.6.1.img bioconvert

    # you can also create an alias
    alias bioconvert="singularity run bioconvert.simg bioconvert"

.. warning:: You will need singularity of course. If you have a conda environment, you are lucky. singularity is there/ 

.. _conda_channels:

Conda channels
--------------

First, you will need to set up the **bioconda** channel if not already done::

    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda config --set channel_priority strict

.. warning:: it is important to add them in this order, as mentionned on bioconda webpage    (https://bioconda.github.io/).

If you have already set the channels, please check that the order is correct.
With the following command::

    conda config --get channels

You should see::

    --add channels 'defaults'
    --add channels 'bioconda'
    --add channels 'conda-forge'# highest priority
