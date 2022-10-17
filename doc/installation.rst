
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

You can now use the environment.yaml file provided in the github repository to create a new environment called
**bioconvert**::

    conda install create -f environment.yml

Docker
------

A Dockerfile (version 0.6.2 of **BioConvert**) is available on dockerhub::

    docker pull bioconvert/bioconvert:0.6.2

Which can be used as follows::

    docker run bioconvert -d /home/user:/home/user bioconvert /home/user/test_file.fastq /home/user/test_file.fasta

Singularity/Apptainer
----------------------

We provide Singulariry image of **BioConvert** within the https://damona.readthedocs.io project.

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

.. warning:: You will need singularity of course. If you have a conda envirnment, you are lucky. singularity is there/ 

