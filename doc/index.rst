Bioconvert
##########

**Bioconvert** is a collaborative project to facilitate the interconversion of life science data from one format to another. **Bioconvert** currently contains more than 40 formats, 90 conversions.

.. image:: https://badge.fury.io/py/bioconvert.svg
    :target: https://pypi.python.org/pypi/bioconvert

.. image:: https://img.shields.io/pypi/pyversions/bioconvert.svg
   :target: https://www.python.org

.. image:: https://secure.travis-ci.org/bioconvert/bioconvert.png
    :target: http://travis-ci.org/bioconvert/bioconvert

.. image:: https://coveralls.io/repos/github/bioconvert/bioconvert/badge.svg?branch=master
   :target: https://coveralls.io/github/bioconvert/bioconvert?branch=master

.. image:: http://readthedocs.org/projects/bioconvert/badge/?version=master
    :target: http://bioconvert.readthedocs.org/en/master/?badge=master
    :alt: Documentation Status

.. .. image:: https://badges.gitter.im/bioconvert/bioconvert.svg
    :target: https://gitter.im/bioconvert/Lobby?source=orgpage

.. image::  https://img.shields.io/github/issues/bioconvert/bioconvert.svg
    :target:  https://github.com/bioconvert/bioconvert/issues

.. image:: https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg
   :target: https://singularity-hub.org/collections/135

.. image:: https://anaconda.org/bioconda/bioconvert/badges/platforms.svg
   :target: https://anaconda.org/bioconda/bioconvert

.. image::  https://anaconda.org/bioconda/bioconvert/badges/installer/conda.svg
    :target: https://conda.anaconda.org/bioconda

:contributions: Please join the team to contribute: https://github.com/bioconvert/bioconvert/issues/1

.. image:: conversion.png
   :width: 80%

Overview
########

Life science uses many different formats. They may be old, or with complex
syntax and converting those formats may be a challenge. **Bioconvert** aims at providing a common tool / interface to convert life science data formats from one to another.

Many conversion tools already exist but they may be dispersed, focused on few
specific formats, difficult to install, or not optimised. With **Bioconvert**, we plan to
cover a wide spectrum of format conversions; we will re-use existing tools when
possible and provide facilities to compare different conversion tools or methods 
via :ref:`benchmarking <benchmarking>`. New implementations are provided when considered 
better than existing ones.



**In Aug 2018, we had 43 formats, 79 direct conversions (129 different methods). More conversions are possible when calling bioconvert several times.**

**In June 2018, we had 66 direct conversions (120 different methods). More conversions are possible when calling bioconvert several times.**


Installation
############

In order to install bioconvert, you can use **pip**::

    pip install bioconvert

This command installs bioconvert and its Python dependencies. Note, however,
that bioconvert may need extra non-Python dependencies as indicated in this
`requirements file <https://raw.githubusercontent.com/bioconvert/bioconvert/master/requirements_tools.txt>`_.

Since Jan 2018 we also provide some versions on bioconda. If you already have
bioconda setup on your system, just type::

    conda install bioconvert

Or if you have never done so, please add those channels before hand (provided
you have installed conda)::

    conda config --add channels r
    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda

Otherwise, please see the instructions in the :ref:`installation` section where
you can find information about our Singularity container as well.


Usage
#####

From the command line, you can convert a :term:`FastQ` file into 
a :term:`FastA` file as follows (compressed or not)::

    bioconvert fastq2fasta input.fastq output.fasta
    bioconvert fastq2fasta input.fq    output.fasta
    bioconvert fastq2fasta input.fq.gz output.fasta.gz
    bioconvert fastq2fasta input.fq.gz output.fasta.bz2

When there is no ambiguity, you can be implicit::    

    bioconvert input.fastq output.fasta

For help, just type::

    bioconvert --help

From Python shell::

    # import a converter
    from bioconvert.fastq2fasta import FASTQ2FASTA

    # Instanciate with infile/outfile names
    convert = FASTQ2FASTA(infile, outfile)

    # the conversion itself
    convert()



User and Developer Guides
#########################

.. toctree::
    :maxdepth: 2
    :numbered:

    installation
    user_guide
    tutorial
    developer_guide
    benchmarking
    auto_examples/index
    references
    formats
    not_provided
    faqs
    glossary
    ChangeLog.rst




