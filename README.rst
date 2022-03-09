Bioconvert
##########

**Bioconvert** is a collaborative project to facilitate the interconversion of life science data from one format to another.

.. image:: https://badge.fury.io/py/bioconvert.svg
    :target: https://pypi.python.org/pypi/bioconvert

.. image:: https://github.com/bioconvert/bioconvert/actions/workflows/main.yml/badge.svg?branch=master
    :target: https://github.com/bioconvert/bioconvert/actions/workflows/main.yml

.. image:: https://github.com/bioconvert/bioconvert/actions/workflows/main.yml/badge.svg?branch=refactoring
    :target: https://github.com/bioconvert/bioconvert/actions/workflows/main.yml

.. image:: https://coveralls.io/repos/github/bioconvert/bioconvert/badge.svg?branch=master
   :target: https://coveralls.io/github/bioconvert/bioconvert?branch=master

.. image:: http://readthedocs.org/projects/bioconvert/badge/?version=master
    :target: http://bioconvert.readthedocs.org/en/master/?badge=master
    :alt: Documentation Status

.. image::  https://img.shields.io/github/issues/bioconvert/bioconvert.svg
    :target:  https://github.com/bioconvert/bioconvert/issues

.. image:: https://anaconda.org/bioconda/bioconvert/badges/platforms.svg
   :target: https://anaconda.org/bioconda/bioconvert

.. image::  https://anaconda.org/bioconda/bioconvert/badges/installer/conda.svg
    :target: https://conda.anaconda.org/bioconda

.. image:: https://mybinder.org/badge_logo.svg
    :target: https://mybinder.org/v2/gh/bioconvert/bioconvert/master


:contributions: Want to add a convertor ? Please join https://github.com/bioconvert/bioconvert/issues/1
:issues: Please use https://github.com/bioconvert/bioconvert/issues

Overview
########


Life science uses many different formats. They may be old, or with complex syntax and converting those formats may be a challenge. Bioconvert aims at providing a common tool / interface to convert life science data formats from one to another.

Many conversion tools already exist but they may be dispersed, focused on few specific formats, difficult to install, or not optimised. With Bioconvert, we plan to cover a wide spectrum of format conversions; we will re-use existing tools when possible and provide facilities to compare different conversion tools or methods via benchmarking. New implementations are provided when considered better than existing ones.

In Aug 2019, we had 46 formats, 98 direct conversions (156 different methods). More conversions are possible when calling bioconvert several times.

In Aug 2018, we had 43 formats, 79 direct conversions (129 different methods). More conversions are possible when calling bioconvert several times.

In June 2018, we had 66 direct conversions (120 different methods). More conversions are possible when calling bioconvert several times.

Installation
###############

In order to install bioconvert, you can use **pip**::

    pip install bioconvert

We also provide releases on bioconda (http://bioconda.github.io/)::

    conda install bioconvert

and Singularity container are available. See
http://bioconvert.readthedocs.io/en/master/user_guide.html#installation for
details.

Usage
##########

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
    bioconvert fastq2fasta --help


From a Python shell::

    # import a converter
    from bioconvert.fastq2fasta import FASTQ2FASTA

    # Instanciate with infile/outfile names
    convert = FASTQ2FASTA(infile, outfile)

    # the conversion itself
    convert()




Available Converters
#######################

.. image:: https://raw.githubusercontent.com/bioconvert/bioconvert/master/doc/conversion.png
    :width: 80%

.. list-table:: Conversion table
    :widths: 20 40 40
    :header-rows: 1

    * - Converters
      - CI testing
      - Benchmarking
    * - fastq2fasta
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fastq2fasta.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fastq2fasta.yml
      - link IMG benchmarking fastq2fasta
    * - bam2sam
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2sam.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2sam.yml
      - link IMG benchmarking bam2sam
    * - abi2fasta
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/abi2fasta.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/abi2fasta.yml
      - link IMG benchmarking abi2fasta
    * - fastq2qual
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fastq2qual.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fastq2qual.yml
      - link IMG benchmarking fastq2qual
    * - fasta2fastq
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2fastq.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2fastq.yml
      - link IMG benchmarking fasta2fastq
    * - sam2bam
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/sam2bam.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/sam2bam.yml
      - link IMG benchmarking sam2bam
    * - fasta2twobit
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2twobit.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2twobit.yml
      - link IMG benchmarking fasta2twobit
    * - twobit2fasta
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/twobit2fasta.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/twobit2fasta.yml
      - link IMG benchmarking twobit2fasta
    * - fasta_qual2fastq
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta_qual2fastq.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta_qual2fastq.yml
      - link IMG benchmarking fasta_qual2fastq
    * - fastq2fasta_qual
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fastq2fasta_qual.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fastq2fasta_qual.yml
      - link IMG benchmarking fastq2fasta_qual
    * - bam2cram
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2cram.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2cram.yml
      - link IMG benchmarking bam2cram
    * - cram2sam
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/cram2sam.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/cram2sam.yml
      - link IMG benchmarking cram2sam
    * - sam2cram
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/sam2cram.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/sam2cram.yml
      - link IMG benchmarking sam2cram
    * - abi2fastq
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/abi2fastq.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/abi2fastq.yml
      - link IMG benchmarking abi2fastq
    * - sra2fastq
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/sra2fastq.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/sra2fastq.yml
      - link IMG benchmarking sra2fastq
    * - fasta2faa
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2faa.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2faa.yml
      - link IMG benchmarking fasta2faa
    * - gz2bz2
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/gz2bz2.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/gz2bz2.yml
      - link IMG benchmarking gz2bz2
    * - bz22gz
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bz22gz.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bz22gz.yml
      - link IMG benchmarking bz22gz
    * - gz2dsrc
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/gz2dsrc.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/gz2dsrc.yml
      - link IMG benchmarking gz2dsrc
    * - dsrc2gz
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/dsrc2gz.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/dsrc2gz.yml
      - link IMG benchmarking dsrc2gz


Contributors
############

Setting up and maintaining Bioconvert has been possible thanks to users and contributors. 
Thanks to all:

.. image:: https://contrib.rocks/image?repo=bioconvert/bioconvert
    :target: https://github.com/bioconvert/bioconvert/graphs/contributors