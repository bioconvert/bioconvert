Bioconvert
==========

**Bioconvert** is a collaborative project to facilitate the interconversion of life science data from one format to another.

.. image:: https://badge.fury.io/py/bioconvert.svg
    :target: https://pypi.python.org/pypi/bioconvert

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



:contributions: Want to add a convertor ? Please join https://github.com/bioconvert/bioconvert/issues/1
:issues: Please use https://github.com/bioconvert/bioconvert/issues


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

::

    bioconvert fastq2fasta input.fastq output.fasta
    bioconvert gz2dsrc input.fq.gz output.dsrc2
    bioconvert bam2bed input.bam output.bed
    bioconvert --help

Available Converters
#######################

.. image:: https://raw.githubusercontent.com/bioconvert/bioconvert/master/doc/conversion.png
    :width: 80%


