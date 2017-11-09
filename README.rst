Bioconvert
==========

Bioconvert is a project to facilitate the interconversion of life science data from one format to another.

.. image:: https://badge.fury.io/py/bioconvert.svg
    :target: https://pypi.python.org/pypi/bioconvert

.. image:: https://secure.travis-ci.org/biokit/bioconvert.png
    :target: http://travis-ci.org/biokit/bioconvert

.. image:: https://coveralls.io/repos/github/biokit/bioconvert/badge.svg?branch=master
    :target: https://coveralls.io/github/biokit/bioconvert?branch=master

.. image:: http://readthedocs.org/projects/bioconvert/badge/?version=master
    :target: http://bioconvert.readthedocs.org/en/latest/?badge=master
    :alt: Documentation Status

.. image:: https://badges.gitter.im/biokit/bioconvert.svg
    :target: https://gitter.im/bioconvert/Lobby?source=orgpage


:note: Bioconvert is tested with Travis for the following Python version: 3.5 and 3.6. Python 2 won't be provided.

:contributions: Please join https://github.com/biokit/bioconvert
:issues: Please use https://github.com/biokit/bioconvert/issues


Installation
###############

In order to install bioconvert, you can use **pip**::

    pip install bioconvert

.. Or using bioconda channel from the Anaconda project::

..    conda install bioconvert

Usage
##########

::

    bioconvert input.fastq output.fasta
    bioconvert input.fq output.fasta
    bioconvert input.mybam output.bed --input-format bam
    bioconvert --formats
    bioconvert --help


