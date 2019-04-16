Bioconvert
==========

Bioconvert is a project to facilitate the interconversion of life science data from one format to another.

.. image:: https://badge.fury.io/py/bioconvert.svg
    :target: https://pypi.python.org/pypi/bioconvert

.. image:: https://secure.travis-ci.org/bioconvert/bioconvert.png
    :target: http://travis-ci.org/bioconvert/bioconvert

.. image:: https://coveralls.io/repos/github/bioconvert/bioconvert/badge.svg?branch=master
    :target: https://coveralls.io/github/bioconvert/bioconvert?branch=master

.. image:: http://readthedocs.org/projects/bioconvert/badge/?version=master
    :target: http://bioconvert.readthedocs.org/en/latest/?badge=master
    :alt: Documentation Status

.. image:: https://badges.gitter.im/bioconvert/bioconvert.svg
    :target: https://gitter.im/bioconvert/Lobby?source=orgpage


:note: Bioconvert is tested with Travis for the following Python version: 3.5 and 3.6. Python 2 won't be provided.

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


