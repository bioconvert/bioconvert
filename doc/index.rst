


Bioconvert
####################

Bioconvert is a project to facilitate the interconversion of life science data from one format to another.


.. image:: https://badge.fury.io/py/bioconvert.svg
    :target: https://pypi.python.org/pypi/bioconvert

.. image:: https://secure.travis-ci.org/biokit/bioconvert.png
    :target: http://travis-ci.org/biokit/bioconvert

.. image:: https://coveralls.io/repos/github/biokit/bioconvert/badge.svg?branch=master
   :target: https://coveralls.io/github/biokit/bioconvert?branch=master

.. image:: http://readthedocs.org/projects/bioconvert/badge/?version=master
    :target: http://bioconvert.readthedocs.org/en/master/?badge=master
    :alt: Documentation Status

.. image:: https://badges.gitter.im/biokit/bioconvert.svg
    :target: https://gitter.im/bioconvert/Lobby?source=orgpage


:note: Bioconvert is tested with Travis for the following Python version: 3.5, 3.6
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

From the command line::

    bioconvert input.fastq output.fasta
    bioconvert input.fq output.fasta
    bioconvert input.mybam output.bed --input-format bam
    bioconvert --formats
    bioconvert --help

From Python shell::

    # import a converter
    from bioconvert.fastq2fasta import Fastq2Fasta

    # Instanciate with infile/outfile names
    convert = Fastq2Fasta(infile, outfile)

    # the conversion itself
    convert()

Available conversion:
---------------------------

.. image:: conversion.png
   :width: 80%

User and Developer Guides
#############################

.. toctree::
    :maxdepth: 2
    :numbered:

    user_guide
    developer_guide
    benchmarking
    auto_examples/index
    references
    formats
    faqs
    glossary
    ChangeLog.rst




