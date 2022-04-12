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


:contributions: Want to add a convertor ? Please join https://github.com/bioconvert/bioconvert/issues/1
:issues: Please use https://github.com/bioconvert/bioconvert/issues

Overview
########


Life science uses many different formats. They may be old, or with complex syntax and converting those formats may be a challenge. Bioconvert aims at providing a common tool / interface to convert life science data formats from one to another.

Many conversion tools already exist but they may be dispersed, focused on few specific formats, difficult to install, or not optimised. With Bioconvert, we plan to cover a wide spectrum of format conversions; we will re-use existing tools when possible and provide facilities to compare different conversion tools or methods via benchmarking. New implementations are provided when considered better than existing ones.

In March 2022, we had 48 formats, 98 direct conversions (125 different methods).

.. image:: https://raw.githubusercontent.com/bioconvert/bioconvert/master/doc/conversion.png
    :width: 80%


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

From the command line, you can convert a `FastQ` file into
a `FastA` file as follows (compressed or not)::

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


.. list-table:: Conversion table
    :widths: 20 40 40
    :header-rows: 1

    * - Converters
      - CI testing
      - Benchmarking
    * - `abi2fasta <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.abi2fasta>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/abi2fasta.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/abi2fasta.yml
      -
    * - `abi2fastq <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.abi2fastq>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/abi2fastq.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/abi2fastq.yml
      -
    * - `abi2qual <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.abi2qual>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/abi2qual.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/abi2qual.yml
      -
    * - `bam2bedgraph <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.bam2bedgraph>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2bedgraph.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2bedgraph.yml
      -
    * - `bam2bigwig <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.bam2bigwig>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2bigwig.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2bigwig.yml
      -
    * - `bam2cov <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.bam2cov>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2cov.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2cov.yml
      -
    * - `bam2cram <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.bam2cram>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2cram.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2cram.yml
      -
    * - `bam2fasta <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.bam2fasta>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2fasta.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2fasta.yml
      -
    * - `bam2fastq <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.bam2fastq>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2fastq.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2fastq.yml
      -
    * - `bam2json <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.bam2json>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2json.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2json.yml
      -
    * - `bam2sam <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.bam2sam>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2sam.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2sam.yml
      -
    * - `bam2tsv <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.bam2tsv>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2tsv.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2tsv.yml
      -
    * - `bam2wiggle <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.bam2wiggle>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2wiggle.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2wiggle.yml
      -
    * - `bcf2vcf <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.bcf2vcf>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bcf2vcf.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bcf2vcf.yml
      -
    * - `bcf2wiggle <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.bcf2wiggle>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bcf2wiggle.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bcf2wiggle.yml
      -
    * - `bed2wiggle <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.bed2wiggle>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bed2wiggle.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bed2wiggle.yml
      -
    * - `bedgraph2bigwig <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.bedgraph2bigwig>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bedgraph2bigwig.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bedgraph2bigwig.yml
      -
    * - `bedgraph2cov <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.bedgraph2cov>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bedgraph2cov.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bedgraph2cov.yml
      -
    * - `bedgraph2wiggle <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.bedgraph2wiggle>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bedgraph2wiggle.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bedgraph2wiggle.yml
      -
    * - `bigbed2bed <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.bigbed2bed>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bigbed2bed.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bigbed2bed.yml
      -
    * - `bigbed2wiggle <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.bigbed2wiggle>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bigbed2wiggle.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bigbed2wiggle.yml
      -
    * - `bigwig2bedgraph <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.bigwig2bedgraph>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bigwig2bedgraph.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bigwig2bedgraph.yml
      -
    * - `bigwig2wiggle <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.bigwig2wiggle>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bigwig2wiggle.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bigwig2wiggle.yml
      -
    * - `bplink2plink <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.bplink2plink>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bplink2plink.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bplink2plink.yml
      -
    * - `bplink2vcf <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.bplink2vcf>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bplink2vcf.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bplink2vcf.yml
      -
    * - `bz22gz <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.bz22gz>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bz22gz.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bz22gz.yml
      -
    * - `clustal2fasta <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.clustal2fasta>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/clustal2fasta.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/clustal2fasta.yml
      -
    * - `clustal2nexus <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.clustal2nexus>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/clustal2nexus.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/clustal2nexus.yml
      -
    * - `clustal2phylip <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.clustal2phylip>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/clustal2phylip.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/clustal2phylip.yml
      -
    * - `clustal2stockholm <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.clustal2stockholm>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/clustal2stockholm.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/clustal2stockholm.yml
      -
    * - `cram2bam <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.cram2bam>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/cram2bam.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/cram2bam.yml
      -
    * - `cram2fasta <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.cram2fasta>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/cram2fasta.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/cram2fasta.yml
      -
    * - `cram2fastq <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.cram2fastq>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/cram2fastq.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/cram2fastq.yml
      -
    * - `cram2sam <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.cram2sam>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/cram2sam.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/cram2sam.yml
      -
    * - `csv2tsv <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.csv2tsv>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/csv2tsv.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/csv2tsv.yml
      -
    * - `csv2xls <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.csv2xls>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/csv2xls.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/csv2xls.yml
      -
    * - `dsrc2gz <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.dsrc2gz>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/dsrc2gz.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/dsrc2gz.yml
      -
    * - `embl2fasta <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.embl2fasta>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/embl2fasta.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/embl2fasta.yml
      -
    * - `embl2genbank <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.embl2genbank>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/embl2genbank.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/embl2genbank.yml
      -
    * - `fasta2clustal <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.fasta2clustal>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2clustal.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2clustal.yml
      -
    * - `fasta2faa <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.fasta2faa>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2faa.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2faa.yml
      -
    * - `fasta2fasta_agp <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.fasta2fasta_agp>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2fasta_agp.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2fasta_agp.yml
      -
    * - `fasta2fastq <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.fasta2fastq>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2fastq.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2fastq.yml
      -
    * - `fasta2genbank <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.fasta2genbank>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2genbank.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2genbank.yml
      -
    * - `fasta2nexus <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.fasta2nexus>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2nexus.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2nexus.yml
      -
    * - `fasta2phylip <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.fasta2phylip>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2phylip.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2phylip.yml
      -
    * - `fasta2twobit <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.fasta2twobit>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2twobit.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2twobit.yml
      -
    * - `fasta_qual2fastq <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.fasta_qual2fastq>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta_qual2fastq.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta_qual2fastq.yml
      -
    * - `fastq2fasta <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.fastq2fasta>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fastq2fasta.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fastq2fasta.yml
      -  `available <_static/benchmark_fastq2fasta.png>`_
    * - `fastq2fasta_qual <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.fastq2fasta_qual>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fastq2fasta_qual.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fastq2fasta_qual.yml
      -
    * - `fastq2qual <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.fastq2qual>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fastq2qual.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fastq2qual.yml
      -
    * - `genbank2embl <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.genbank2embl>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/genbank2embl.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/genbank2embl.yml
      -
    * - `genbank2fasta <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.genbank2fasta>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/genbank2fasta.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/genbank2fasta.yml
      -
    * - `genbank2gff3 <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.genbank2gff3>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/genbank2gff3.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/genbank2gff3.yml
      -
    * - `gfa2fasta <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.gfa2fasta>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/gfa2fasta.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/gfa2fasta.yml
      -
    * - `gff22gff3 <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.gff22gff3>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/gff22gff3.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/gff22gff3.yml
      -
    * - `gff32gff2 <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.gff32gff2>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/gff32gff2.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/gff32gff2.yml
      -
    * - `gz2bz2 <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.gz2bz2>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/gz2bz2.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/gz2bz2.yml
      -
    * - `gz2dsrc <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.gz2dsrc>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/gz2dsrc.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/gz2dsrc.yml
      -
    * - `json2yaml <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.json2yaml>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/json2yaml.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/json2yaml.yml
      -
    * - `maf2sam <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.maf2sam>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/maf2sam.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/maf2sam.yml
      -
    * - `newick2nexus <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.newick2nexus>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/newick2nexus.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/newick2nexus.yml
      -
    * - `newick2phyloxml <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.newick2phyloxml>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/newick2phyloxml.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/newick2phyloxml.yml
      -
    * - `nexus2clustal <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.nexus2clustal>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/nexus2clustal.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/nexus2clustal.yml
      -
    * - `nexus2fasta <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.nexus2fasta>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/nexus2fasta.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/nexus2fasta.yml
      -
    * - `nexus2newick <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.nexus2newick>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/nexus2newick.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/nexus2newick.yml
      -
    * - `nexus2phylip <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.nexus2phylip>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/nexus2phylip.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/nexus2phylip.yml
      -
    * - `nexus2phyloxml <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.nexus2phyloxml>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/nexus2phyloxml.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/nexus2phyloxml.yml
      -
    * - `ods2csv <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.ods2csv>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/ods2csv.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/ods2csv.yml
      -
    * - `phylip2clustal <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.phylip2clustal>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2clustal.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2clustal.yml
      -
    * - `phylip2fasta <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.phylip2fasta>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2fasta.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2fasta.yml
      -
    * - `phylip2nexus <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.phylip2nexus>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2nexus.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2nexus.yml
      -
    * - `phylip2stockholm <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.phylip2stockholm>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2stockholm.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2stockholm.yml
      -
    * - `phylip2xmfa <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.phylip2xmfa>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2xmfa.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2xmfa.yml
      -
    * - `phyloxml2newick <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.phyloxml2newick>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/phyloxml2newick.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/phyloxml2newick.yml
      -
    * - `phyloxml2nexus <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.phyloxml2nexus>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/phyloxml2nexus.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/phyloxml2nexus.yml
      -
    * - `plink2bplink <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.plink2bplink>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/plink2bplink.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/plink2bplink.yml
      -
    * - `plink2vcf <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.plink2vcf>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/plink2vcf.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/plink2vcf.yml
      -
    * - `sam2bam <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.sam2bam>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/sam2bam.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/sam2bam.yml
      -
    * - `sam2cram <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.sam2cram>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/sam2cram.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/sam2cram.yml
      -
    * - `sam2paf <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.sam2paf>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/sam2paf.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/sam2paf.yml
      -
    * - `scf2fasta <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.scf2fasta>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/scf2fasta.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/scf2fasta.yml
      -
    * - `scf2fastq <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.scf2fastq>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/scf2fastq.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/scf2fastq.yml
      -
    * - `sra2fastq <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.sra2fastq>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/sra2fastq.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/sra2fastq.yml
      -
    * - `stockholm2clustal <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.stockholm2clustal>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/stockholm2clustal.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/stockholm2clustal.yml
      -
    * - `stockholm2phylip <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.stockholm2phylip>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/stockholm2phylip.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/stockholm2phylip.yml
      -
    * - `tsv2csv <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.tsv2csv>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/tsv2csv.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/tsv2csv.yml
      -
    * - `twobit2fasta <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.twobit2fasta>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/twobit2fasta.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/twobit2fasta.yml
      -
    * - `vcf2bcf <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.vcf2bcf>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/vcf2bcf.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/vcf2bcf.yml
      -
    * - `vcf2bed <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.vcf2bed>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/vcf2bed.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/vcf2bed.yml
      -
    * - `vcf2bplink <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.vcf2bplink>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/vcf2bplink.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/vcf2bplink.yml
      -
    * - `vcf2plink <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.vcf2plink>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/vcf2plink.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/vcf2plink.yml
      -
    * - `vcf2wiggle <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.vcf2wiggle>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/vcf2wiggle.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/vcf2wiggle.yml
      -
    * - `wig2bed <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.wig2bed>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/wig2bed.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/wig2bed.yml
      -
    * - `xls2csv <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.xls2csv>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/xls2csv.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/xls2csv.yml
      -
    * - `xlsx2csv <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.xlsx2csv>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/xlsx2csv.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/xlsx2csv.yml
      -
    * - `xmfa2phylip <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.xmfa2phylip>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/xmfa2phylip.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/xmfa2phylip.yml
      -
    * - `yaml2json <https://bioconvert.readthedocs.io/en/master/ref_converters.html#module-bioconvert.yaml2json>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/yaml2json.yml/badge.svg?branch=master
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/yaml2json.yml
      -



Contributors
############

Setting up and maintaining Bioconvert has been possible thanks to users and contributors.
Thanks to all:

.. image:: https://contrib.rocks/image?repo=bioconvert/bioconvert
    :target: https://github.com/bioconvert/bioconvert/graphs/contributors


Changes
########

========= ==============================================================================
Version   Description
========= ==============================================================================
0.5.0     * Add CI actions for all converters
          * remove sniffer (now in biosniff on pypi https://pypi.org/project/biosniff/)
          * A complete benchmarking suite (see doc/Snakefile_benchmark file and
            `benchmarking`)
          * documentation and tests for all converters
          * removed the validators (we assume intputs are correct)
========= ==============================================================================

