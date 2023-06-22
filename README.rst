Bioconvert
##########

**Bioconvert** is a collaborative project to facilitate the interconversion of life science data from one format to another.

.. image:: https://badge.fury.io/py/bioconvert.svg
    :target: https://pypi.python.org/pypi/bioconvert

.. image:: https://github.com/bioconvert/bioconvert/actions/workflows/main.yml/badge.svg?branch=main
    :target: https://github.com/bioconvert/bioconvert/actions/workflows/main.yml

.. image:: https://coveralls.io/repos/github/bioconvert/bioconvert/badge.svg?branch=main
   :target: https://coveralls.io/github/bioconvert/bioconvert?branch=main

.. image:: http://readthedocs.org/projects/bioconvert/badge/?version=main
    :target: http://bioconvert.readthedocs.org/en/main/?badge=main
    :alt: Documentation Status

.. image::  https://img.shields.io/github/issues/bioconvert/bioconvert.svg
    :target:  https://github.com/bioconvert/bioconvert/issues

.. image:: https://anaconda.org/bioconda/bioconvert/badges/platforms.svg
   :target: https://anaconda.org/bioconda/bioconvert/badges/platforms.svg

.. image::  https://anaconda.org/bioconda/bioconvert/badges/version.svg
    :target: https://anaconda.org/bioconda/bioconvert

.. image:: https://anaconda.org/bioconda/bioconvert/badges/downloads.svg
    :target: https://github.com/bioconvert/bioconvert/releases

.. image:: https://zenodo.org/badge/106598809.svg
   :target: https://zenodo.org/badge/latestdoi/106598809

.. image:: https://static.pepy.tech/personalized-badge/bioconvert?period=month&units=international_system&left_color=black&right_color=blue&left_text=Downloads/months
    :target: https://pepy.tech/project/bioconvert

.. image:: https://raw.githubusercontent.com/bioconvert/bioconvert/main/doc/_static/logo_300x200.png
    :target: https://raw.githubusercontent.com/bioconvert/bioconvert/main/doc/_static/logo_300x200.png


:contributions: Want to add a convertor ? Please join https://github.com/bioconvert/bioconvert/issues/1

Overview
########


Life science uses many different formats. They may be old, or with complex syntax and converting those formats may be a challenge. Bioconvert aims at providing a common tool / interface to convert life science data formats from one to another.

Many conversion tools already exist but they may be dispersed, focused on few specific formats, difficult to install, or not optimised. With Bioconvert, we plan to cover a wide spectrum of format conversions; we will re-use existing tools when possible and provide facilities to compare different conversion tools or methods via benchmarking. New implementations are provided when considered better than existing ones.

In Jan 2023, we had 50 formats, 100 direct conversions available.

.. image:: https://raw.githubusercontent.com/bioconvert/bioconvert/main/doc/conversion.png
    :width: 80%


Installation
###############

**BioConvert** is developped in Python. Please use conda or any Python environment manager to install **BioConvert** using the **pip** command::

    pip install bioconvert

50% of the conversions should work out of the box. However, many conversions require external tools. This is why we
recommend to use a **conda** environment. In particular, most external tools are available on the **bioconda** channel. 
For instance if you want to convert a SAM file to a BAM file you would need to install **samtools** as follow::

    conda install -c bioconda samtools

Since **bioconvert** is available on `bioconda <https://bioconda.github.io>`_ on solution that installs **BioConvert** and all its dependencies is to use conda/mamba::

    conda env create --name bioconvert mamba
    conda activate bioconvert
    mamba install bioconvert
    bioconvert --help

See the Installation section for more details and alternative solutions (docker, singularity).

Quick Start
##############
There are many conversions available. Type::

    bioconvert --help 

to get a list of valid method of conversions. Taking the example of a conversion from a `FastQ` file into
a `FastA` file, you could do the conversion as follows::

    bioconvert fastq2fasta input.fastq output.fasta
    bioconvert fastq2fasta input.fq    output.fasta
    bioconvert fastq2fasta input.fq.gz output.fasta.gz
    bioconvert fastq2fasta input.fq.gz output.fasta.bz2

When there is no ambiguity, you can be implicit::

     bioconvert input.fastq output.fasta

The default method of conversion is used but you may use another one. Checkout the available methods with::

    bioconvert fastq2fasta --show-methods

For more help about a conversion, just type::

    bioconvert fastq2fasta --help

and more generally::

    bioconvert --help


You may also call **BioConvert** from a Python shell::

    # import a converter
    from bioconvert.fastq2fasta import FASTQ2FASTA

    # Instanciate with infile/outfile names
    convert = FASTQ2FASTA(infile, outfile)

    # the conversion itself:
    convert()


Available Converters
#######################


.. list-table:: Conversion table
    :widths: 20 40 40
    :header-rows: 1

    * - Converters
      - CI testing
      - Default method
    * - `abi2fasta <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.abi2fasta>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/abi2fasta.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/abi2fasta.yml
      - `BIOPYTHON <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `abi2fastq <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.abi2fastq>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/abi2fastq.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/abi2fastq.yml
      - `BIOPYTHON <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `abi2qual <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.abi2qual>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/abi2qual.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/abi2qual.yml
      - `BIOPYTHON <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `bam2bedgraph <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.bam2bedgraph>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2bedgraph.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2bedgraph.yml
      - `BEDTOOLS <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `bam2bigwig <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.bam2bigwig>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2bigwig.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2bigwig.yml
      - `DEEPTOOLS <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `bam2cov <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.bam2cov>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2cov.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2cov.yml
      - `BEDTOOLS <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `bam2cram <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.bam2cram>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2cram.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2cram.yml
      - `SAMTOOLS <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `bam2fasta <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.bam2fasta>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2fasta.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2fasta.yml
      - `SAMTOOLS <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `bam2fastq <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.bam2fastq>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2fastq.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2fastq.yml
      - `SAMTOOLS <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `bam2json <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.bam2json>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2json.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2json.yml
      - `BAMTOOLS <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `bam2sam <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.bam2sam>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2sam.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2sam.yml
      - `SAMBAMBA <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `bam2tsv <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.bam2tsv>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2tsv.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2tsv.yml
      - `SAMTOOLS <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `bam2wiggle <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.bam2wiggle>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2wiggle.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2wiggle.yml
      - `WIGGLETOOLS <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `bcf2vcf <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.bcf2vcf>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bcf2vcf.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bcf2vcf.yml
      - `BCFTOOLS <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `bcf2wiggle <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.bcf2wiggle>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bcf2wiggle.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bcf2wiggle.yml
      - `WIGGLETOOLS <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `bed2wiggle <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.bed2wiggle>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bed2wiggle.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bed2wiggle.yml
      - `WIGGLETOOLS <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `bedgraph2bigwig <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.bedgraph2bigwig>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bedgraph2bigwig.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bedgraph2bigwig.yml
      - `UCSC <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `bedgraph2cov <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.bedgraph2cov>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bedgraph2cov.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bedgraph2cov.yml
      - `BIOCONVERT <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `bedgraph2wiggle <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.bedgraph2wiggle>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bedgraph2wiggle.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bedgraph2wiggle.yml
      - `WIGGLETOOLS <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `bigbed2bed <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.bigbed2bed>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bigbed2bed.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bigbed2bed.yml
      - `DEEPTOOLS <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `bigbed2wiggle <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.bigbed2wiggle>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bigbed2wiggle.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bigbed2wiggle.yml
      - `WIGGLETOOLS <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `bigwig2bedgraph <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.bigwig2bedgraph>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bigwig2bedgraph.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bigwig2bedgraph.yml
      - `DEEPTOOLS <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `bigwig2wiggle <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.bigwig2wiggle>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bigwig2wiggle.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bigwig2wiggle.yml
      - `WIGGLETOOLS <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `bplink2plink <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.bplink2plink>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bplink2plink.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bplink2plink.yml
      - `PLINK <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `bplink2vcf <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.bplink2vcf>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bplink2vcf.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bplink2vcf.yml
      - `PLINK <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `bz22gz <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.bz22gz>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bz22gz.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bz22gz.yml
      - Unix commands
    * - `clustal2fasta <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.clustal2fasta>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/clustal2fasta.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/clustal2fasta.yml
      - `BIOPYTHON <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `clustal2nexus <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.clustal2nexus>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/clustal2nexus.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/clustal2nexus.yml
      - `GOALIGN <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `clustal2phylip <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.clustal2phylip>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/clustal2phylip.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/clustal2phylip.yml
      - `BIOPYTHON <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `clustal2stockholm <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.clustal2stockholm>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/clustal2stockholm.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/clustal2stockholm.yml
      - `BIOPYTHON <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `cram2bam <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.cram2bam>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/cram2bam.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/cram2bam.yml
      - `SAMTOOLS <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `cram2fasta <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.cram2fasta>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/cram2fasta.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/cram2fasta.yml
      - `SAMTOOLS <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `cram2fastq <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.cram2fastq>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/cram2fastq.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/cram2fastq.yml
      - `SAMTOOLS <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `cram2sam <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.cram2sam>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/cram2sam.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/cram2sam.yml
      - `SAMTOOLS <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `csv2tsv <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.csv2tsv>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/csv2tsv.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/csv2tsv.yml
      - `BIOCONVERT <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `csv2xls <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.csv2xls>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/csv2xls.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/csv2xls.yml
      - Pandas
    * - `dsrc2gz <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.dsrc2gz>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/dsrc2gz.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/dsrc2gz.yml
      - DSRC software
    * - `embl2fasta <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.embl2fasta>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/embl2fasta.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/embl2fasta.yml
      - `BIOPYTHON <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `embl2genbank <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.embl2genbank>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/embl2genbank.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/embl2genbank.yml
      - `BIOPYTHON <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `fasta2clustal <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.fasta2clustal>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2clustal.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2clustal.yml
      - `BIOPYTHON <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `fasta2faa <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.fasta2faa>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2faa.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2faa.yml
      - `BIOCONVERT <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `fasta2fasta_agp <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.fasta2fasta_agp>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2fasta_agp.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2fasta_agp.yml
      - `BIOCONVERT <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `fasta2fastq <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.fasta2fastq>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2fastq.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2fastq.yml
      - `PYSAM <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `fasta2genbank <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.fasta2genbank>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2genbank.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2genbank.yml
      - `BIOCONVERT <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `fasta2nexus <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.fasta2nexus>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2nexus.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2nexus.yml
      - `GOALIGN <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `fasta2phylip <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.fasta2phylip>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2phylip.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2phylip.yml
      - `BIOPYTHON <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `fasta2twobit <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.fasta2twobit>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2twobit.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2twobit.yml
      - `UCSC <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `fasta_qual2fastq <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.fasta_qual2fastq>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta_qual2fastq.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta_qual2fastq.yml
      - `PYSAM <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `fastq2fasta <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.fastq2fasta>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fastq2fasta.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fastq2fasta.yml
      -  `BIOCONVERT <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_  `available <_static/benchmark_fastq2fasta.png>`_
    * - `fastq2fasta_qual <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.fastq2fasta_qual>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fastq2fasta_qual.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fastq2fasta_qual.yml
      - `BIOCONVERT <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `fastq2qual <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.fastq2qual>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fastq2qual.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fastq2qual.yml
      - `READFQ <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `genbank2embl <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.genbank2embl>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/genbank2embl.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/genbank2embl.yml
      - `BIOPYTHON <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `genbank2fasta <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.genbank2fasta>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/genbank2fasta.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/genbank2fasta.yml
      - `BIOPYTHON <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `genbank2gff3 <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.genbank2gff3>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/genbank2gff3.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/genbank2gff3.yml
      - `BIOCODE <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `gfa2fasta <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.gfa2fasta>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/gfa2fasta.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/gfa2fasta.yml
      - `BIOCONVERT <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `gff22gff3 <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.gff22gff3>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/gff22gff3.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/gff22gff3.yml
      - `BIOCONVERT <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `gff32gff2 <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.gff32gff2>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/gff32gff2.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/gff32gff2.yml
      - `BIOCONVERT <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `gff32gtf <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.gff32gtf>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/gff32gtf.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/gff32gtf.yml
      - `BIOCONVERT <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `gz2bz2 <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.gz2bz2>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/gz2bz2.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/gz2bz2.yml
      - pigz/pbzip2 software
    * - `gz2dsrc <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.gz2dsrc>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/gz2dsrc.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/gz2dsrc.yml
      - DSRC software
    * - `json2yaml <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.json2yaml>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/json2yaml.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/json2yaml.yml
      - Python
    * - `maf2sam <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.maf2sam>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/maf2sam.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/maf2sam.yml
      - `BIOCONVERT <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `newick2nexus <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.newick2nexus>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/newick2nexus.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/newick2nexus.yml
      - `GOTREE <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `newick2phyloxml <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.newick2phyloxml>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/newick2phyloxml.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/newick2phyloxml.yml
      - `GOTREE <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `nexus2clustal <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.nexus2clustal>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/nexus2clustal.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/nexus2clustal.yml
      - `GOALIGN <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `nexus2fasta <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.nexus2fasta>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/nexus2fasta.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/nexus2fasta.yml
      - `BIOPYTHON <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `nexus2newick <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.nexus2newick>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/nexus2newick.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/nexus2newick.yml
      - `GOTREE <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `nexus2phylip <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.nexus2phylip>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/nexus2phylip.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/nexus2phylip.yml
      - `GOALIGN <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `nexus2phyloxml <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.nexus2phyloxml>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/nexus2phyloxml.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/nexus2phyloxml.yml
      - `GOTREE <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `ods2csv <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.ods2csv>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/ods2csv.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/ods2csv.yml
      - pyexcel library
    * - `pdb2faa <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.pdb2faa>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/pdb2faa.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/pdb2faa.yml
      - `BIOCONVERT <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `phylip2clustal <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.phylip2clustal>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2clustal.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2clustal.yml
      - `BIOPYTHON <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `phylip2fasta <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.phylip2fasta>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2fasta.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2fasta.yml
      - `BIOPYTHON <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `phylip2nexus <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.phylip2nexus>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2nexus.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2nexus.yml
      - `GOALIGN <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `phylip2stockholm <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.phylip2stockholm>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2stockholm.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2stockholm.yml
      - `BIOPYTHON <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `phylip2xmfa <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.phylip2xmfa>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2xmfa.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2xmfa.yml
      - `BIOPYTHON <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `phyloxml2newick <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.phyloxml2newick>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/phyloxml2newick.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/phyloxml2newick.yml
      - `GOTREE <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `phyloxml2nexus <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.phyloxml2nexus>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/phyloxml2nexus.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/phyloxml2nexus.yml
      - `GOTREE <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `plink2bplink <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.plink2bplink>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/plink2bplink.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/plink2bplink.yml
      - `PLINK <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `plink2vcf <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.plink2vcf>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/plink2vcf.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/plink2vcf.yml
      - `PLINK <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `sam2bam <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.sam2bam>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/sam2bam.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/sam2bam.yml
      - `SAMTOOLS <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `sam2cram <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.sam2cram>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/sam2cram.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/sam2cram.yml
      - `SAMTOOLS <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `sam2paf <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.sam2paf>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/sam2paf.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/sam2paf.yml
      - `BIOCONVERT <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `scf2fasta <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.scf2fasta>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/scf2fasta.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/scf2fasta.yml
      - `BIOCONVERT <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `scf2fastq <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.scf2fastq>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/scf2fastq.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/scf2fastq.yml
      - `BIOCONVERT <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `sra2fastq <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.sra2fastq>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/sra2fastq.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/sra2fastq.yml
      - `FASTQDUMP <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `stockholm2clustal <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.stockholm2clustal>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/stockholm2clustal.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/stockholm2clustal.yml
      - `BIOPYTHON <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `stockholm2phylip <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.stockholm2phylip>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/stockholm2phylip.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/stockholm2phylip.yml
      - `BIOPYTHON <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `tsv2csv <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.tsv2csv>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/tsv2csv.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/tsv2csv.yml
      - `BIOCONVERT <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `twobit2fasta <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.twobit2fasta>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/twobit2fasta.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/twobit2fasta.yml
      - `DEEPTOOLS <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `vcf2bcf <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.vcf2bcf>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/vcf2bcf.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/vcf2bcf.yml
      - `BCFTOOLS <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `vcf2bed <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.vcf2bed>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/vcf2bed.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/vcf2bed.yml
      - `BIOCONVERT <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `vcf2bplink <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.vcf2bplink>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/vcf2bplink.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/vcf2bplink.yml
      - `PLINK <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `vcf2plink <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.vcf2plink>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/vcf2plink.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/vcf2plink.yml
      - `PLINK <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `vcf2wiggle <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.vcf2wiggle>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/vcf2wiggle.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/vcf2wiggle.yml
      - `WIGGLETOOLS <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `wig2bed <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.wig2bed>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/wig2bed.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/wig2bed.yml
      - `BEDOPS <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `xls2csv <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.xls2csv>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/xls2csv.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/xls2csv.yml
      -
    * - `xlsx2csv <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.xlsx2csv>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/xlsx2csv.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/xlsx2csv.yml
      - Pandas library
    * - `xmfa2phylip <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.xmfa2phylip>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/xmfa2phylip.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/xmfa2phylip.yml
      - `BIOPYTHON <https://bioconvert.readthedocs.io/en/main/bibliography.html>`_
    * - `yaml2json <https://bioconvert.readthedocs.io/en/main/ref_converters.html#module-bioconvert.yaml2json>`_
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/yaml2json.yml/badge.svg
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/yaml2json.yml
      - Pandas library



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
1.1.0     * Implement ability to benchmark the CPU and memory usage (not just time)
            benchmark incorporates CPU/memory usage
          * add missing cython in requirements (and psutil for benchmarking)
1.0.0     * Fix bam2fastq for paired data that computed useless intermediate file
            https://github.com/bioconvert/bioconvert/issues/325
          * more realistic fastq simulator
          * pin openpyxl to <=3.0.10 to prevent regression error in v3.1.0
0.6.3     * add picard method in bam2sam
          * Fixed all CI workflows to use mamba
          * drop python3.7 support and add 3.10 support
          * update bedops test file to fit the latest bedops 2.4.41 version
          * revisit logging system
0.6.2     * added gff3 to gtf conversion. 
          * Added pdb to faa conversion
          * Added missing --reference argument to the cram2sam conversion
0.6.1     * output file can be in sub-directories allowing syntax such as 
            'bioconvert fastq2fasta test.fastq outputs/test.fasta
          * fix all CI actions
          * add more examples as notebooks in ./examples
          * add a Snakefile for the paper in ./doc/Snakefile_paper
0.6.0     * Fix bug in bam2sam (method sambamba)
          * Fix graph layout
          * add threading in fastq2fasta (seqkit method)
          * multibenchmark feature added
          * stable version used for web interface
0.5.2     * Update requirements and environment.yml and add a conda spec-file.txt file
0.5.1     * add genbank2gff3 requirement material in bioconvert.utils.biocode
0.5.0     * Add CI actions for all converters
          * remove sniffer (now in biosniff on pypi https://pypi.org/project/biosniff/)
          * A complete benchmarking suite (see doc/Snakefile_benchmark file and
            `benchmarking`)
          * documentation and tests for all converters
          * removed the validators (we assume intputs are correct)
0.4.X     * (aug 2019) added nexus2fasta, cram2fasta, fasta2faa ... ; 1-to-many and 
            many-to-one converters are now part of the API.
0.3.X       may 2019. new methods abi2qual, bigbed2bed, etc. added --threads option
0.2.X       aug 2018. abi2fastx, bioconvert_stats tool added
0.1.X       major refactoring to  have subcommands with implicit/explicit mode
========= ==============================================================================

