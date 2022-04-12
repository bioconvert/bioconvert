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


.. list-table:: Conversion table
    :widths: 20 40 40
    :header-rows: 1

    * - Converters
      - CI testing
      - Benchmarking
    * - :mod:`~bioconvert.abi2fasta`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/abi2fasta.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/abi2fasta.yml
      - link IMG benchmarking abi2fasta
    * - :mod:`~bioconvert.abi2fastq`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/abi2fastq.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/abi2fastq.yml
      - link IMG benchmarking abi2fastq
    * - :mod:`~bioconvert.abi2qual`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/abi2qual.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/abi2qual.yml
      - link IMG benchmarking abi2qual
    * - :mod:`~bioconvert.bam2bedgraph`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2bedgraph.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2bedgraph.yml
      - link IMG benchmarking bam2bedgraph
    * - :mod:`~bioconvert.bam2bigwig`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2bigwig.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2bigwig.yml
      - link IMG benchmarking bam2bigwig
    * - :mod:`~bioconvert.bam2cov`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2cov.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2cov.yml
      - link IMG benchmarking bam2cov
    * - :mod:`~bioconvert.bam2cram`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2cram.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2cram.yml
      - link IMG benchmarking bam2cram
    * - :mod:`~bioconvert.bam2fasta`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2fasta.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2fasta.yml
      - link IMG benchmarking bam2fasta
    * - :mod:`~bioconvert.bam2fastq`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2fastq.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2fastq.yml
      - link IMG benchmarking bam2fastq
    * - :mod:`~bioconvert.bam2json`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2json.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2json.yml
      - link IMG benchmarking bam2json
    * - :mod:`~bioconvert.bam2sam`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2sam.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2sam.yml
      - link IMG benchmarking bam2sam
    * - :mod:`~bioconvert.bam2tsv`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2tsv.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2tsv.yml
      - link IMG benchmarking bam2tsv
    * - :mod:`~bioconvert.bam2wiggle`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2wiggle.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2wiggle.yml
      - link IMG benchmarking bam2wiggle
    * - :mod:`~bioconvert.bcf2vcf`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bcf2vcf.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bcf2vcf.yml
      - link IMG benchmarking bcf2vcf
    * - :mod:`~bioconvert.bcf2wiggle`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bcf2wiggle.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bcf2wiggle.yml
      - link IMG benchmarking bcf2wiggle
    * - :mod:`~bioconvert.bed2wiggle`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bed2wiggle.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bed2wiggle.yml
      - link IMG benchmarking bed2wiggle
    * - :mod:`~bioconvert.bedgraph2bigwig`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bedgraph2bigwig.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bedgraph2bigwig.yml
      - link IMG benchmarking bedgraph2bigwig
    * - :mod:`~bioconvert.bedgraph2cov`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bedgraph2cov.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bedgraph2cov.yml
      - link IMG benchmarking bedgraph2cov
    * - :mod:`~bioconvert.bedgraph2wiggle`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bedgraph2wiggle.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bedgraph2wiggle.yml
      - link IMG benchmarking bedgraph2wiggle
    * - :mod:`~bioconvert.bigbed2bed`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bigbed2bed.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bigbed2bed.yml
      - link IMG benchmarking bigbed2bed
    * - :mod:`~bioconvert.bigbed2wiggle`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bigbed2wiggle.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bigbed2wiggle.yml
      - link IMG benchmarking bigbed2wiggle
    * - :mod:`~bioconvert.bigwig2bedgraph`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bigwig2bedgraph.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bigwig2bedgraph.yml
      - link IMG benchmarking bigwig2bedgraph
    * - :mod:`~bioconvert.bigwig2wiggle`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bigwig2wiggle.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bigwig2wiggle.yml
      - link IMG benchmarking bigwig2wiggle
    * - :mod:`~bioconvert.bplink2plink`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bplink2plink.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bplink2plink.yml
      - link IMG benchmarking bplink2plink
    * - :mod:`~bioconvert.bplink2vcf`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bplink2vcf.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bplink2vcf.yml
      - link IMG benchmarking bplink2vcf
    * - :mod:`~bioconvert.bz22gz`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bz22gz.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bz22gz.yml
      - link IMG benchmarking bz22gz
    * - :mod:`~bioconvert.clustal2fasta`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/clustal2fasta.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/clustal2fasta.yml
      - link IMG benchmarking clustal2fasta
    * - :mod:`~bioconvert.clustal2nexus`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/clustal2nexus.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/clustal2nexus.yml
      - link IMG benchmarking clustal2nexus
    * - :mod:`~bioconvert.clustal2phylip`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/clustal2phylip.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/clustal2phylip.yml
      - link IMG benchmarking clustal2phylip`
    * - :mod:`~bioconvert.clustal2stockholm`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/clustal2stockholm.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/clustal2stockholm.yml
      - link IMG benchmarking clustal2stockholm
    * - :mod:`~bioconvert.cram2bam`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/cram2bam.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/cram2bam.yml
      - link IMG benchmarking cram2bam
    * - :mod:`~bioconvert.cram2fasta`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/cram2fasta.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/cram2fasta.yml
      - link IMG benchmarking cram2fasta
    * - :mod:`~bioconvert.cram2fastq`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/cram2fastq.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/cram2fastq.yml
      - link IMG benchmarking cram2fastq
    * - :mod:`~bioconvert.cram2sam`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/cram2sam.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/cram2sam.yml
      - link IMG benchmarking cram2sam
    * - :mod:`~bioconvert.csv2tsv`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/csv2tsv.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/csv2tsv.yml
      - link IMG benchmarking csv2tsv
    * - :mod:`~bioconvert.csv2xls`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/csv2xls.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/csv2xls.yml
      - link IMG benchmarking csv2xls
    * - :mod:`~bioconvert.dsrc2gz`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/dsrc2gz.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/dsrc2gz.yml
      - link IMG benchmarking dsrc2gz
    * - :mod:`~bioconvert.embl2fasta`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/embl2fasta.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/embl2fasta.yml
      - link IMG benchmarking embl2fasta
    * - :mod:`~bioconvert.embl2genbank`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/embl2genbank.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/embl2genbank.yml
      - link IMG benchmarking embl2genbank
    * - :mod:`~bioconvert.fasta_qual2fastq`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta_qual2fastq.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta_qual2fastq.yml
      - link IMG benchmarking fasta_qual2fastq
    * - :mod:`~bioconvert.fasta2clustal`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2clustal.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2clustal.yml
      - link IMG benchmarking fasta2clustal
    * - :mod:`~bioconvert.fasta2faa`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2faa.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2faa.yml
      - link IMG benchmarking fasta2faa
    * - :mod:`~bioconvert.fasta2fasta_agp`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2fasta_agp.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2fasta_agp.yml
      - link IMG benchmarking fasta2fasta_agp
    * - :mod:`~bioconvert.fasta2fastq`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2fastq.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2fastq.yml
      - link IMG benchmarking fasta2fastq
    * - :mod:`~bioconvert.fasta2genbank`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2genbank.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2genbank.yml
      - link IMG benchmarking fasta2genbank
    * - :mod:`~bioconvert.fasta2nexus`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2nexus.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2nexus.yml
      - link IMG benchmarking fasta2nexus
    * - :mod:`~bioconvert.fasta2phylip`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2phylip.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2phylip.yml
      - link IMG benchmarking fasta2phylip
    * - :mod:`~bioconvert.fasta2twobit`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2twobit.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2twobit.yml
      - link IMG benchmarking fasta2twobit
    * - :mod:`~bioconvert.fastq2fasta_qual`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fastq2fasta_qual.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fastq2fasta_qual.yml
      - link IMG benchmarking fastq2fasta_qual
    * - :mod:`~bioconvert.fastq2fasta`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fastq2fasta.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fastq2fasta.yml
      - link IMG benchmarking fastq2fasta
    * - :mod:`~bioconvert.fastq2qual`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fastq2qual.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fastq2qual.yml
      - link IMG benchmarking fastq2qual
    * - :mod:`~bioconvert.genbank2embl`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/genbank2embl.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/genbank2embl.yml
      - link IMG benchmarking genbank2embl
    * - :mod:`~bioconvert.genbank2fasta`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/genbank2fasta.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/genbank2fasta.yml
      - link IMG benchmarking genbank2fasta
    * - :mod:`~bioconvert.genbank2gff3`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/genbank2gff3.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/genbank2gff3.yml
      - link IMG benchmarking genbank2gff3
    * - :mod:`~bioconvert.gfa2fasta`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/gfa2fasta.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/gfa2fasta.yml
      - link IMG benchmarking gfa2fasta
    * - :mod:`~bioconvert.gff22gff3`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/gff22gff3.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/gff22gff3.yml
      - link IMG benchmarking gff22gff3
    * - :mod:`~bioconvert.gff32gff2`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/gff32gff2.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/gff32gff2.yml
      - link IMG benchmarking gff32gff2
    * - :mod:`~bioconvert.gz2bz2`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/gz2bz2.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/gz2bz2.yml
      - link IMG benchmarking gz2bz2
    * - :mod:`~bioconvert.gz2dsrc`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/gz2dsrc.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/gz2dsrc.yml
      - link IMG benchmarking gz2dsrc
    * - :mod:`~bioconvert.json2yaml`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/json2yaml.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/json2yaml.yml
      - link IMG benchmarking json2yaml
    * - :mod:`~bioconvert.maf2sam`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/maf2sam.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/maf2sam.yml
      - link IMG benchmarking maf2sam
    * - :mod:`~bioconvert.newick2nexus`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/newick2nexus.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/newick2nexus.yml
      - link IMG benchmarking newick2nexus
    * - :mod:`~bioconvert.newick2phyloxml`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/newick2phyloxml.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/newick2phyloxml.yml
      - link IMG benchmarking newick2phyloxml
    * - :mod:`~bioconvert.nexus2clustal`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/nexus2clustal.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/nexus2clustal.yml
      - link IMG benchmarking nexus2clustal
    * - :mod:`~bioconvert.nexus2fasta`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/nexus2fasta.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/nexus2fasta.yml
      - link IMG benchmarking nexus2fasta
    * - :mod:`~bioconvert.nexus2newick`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/nexus2newick.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/nexus2newick.yml
      - link IMG benchmarking nexus2newick
    * - :mod:`~bioconvert.nexus2phylip`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/nexus2phylip.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/nexus2phylip.yml
      - link IMG benchmarking nexus2phylip
    * - :mod:`~bioconvert.nexus2phyloxml`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/nexus2phyloxml.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/nexus2phyloxml.yml
      - link IMG benchmarking nexus2phyloxml
    * - :mod:`~bioconvert.ods2csv`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/ods2csv.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/ods2csv.yml
      - link IMG benchmarking ods2csv
    * - :mod:`~bioconvert.phylip2clustal`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2clustal.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2clustal.yml
      - link IMG benchmarking phylip2clustal
    * - :mod:`~bioconvert.phylip2fasta`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2fasta.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2fasta.yml
      - link IMG benchmarking phylip2fasta
    * - :mod:`~bioconvert.phylip2nexus`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2nexus.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2nexus.yml
      - link IMG benchmarking phylip2nexus
    * - :mod:`~bioconvert.phylip2stockholm`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2stockholm.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2stockholm.yml
      - link IMG benchmarking phylip2stockholm
    * - :mod:`~bioconvert.phylip2xmfa`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2xmfa.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2xmfa.yml
      - link IMG benchmarking phylip2xmfa
    * - :mod:`~bioconvert.phyloxml2newick`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/phyloxml2newick.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/phyloxml2newick.yml
      - link IMG benchmarking phyloxml2newick
    * - :mod:`~bioconvert.phyloxml2nexus`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/phyloxml2nexus.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/phyloxml2nexus.yml
      - link IMG benchmarking phyloxml2nexus
    * - :mod:`~bioconvert.plink2bplink`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/plink2bplink.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/plink2bplink.yml
      - link IMG benchmarking plink2bplink
    * - :mod:`~bioconvert.plink2vcf`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/plink2vcf.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/plink2vcf.yml
      - link IMG benchmarking plink2vcf
    * - :mod:`~bioconvert.sam2bam`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/sam2bam.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/sam2bam.yml
      - link IMG benchmarking sam2bam
    * - :mod:`~bioconvert.sam2cram`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/sam2cram.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/sam2cram.yml
      - link IMG benchmarking sam2cram
    * - :mod:`~bioconvert.sam2paf`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/sam2paf.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/sam2paf.yml
      - link IMG benchmarking sam2paf
    * - :mod:`~bioconvert.scf2fasta`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/scf2fasta.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/scf2fasta.yml
      - link IMG benchmarking scf2fasta
    * - :mod:`~bioconvert.scf2fastq`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/scf2fastq.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/scf2fastq.yml
      - link IMG benchmarking scf2fastq
    * - :mod:`~bioconvert.sra2fastq`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/sra2fastq.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/sra2fastq.yml
      - link IMG benchmarking sra2fastq
    * - :mod:`~bioconvert.stockholm2clustal`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/stockholm2clustal.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/stockholm2clustal.yml
      - link IMG benchmarking stockholm2clustal
    * - :mod:`~bioconvert.stockholm2phylip`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/stockholm2phylip.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/stockholm2phylip.yml
      - link IMG benchmarking stockholm2phylip
    * - :mod:`~bioconvert.tsv2csv`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/tsv2csv.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/tsv2csv.yml
      - link IMG benchmarking tsv2csv
    * - :mod:`~bioconvert.twobit2fasta`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/twobit2fasta.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/twobit2fasta.yml
      - link IMG benchmarking twobit2fasta
    * - :mod:`~bioconvert.vcf2bcf`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/vcf2bcf.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/vcf2bcf.yml
      - link IMG benchmarking vcf2bcf
    * - :mod:`~bioconvert.vcf2bed`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/vcf2bed.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/vcf2bed.yml
      - link IMG benchmarking vcf2bed
    * - :mod:`~bioconvert.vcf2bplink`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/vcf2bplink.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/vcf2bplink.yml
      - link IMG benchmarking vcf2bplink
    * - :mod:`~bioconvert.vcf2plink`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/vcf2plink.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/vcf2plink.yml
      - link IMG benchmarking vcf2plink
    * - :mod:`~bioconvert.vcf2wiggle`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/vcf2wiggle.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/vcf2wiggle.yml
      - link IMG benchmarking vcf2wiggle
    * - :mod:`~bioconvert.wig2bed`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/wig2bed.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/wig2bed.yml
      - link IMG benchmarking wig2bed
    * - :mod:`~bioconvert.xls2csv`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/xls2csv.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/xls2csv.yml
      - link IMG benchmarking xls2csv    
    * - :mod:`~bioconvert.xlsx2csv`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/xlsx2csv.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/xlsx2csv.yml
      - link IMG benchmarking xlsx2csv
    * - :mod:`~bioconvert.xmfa2phylip`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/xmfa2phylip.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/xmfa2phylip.yml
      - link IMG benchmarking xmfa2phylip
    * - :mod:`~bioconvert.yaml2json`
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/yaml2json.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/yaml2json.yml
      - link IMG benchmarking yaml2json
    


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
            :ref:`benchmarking`)
          * documentation and tests for all converters
          * removed the validators (we assume intputs are correct)
========= ==============================================================================

