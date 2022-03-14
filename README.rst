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
    * - abi2fasta
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/abi2fasta.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/abi2fasta.yml
      - link IMG benchmarking abi2fasta
    * - abi2fastq
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/abi2fastq.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/abi2fastq.yml
      - link IMG benchmarking abi2fastq
    * - abi2qual
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/abi2qual.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/abi2qual.yml
      - link IMG benchmarking abi2qual
    * - bam2bedgraph
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2bedgraph.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2bedgraph.yml
      - link IMG benchmarking bam2bedgraph
    * - bam2bigwig
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2bigwig.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2bigwig.yml
      - link IMG benchmarking bam2bigwig
    * - bam2cov
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2cov.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2cov.yml
      - link IMG benchmarking bam2cov
    * - bam2cram
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2cram.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2cram.yml
      - link IMG benchmarking bam2cram
    * - bam2fasta
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2fasta.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2fasta.yml
      - link IMG benchmarking bam2fasta
    * - bam2fastq
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2fastq.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2fastq.yml
      - link IMG benchmarking bam2fastq
    * - bam2json
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2json.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2json.yml
      - link IMG benchmarking bam2json
    * - bam2sam
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2sam.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2sam.yml
      - link IMG benchmarking bam2sam
    * - bam2tsv
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2tsv.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2tsv.yml
      - link IMG benchmarking bam2tsv
    * - bam2wiggle
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bam2wiggle.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bam2wiggle.yml
      - link IMG benchmarking bam2wiggle
    * - bcf2vcf
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bcf2vcf.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bcf2vcf.yml
      - link IMG benchmarking bcf2vcf
    * - bcf2wiggle
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bcf2wiggle.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bcf2wiggle.yml
      - link IMG benchmarking bcf2wiggle
    * - bed2wiggle
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bed2wiggle.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bed2wiggle.yml
      - link IMG benchmarking bed2wiggle
    * - bedgraph2bigwig
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bedgraph2bigwig.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bedgraph2bigwig.yml
      - link IMG benchmarking bedgraph2bigwig
    * - bedgraph2cov
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bedgraph2cov.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bedgraph2cov.yml
      - link IMG benchmarking bedgraph2cov
    * - bedgraph2wiggle
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bedgraph2wiggle.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bedgraph2wiggle.yml
      - link IMG benchmarking bedgraph2wiggle
    * - bigbed2bed
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bigbed2bed.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bigbed2bed.yml
      - link IMG benchmarking bigbed2bed
    * - bigbed2wiggle
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bigbed2wiggle.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bigbed2wiggle.yml
      - link IMG benchmarking bigbed2wiggle
    * - bigwig2bedgraph
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bigwig2bedgraph.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bigwig2bedgraph.yml
      - link IMG benchmarking bigwig2bedgraph
    * - bigwig2wiggle
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bigwig2wiggle.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bigwig2wiggle.yml
      - link IMG benchmarking bigwig2wiggle
    * - bplink2plink
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bplink2plink.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bplink2plink.yml
      - link IMG benchmarking bplink2plink
    * - bplink2vcf
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bplink2vcf.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bplink2vcf.yml
      - link IMG benchmarking bplink2vcf
    * - bz22gz
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/bz22gz.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/bz22gz.yml
      - link IMG benchmarking bz22gz
    * - clustal2fasta
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/clustal2fasta.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/clustal2fasta.yml
      - link IMG benchmarking clustal2fasta
    * - clustal2nexus
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/clustal2nexus.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/clustal2nexus.yml
      - link IMG benchmarking clustal2nexus
    * - clustal2phylip
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/clustal2phylip.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/clustal2phylip.yml
      - link IMG benchmarking clustal2phylip
    * - clustal2stockholm
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/clustal2stockholm.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/clustal2stockholm.yml
      - link IMG benchmarking clustal2stockholm
    * - cram2bam
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/cram2bam.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/cram2bam.yml
      - link IMG benchmarking cram2bam
    * - cram2fasta
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/cram2fasta.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/cram2fasta.yml
      - link IMG benchmarking cram2fasta
    * - cram2fastq
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/cram2fastq.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/cram2fastq.yml
      - link IMG benchmarking cram2fastq
    * - cram2sam
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/cram2sam.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/cram2sam.yml
      - link IMG benchmarking cram2sam
    * - csv2tsv
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/csv2tsv.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/csv2tsv.yml
      - link IMG benchmarking csv2tsv
    * - csv2xls
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/csv2xls.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/csv2xls.yml
      - link IMG benchmarking csv2xls
    * - dsrc2gz
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/dsrc2gz.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/dsrc2gz.yml
      - link IMG benchmarking dsrc2gz
    * - embl2fasta
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/embl2fasta.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/embl2fasta.yml
      - link IMG benchmarking embl2fasta
    * - embl2genbank
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/embl2genbank.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/embl2genbank.yml
      - link IMG benchmarking embl2genbank
    * - fasta_qual2fastq
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta_qual2fastq.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta_qual2fastq.yml
      - link IMG benchmarking fasta_qual2fastq
    * - fasta2clustal
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2clustal.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2clustal.yml
      - link IMG benchmarking fasta2clustal
    * - fasta2faa
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2faa.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2faa.yml
      - link IMG benchmarking fasta2faa
    * - fasta2fasta_agp
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2fasta_agp.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2fasta_agp.yml
      - link IMG benchmarking fasta2fasta_agp
    * - fasta2fastq
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2fastq.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2fastq.yml
      - link IMG benchmarking fasta2fastq
    * - fasta2genbank
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2genbank.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2genbank.yml
      - link IMG benchmarking fasta2genbank
    * - fasta2nexus
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2nexus.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2nexus.yml
      - link IMG benchmarking fasta2nexus
    * - fasta2phylip
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2phylip.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2phylip.yml
      - link IMG benchmarking fasta2phylip
    * - fasta2twobit
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2twobit.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fasta2twobit.yml
      - link IMG benchmarking fasta2twobit
    * - fastq2fasta_qual
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fastq2fasta_qual.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fastq2fasta_qual.yml
      - link IMG benchmarking fastq2fasta_qual
    * - fastq2fasta
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fastq2fasta.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fastq2fasta.yml
      - link IMG benchmarking fastq2fasta
    * - fastq2qual
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/fastq2qual.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/fastq2qual.yml
      - link IMG benchmarking fastq2qual
    * - genbank2embl
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/genbank2embl.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/genbank2embl.yml
      - link IMG benchmarking genbank2embl
    * - genbank2fasta
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/genbank2fasta.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/genbank2fasta.yml
      - link IMG benchmarking genbank2fasta
    * - genbank2gff3
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/genbank2gff3.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/genbank2gff3.yml
      - link IMG benchmarking genbank2gff3
    * - gfa2fasta
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/gfa2fasta.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/gfa2fasta.yml
      - link IMG benchmarking gfa2fasta
    * - gff22gff3
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/gff22gff3.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/gff22gff3.yml
      - link IMG benchmarking gff22gff3
    * - gff32gff2
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/gff32gff2.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/gff32gff2.yml
      - link IMG benchmarking gff32gff2
    * - gz2bz2
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/gz2bz2.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/gz2bz2.yml
      - link IMG benchmarking gz2bz2
    * - gz2dsrc
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/gz2dsrc.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/gz2dsrc.yml
      - link IMG benchmarking gz2dsrc
    * - json2yaml
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/json2yaml.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/json2yaml.yml
      - link IMG benchmarking json2yaml
    * - maf2sam
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/maf2sam.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/maf2sam.yml
      - link IMG benchmarking maf2sam
    * - newick2nexus
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/newick2nexus.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/newick2nexus.yml
      - link IMG benchmarking newick2nexus
    * - newick2phyloxml
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/newick2phyloxml.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/newick2phyloxml.yml
      - link IMG benchmarking newick2phyloxml
    * - nexus2clustal
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/nexus2clustal.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/nexus2clustal.yml
      - link IMG benchmarking nexus2clustal
    * - nexusfasta
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/nexusfasta.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/nexusfasta.yml
      - link IMG benchmarking nexusfasta
    * - nexus2newick
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/nexus2newick.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/nexus2newick.yml
      - link IMG benchmarking nexus2newick
    * - nexus2phylip
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/nexus2phylip.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/nexus2phylip.yml
      - link IMG benchmarking nexus2phylip
    * - nexus2phyloxml
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/nexus2phyloxml.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/nexus2phyloxml.yml
      - link IMG benchmarking nexus2phyloxml
    * - ods2csv
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/ods2csv.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/ods2csv.yml
      - link IMG benchmarking ods2csv
    * - phylip2clustal
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2clustal.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2clustal.yml
      - link IMG benchmarking phylip2clustal
    * - phylip2fasta
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2fasta.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2fasta.yml
      - link IMG benchmarking phylip2fasta
    * - phylip2nexus
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2nexus.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2nexus.yml
      - link IMG benchmarking phylip2nexus
    * - phylip2stockholm
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2stockholm.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2stockholm.yml
      - link IMG benchmarking phylip2stockholm
    * - phylip2xmfa
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2xmfa.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/phylip2xmfa.yml
      - link IMG benchmarking phylip2xmfa
    * - phyloxml2newick
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/phyloxml2newick.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/phyloxml2newick.yml
      - link IMG benchmarking phyloxml2newick
    * - phyloxml2nexus
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/phyloxml2nexus.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/phyloxml2nexus.yml
      - link IMG benchmarking phyloxml2nexus
    * - plink2bplink
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/plink2bplink.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/plink2bplink.yml
      - link IMG benchmarking plink2bplink
    * - plink2vcf
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/plink2vcf.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/plink2vcf.yml
      - link IMG benchmarking plink2vcf
    * - sam2bam
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/sam2bam.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/sam2bam.yml
      - link IMG benchmarking sam2bam
    * - sam2cram
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/sam2cram.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/sam2cram.yml
      - link IMG benchmarking sam2cram
    * - sra2fastq
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/sra2fastq.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/sra2fastq.yml
      - link IMG benchmarking sra2fastq
    * - tsv2csv
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/tsv2csv.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/tsv2csv.yml
      - link IMG benchmarking tsv2csv
    * - twobit2fasta
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/twobit2fasta.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/twobit2fasta.yml
      - link IMG benchmarking twobit2fasta
    * - xls2csv
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/xls2csv.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/xls2csv.yml
      - link IMG benchmarking xls2csv    
    * - xlsx2csv
      - .. image:: https://github.com/bioconvert/bioconvert/actions/workflows/xlsx2csv.yml/badge.svg?branch=refactoring
            :target: https://github.com/bioconvert/bioconvert/actions/workflows/xlsx2csv.yml
      - link IMG benchmarking xlsx2csv
    


Contributors
############

Setting up and maintaining Bioconvert has been possible thanks to users and contributors. 
Thanks to all:

.. image:: https://contrib.rocks/image?repo=bioconvert/bioconvert
    :target: https://github.com/bioconvert/bioconvert/graphs/contributors


Changes
########

========= ====================================================================
Version   Description
========= ====================================================================
0.5.0     * Add CI actions for all converters
          * remove sniffer (now in biosniff)
========= ====================================================================

