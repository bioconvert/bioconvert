User Guide
============

.. contents::

Quick Start
-----------

If you want to convert a format to another and you knwo the extensions of the
output format, just try bioconvert naively::

    bioconvert test.fastq test.fasta

If the converter fastq to fasta exists in **Bioconvert**, it will work out of
the box. In order to get a list of all possible conversions, just type::

    bioconvert

or for more details::

    bioconvert --help

To obtain more specific help about a converter that you found in the list::

    bioconvert fastq2fasta --help

.. note:: All converters are named as <INPUT_EXTENSION>2<OUTPUT_EXTENSION>


Explicit conversion
--------------------


You can use **bioconvert** from a developer point of view, or as an end-user.
Here we describe the standalone application that is::

    bioconvert

You can obtain help by using::

    bioconvert --help

To convert a format into another you need to provide the name of the conversion.
For instance, to convert a fastq file into a fasta, you need to use the sub
command **fastq2fasta**::

    bioconvert fastq2fasta  input.fastq output.fasta

The rationale behind the subcommand choice is manyfold. First, you may have dedicated help
for a given conversion, which may be different from one conversion to the other::

    bioconvert fastq2fasta --help

Second, the extensions of your input and output may be non-standard or different
from the choice made by the **bioconvert** developers. So, using the subcommand you can do::

    bioconvert fastq2fasta  input.fq output.fa

where the extensions can actually be whatever you want.

If you do not provide the output file, it will be created based on the input
filename by replacing the extension automatically. So this command::

    bioconvert fastq2fasta input.fq

generates an output file called *input.fasta*. Note that it will be placed in
the same directory as the input file, not locally. So::

    bioconvert fastq2fasta ~/test/input.fq

will create the *input.fasta* file in the ~/test directory.

If an output file exists, it will not be overwritten. If you want to do so, use
the --force argument::

    bioconvert fastq2fasta  input.fq output.fa --force

Implicit conversion
-------------------

If the extensions match the conversion name, you can perform implicit
conversion::

    bioconvert input.fastq output.fasta

Internally, a format  may be registered with several extensions. For instance
the extensions possible for a FastA file are ``fasta`` and ``fa`` so you can
also write::

    bioconvert input.fastq output.fa

Compression
-----------

.. todo:: this section will be coming soon


Parallelization
--------------------


Some converters can use several threads, but if you have hundreds of files and
wish to use several CPUs, or run bioconvert on a cluster, we provide here below
a simple Snakefile (snakemake) that can be run easily as follows.

You can download the following file :download:`Snakefile` 

.. literalinclude:: Snakefile
    :language: python
    :lines: 12-28

and execute it locally as follows (assuming you have 4 CPUs)::

    snakemake -s Snakefile --cores 4

or on a cluster::

    snakemake -s Snakefile --cluster "--mem=1000 -j 10"








