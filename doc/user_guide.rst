User Guide
==========

.. contents::

Quick Start
-----------

Most of the time, Bioconvert simply requires the input and output filenames.
If there is no ambiguity, the extension are used to infer the type of conversion
you wish to perform.

For instance, to convert a FASTQ to a FASTA file, use this type of command::

    bioconvert test.fastq test.fasta

If the converter *fastq* to *fasta*² exists in **Bioconvert**, it will work out of
the box. In order to get a list of all possible conversions, just type::

    bioconvert --help

To obtain more specific help about a converter that you found in the list::

    bioconvert fastq2fasta --help

.. note:: All converters are named as <input_extension>2<output_extension>


Explicit conversion
-------------------

Sometimes, Bioconvert won't be able to know what you want solely based on
the input and ouput extensions. So, you may need to be explicit and use a
subcommand. For instance to use the converter *fastq2fasta*, type::


    bioconvert fastq2fasta  input.fastq output.fasta

The rationale behind the subcommand choice is manyfold. First, you may have dedicated help
for a given conversion, which may be different from one conversion to the other::

    bioconvert fastq2fasta --help

Second, the extensions of your input and output may be non-standard or different
from the choice made by the **bioconvert** developers. So, using the subcommand you can do::

    bioconvert fastq2fasta  input.fq output.fas

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

Input files may be compressed. For instance, most FASTQ are compressed in GZ
format. Compression are handled in some converters. Basically, most of the
humand-readable files handle compression. For instance, all those commands
should work and can be used to compress output files, or handle input compressed
files::

    bioconvert test.fastq.gz test.fasta
    bioconvert test.fastq.gz test.fasta.gz
    bioconvert test.fastq.gz test.fasta.bz2

Note that you can also decompress and compress into another compression keeping
without doing any conversion (note the fastq extension in both input and output
files)::

    bioconvert test.fastq.gz test.fastq.dsrc



Parallelization
---------------

In Bioconvert, if the input contains a wildcard such as ``*`` or ``?`` characters, then, input filenames are treated separately and converted sequentially::

    bioconvert fastq2fasta "*.fastq"

Note, however, that the files are processed sequentially one by one. So, we may
want to parallelise the computation.

Iteration with unix commands
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can use a bash script under unix to run Bioconvert on a set of files. For
instance the following script takes all files with the *.fastq* extension and
convert them to fasta:


.. literalinclude:: script.sh
    :language: shell

Note, however, that this is still a sequential computation. Yet, you may now
change it slightly to run the commands on a cluster. For instance, on a SLURM
scheduler, you can use:

.. literalinclude:: script2.sh
    :language: shell

Snakemake option
~~~~~~~~~~~~~~~~

If you have lots of files to convert, a snakemake pipeline is available in the `Sequana <https://github.com/sequa    na/sequana>`_ project and can be installed using **pip install sequana_bioconvert**. It also installs bioconvert with an ap    ptainer image that contains all dependencies for you.



Here is another way of running your jobs in parallel using a 
simple Snakefile (`snakemake <https://snakemake.readthedocs.io/en/stable/>`_)
that can be run easily either locally or on a cluster.

You can download the following file :download:`Snakefile` 

.. literalinclude:: Snakefile
    :language: python
    :lines: 12-28

and execute it locally as follows (assuming you have 4 CPUs)::

    snakemake -s Snakefile --cores 4

or on a cluster::

    snakemake -s Snakefile --cluster "--mem=1000 -j 10"


