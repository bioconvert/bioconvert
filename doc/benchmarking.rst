Benchmarking
============


Converters (e.g. :class:`Fastq2Fasta`) may have several methods. A developr may
also want to compare his/her method with those available in bioconvert. 

In order to help developers comparing their methods, we provide a benchmark
framework. One simply need to add his method inside the converter (see :ref:`developer_guide`).


In practice, you would need a large test file. We do not provide such files
inside bioconverter. However, we provide some utilities to generate dummy data. 
For the moment, only FastQ dummy generator is provided. 

.. plot::
    :include-source: 

    # Generate the dummy data, saving the results in a temporary file
    from easydev import TempFile
    from bioconvert.simulator.fastq import FastqSim

    infile = TempFile(suffix=".fastq")
    outfile = TempFile(suffix=".fasta")
    fs = FastqSim(infile.name)
    fs.nreads = 2000 # 1,000,000 by default
    fs.simulate()

    # Perfrm the benchmarking
    from bioconvert.fastq2fasta import Fastq2Fasta
    c = Fastq2Fasta(infile.name, outfile.name)
    c.boxplot_benchmark(to_exclude=["gatb"])

    infile.delete()
    outfile.delete()

Here, the boxplot_benchmark methods call 5 times each method available. 

It looks like GATB and biopython are the fastest. In reality, this
benchmark is biased. Indeed, GATB and biopython being in pure Python, there is
no starting cost. While the seqtk uses a subprocess, which currently has a 1
second cost. Therefore, it is essntial to use large data sets for testing. 

We have a dummy method, which does nothing else that calling a subprocess that
does ... nothing. You could use it that way::

    c.boxplot_benchmark(include_dummy=True)

and see that indeed this takes 1 second like seqtk. So we need more data.


If we use 10,000,000 reads instead of just 2,000, we would get this results
(which may change significantly depeding on your system and IO performance):


.. image:: benchmark.png

Here, the times start at about 15-20 seconds that is much larger than 1 second.
So we are more confident in  this benchmark.



