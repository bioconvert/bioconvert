.. _benchmarking:

Benchmarking
============

Introduction
------------

Converters (e.g. :class:`~bioconvert.fastq2fasta.FASTQ2FASTA`) may have several
methods implemented. A developer may also want to compare his/her method with 
those available in bioconvert.

In order to help developers comparing their methods, we provide a benchmark
framework. One simply need to add his method inside the converter (see :ref:`developer_guide`) and use the method :meth:`~bioconvert.core.ConvBase.boxplot_benchmark`.

In practice, you would need a large test file. We do not provide such files
inside bioconverter.

In practice, you could use the following code to generate the boxplot:


.. plot::
    :include-source: 

    # Generate the dummy data, saving the results in a temporary file
    from easydev import TempFile
    from bioconvert.simulator.fastq import FastqSim

    infile = TempFile(suffix=".fastq")
    outfile = TempFile(suffix=".fasta")
    fs = FastqSim(infile.name)
    fs.nreads = 20000 # 1,000,000 by default
    fs.simulate()

    # Perform the benchmarking
    from bioconvert.fastq2fasta import FASTQ2FASTA
    c = FASTQ2FASTA(infile.name, outfile.name)
    c.boxplot_benchmark(to_exclude=["GATB"])

    infile.delete()
    outfile.delete()

Here, the boxplot_benchmark methods s called 5 times for each available method.

Be aware that the pure Python methods may be faster for small data set due to
a non-negligeable cost of using subprocess. 

You can try to provide a FastQ file with only 1 read and realise that it takes
about a second in all method. This is a incompressible delay so benchmarking needs 
large files to be meaningful !

If we use 1,000,000 reads instead of just 20,000, we would get different results
(which may change depending on your system and IO performance):

.. image:: benchmark.png

If we substract the subprocess cost from methods that use it, then we would get this result

.. image:: benchmark2.png

Here, what you see is that readfq and GATB methods (pure python) in this image
and the previous image give similar results. On the contrary, other methods (e.g., seqtk) 
are systematically shifted by about 0.5 to 1 second (cost of the subprocess).

So, your benchmarks should last  more than a few seconds.


Multiple benchmarking for more robustness
-------------------------------------------

With th previous method, even though you can decrease the error bars using more trials per method, we still suffer from
local computation or IO access that may bias the results. We provide a Snakefile here: :download:`Snakefile_benchmark`
that allows to run the previous benchmarking several times. So at the end you have a benchmark ... of benchmarks
somehow. We found it far more robust. Here is an example for the fastq2fasta case where each method was run 3 times and
in each case, 10 instances of conversion were performed. The orange vertical lines give the median and a final statement
indicates whther the final best method is significantly better than the others. 

.. image:: multi_benchmark.png


