"""
Converter benchmarking
===========================

Converter have a default method.

Notem however, that several methods may be available.
Moreover, you may have a method that you want to compare
with the implemented one. To do so you will need to 
implement your method first. Then, simply use our benchmarking
framework as follows.

"""
#################################################
#
from bioconvert import Benchmark
from bioconvert import bioconvert_data
from bioconvert.bam2bed import BAM2BED

#####################################################
# Get the convert you wish to benchmark
input_file = bioconvert_data("test_measles.sorted.bam")
conv = BAM2BED(input_file, "test.bed")

#####################################################
# Get the Benchmark instance
bench = Benchmark(conv)
bench.plot()

# You can now see the different methods implemented in this
# converter and which one is the fastest.
