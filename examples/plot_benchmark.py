# -*- coding: utf-8 -*-

###########################################################################
# Bioconvert is a project to facilitate the interconversion               #
# of life science data from one format to another.                        #
#                                                                         #
# Authors: see CONTRIBUTORS.rst                                           #
# Copyright © 2018  Institut Pasteur, Paris and CNRS.                     #
# See the COPYRIGHT file for details                                      #
#                                                                         #
# bioconvert is free software: you can redistribute it and/or modify      #
# it under the terms of the GNU General Public License as published by    #
# the Free Software Foundation, either version 3 of the License, or       #
# (at your option) any later version.                                     #
#                                                                         #
# bioconvert is distributed in the hope that it will be useful,           #
# but WITHOUT ANY WARRANTY; without even the implied warranty of          #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
# GNU General Public License for more details.                            #
#                                                                         #
# You should have received a copy of the GNU General Public License       #
# along with this program (COPYING file).                                 #
# If not, see <http://www.gnu.org/licenses/>.                             #
###########################################################################

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
from bioconvert.bam2cov import BAM2COV
from bioconvert.fastq2fasta import FASTQ2FASTA

#####################################################
# Get the convert you wish to benchmark
input_file = bioconvert_data("test_measles.sorted.bam")
conv = BAM2COV(input_file, "test.cov")
#input_file = bioconvert_data("test_fastq2fasta_v1.fastq")
#conv = FASTQ2FASTA(input_file, "test.fasta")


#####################################################
# Get the Benchmark instance
bench = Benchmark(conv)
bench.plot()

# You can now see the different methods implemented in this
# converter and which one is the fastest.
