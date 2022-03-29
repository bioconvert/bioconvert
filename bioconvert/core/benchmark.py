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
"""Tools for benchmarking"""
from asyncore import write
from collections import defaultdict
from itertools import chain

import numpy as np
from easydev import Timer
from tqdm import tqdm

import colorlog

### Import for creating the CSV file ###
import csv
from pathlib import Path
########################################

_log = colorlog.getLogger(__name__)


__all__ = ["Benchmark"]


class Benchmark:
    """Convenient class to benchmark several methods for a given converter

    ::

        c = BAM2COV(infile, outfile)
        b = Benchmark(c, N=5)
        b.run_methods()
        b.plot()

    """

    def __init__(self, obj, N=5, to_exclude=None, to_include=None):
        """.. rubric:: Constructor

        :param obj: can be an instance of a converter class or a class name
        :param int N: number of replicates
        :param list to_exclude: methods to exclude from the benchmark
        :param list to_include: methods to include ONLY

        Use one of *to_exclude* or *to_include*.
        If both are provided, only the *to_include* one is used.

        """
        if isinstance(obj, str):
            raise NotImplementedError

        self.converter = obj
        self.N = N
        self.results = None
        if to_exclude is None:
            self.to_exclude = []
        else:
            self.to_exclude = to_exclude
        if to_include is None:
            self.to_include = []
        else:
            self.to_include = to_include

    def run_methods(self):
        """Runs the benchmarks, and stores the timings in *self.results*."""
        results = {}
        methods = self.converter.available_methods[:]  # a copy !

        if self.to_include:
            methods = [x for x in methods if x in self.to_include]
        elif self.to_exclude:
            methods = [x for x in methods if x not in self.to_exclude]

        for method in methods:
            times = []
            for i in tqdm(range(self.N), desc="Evaluating method {}".format(method)):
                with Timer(times):
                    # Need to get all extra arguments for specify method e.g Bam2BIGWIG.uscs method
                    kwargs = {"method": method}
                    for k, v in self.converter.others.items():
                        kwargs[k] = v
                    self.converter(**kwargs)
            results[method] = times
        self.results = results

    def plot(self, rerun=False, ylabel="Time (seconds)", rot_xticks=0, boxplot_args={}):
        """Plots the benchmark results, running the benchmarks
        if needed or if *rerun* is True.

        :param rot_xlabel: rotation of the xticks function
        :param boxplot_args: dictionary with any of the pylab.boxplot arguments
        :return: dataframe with all results
        """
        import pylab

        if self.results is None or rerun is True:
            self.run_methods()

        # an alias.
        data = self.results.copy()

        methods = sorted(data, key=lambda x: pylab.mean(data[x]))

        ########### Data recovery for the creation of the CSV file ###########
        
        converter_name = str(self.converter).split('.')[1]
        benchmark = 1
        fileObj = Path(f"{converter_name}.csv")
        
        if(fileObj.exists()):
            with open(f"{converter_name}.csv",'a',newline='') as f:  # Opening the CSV file by adding a line
                with open(f"{converter_name}.csv",'r',newline='') as read:  # Opening the CSV file for reading
                    my_reader = csv.reader(read)
                    for line in my_reader: # line is a list of column values
                        if(line[0] != "Benchmark"):
                            if(benchmark < int(line[0])): # line[0] is the value of the 1st column of the considered row
                                benchmark = int(line[0])
                    benchmark += 1
                write=csv.writer(f)  
                for key, value in data.items():
                    for i in value:
                        l = [benchmark, key, i]
                        # l.append(benchmark)
                        # l.append(key)
                        # l.append(i)
                        write.writerow(l)
        else : 
            with open(f"{converter_name}.csv",'w',newline='') as f:  # Opening the CSV file for writing
                write=csv.writer(f)  
                write.writerow(["Benchmark", "Method", "Value"])
                for key, value in data.items():
                    for i in value:
                        l = [benchmark, key, i]
                        # l.append(benchmark)
                        # l.append(key)
                        # l.append(i)
                        write.writerow(l)

        ######################################################################

        pylab.boxplot([data[x] for x in methods], **boxplot_args)
        # pylab.xticks([1+this for this in range(len(methods))], methods)
        if "vert" in boxplot_args and boxplot_args["vert"] is False:
            pylab.yticks(*zip(*enumerate(methods, start=1)), rotation=rot_xticks)
            pylab.xlabel(ylabel)
            # pylab.xlim([0, len(methods)+1])
        else:
            pylab.xticks(*zip(*enumerate(methods, start=1)), rotation=rot_xticks)
            pylab.ylabel(ylabel)
            pylab.xlim([0, len(methods) + 1])

        pylab.grid(True)
        pylab.tight_layout()

        return data
