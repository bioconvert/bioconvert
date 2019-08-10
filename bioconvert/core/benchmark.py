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
"""Tools for benchmarking"""
from collections import defaultdict
from itertools import chain
from pandas import np
from easydev import Timer, Progress

import colorlog
_log = colorlog.getLogger(__name__)


__all__ = ["Benchmark", "BenchmarkMulticonvert"]


def gmean(a, axis=0, dtype=None):
    # A copy/paste of scipy.stats.mstats.gmean function to 
    # avoid the scipy dependency
    if not isinstance(a, np.ndarray):
        # if not an ndarray object attempt to convert it
        log_a = np.log(np.array(a, dtype=dtype))
    elif dtype:
        # Must change the default dtype allowing array type
        if isinstance(a, np.ma.MaskedArray):
            log_a = np.log(np.ma.asarray(a, dtype=dtype))
        else:
            log_a = np.log(np.asarray(a, dtype=dtype))
    else:
        log_a = np.log(a)
    return np.exp(log_a.mean(axis=axis))


class Benchmark():
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
        self.include_dummy = False
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

        if self.include_dummy:
            methods += ['dummy']

        if self.to_include:
            methods = [x for x in methods if x in self.to_include]
        elif self.to_exclude:
            methods = [x for x in methods if x not in self.to_exclude]

        for method in methods:
            print("\nEvaluating method {}".format(method))
            times = []
            pb = Progress(self.N)
            for i in range(self.N):
                with Timer(times):
                    self.converter(method=method)
                pb.animate(i+1)
            results[method] = times
        self.results = results

    def plot(self, rerun=False, ylabel="Time (seconds)", rot_xticks=0, 
             boxplot_args={}):
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
        pylab.boxplot([data[x] for x in methods], **boxplot_args)
        # pylab.xticks([1+this for this in range(len(methods))], methods)
        if "vert" in boxplot_args and boxplot_args['vert'] is False:
            pylab.yticks(*zip(*enumerate(methods, start=1)), rotation=rot_xticks)
            pylab.xlabel(ylabel)
            #pylab.xlim([0, len(methods)+1])
        else:
            pylab.xticks(*zip(*enumerate(methods, start=1)), rotation=rot_xticks)
            pylab.ylabel(ylabel)
            pylab.xlim([0, len(methods)+1])

        pylab.grid(True)
        pylab.tight_layout()

        return data


class BenchmarkMulticonvert(Benchmark):
    """Convenient class to benchmark several methods for a series of converters

    ::

        from bioconvert.bam2cov import BAM2COV
        from bioconvert import BenchmarkMulticonvert
        c1 = BAM2COV(infile1, outfile1)
        c2 = BAM2COV(infile2, outfile2)
        b = BenchmarkMulticonvert([c1, c2], N=5)
        b.run_methods()
        b.plot()

    """
    def __init__(self, objs, **kwargs):
        """.. rubric:: constructor

        :param list objs: a list of converters
        :param int N: number of replicates
        :param list to_exclude: method to exclude from the benchmark

        """

        # Set self.converter to None
        super().__init__(None, **kwargs)
        self.converters = objs

    def run_methods(self):
        results = defaultdict(list)
        # We only test the methods common to all converters
        # (The intended use is with a list of converters all
        # having the same methods, but different input files)
        methods = set(self.converters[0].available_methods[:])  # a copy !
        for converter in self.converters[1:]:
            methods &= set(converter.available_methods[:])
        methods = sorted(methods)

        if self.include_dummy:
            methods += ['dummy']

        if self.to_include:
            methods = [x for x in methods if x in self.to_include]
        elif self.to_exclude:
            methods = [x for x in methods if x not in self.to_exclude]

        for method in methods:
            print("\nEvaluating method {}".format(method))
            # key: converter.infile
            # value: list of times
            times = defaultdict(list)
            pb = Progress(self.N)
            for i in range(self.N):
                for converter in self.converters:
                    with Timer(times[converter.infile]):
                        converter(method=method)
                pb.animate(i+1)
            # Normalize times so that each converter has comparable times
            mean_time = gmean(np.fromiter(chain(*times.values()), dtype=float))
            # median of ratios to geometric mean (c.f. DESeq normalization)
            scales = {conv: np.median(np.asarray(conv_times) / mean_time)
                      for conv, conv_times in times.items()}
            for (conv, conv_times) in times.items():
                scale = scales[conv]
                results[method].extend(
                    [conv_time / scale for conv_time in conv_times])
        self.results = results

    def plot(self, rerun=False, ylabel="Time (normalized seconds)"):
        super().plot(rerun, ylabel)
