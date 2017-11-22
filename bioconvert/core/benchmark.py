# -*- coding: utf-8 -*-
#
#  This file is part of Bioconvert software
#
#  Copyright (c) 2017 - Bioconvert Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/biokit/bioconvert
#  documentation: http://bioconvert.readthedocs.io
#
##############################################################################
"""Tools for benchmarking"""
from easydev import Timer, Progress
import pylab
from itertools import cycle

import colorlog
_log = colorlog.getLogger(__name__)


__all__ = ["Benchmark", "Benchmark_multiconvert"]


class Benchmark():
    """Convenient class to benchmark several methods for a given converter

    ::

        c = Bam2Bed(infile, outfile)
        b = Benchmark(c, N=5)
        b.run_methods()
        b.plot()

    """
    def __init__(self, obj, N=5, to_exclude=[], to_include=[]):
        """.. rubric:: constructor

        :param obj: can be an instance of a converter class or a class name
        :param int N: number of replicates
        :param list to_exclude: methods to exclude from the benchmark
        :param list to_include: methods to include ONLY

        Use one of to_exclude or to_include. If both are provided, only, to_include
        is used.

        """
        if isinstance(obj, str):
            raise NotImplementedError

        self.converter = obj
        self.N = N
        self.results = None
        self.include_dummy = False
        self.to_exclude = to_exclude
        self.to_include = to_include


    def run_methods(self):
        results = {}
        methods = self.converter.available_methods[:]  # a copy !

        if self.include_dummy:
            methods += ['dummy']

        if self.to_include:
            methods = [x for x in methods if x in self.to_include]
        elif self.to_exclude:
            methods = [x for x in methods if x not in self.to_exclude]

        for method in methods:
            print("\nEvaluating method %s" % method)
            times = []
            pb = Progress(self.N)
            for i in range(self.N):
                with Timer(times):
                    self.converter(method=method)
                pb.animate(i+1)
            results[method] = times
        self.results = results

    def plot(self, rerun=False):
        if self.results is None or rerun is True:
            self.run_methods()
        # an alias
        data = self.results

        methods = sorted(data, key=lambda x: pylab.mean(data[x]))
        pylab.boxplot([data[x] for x in methods])
        pylab.xticks([1+this for this in range(len(methods))], methods)
        pylab.grid(True)
        pylab.ylabel("Time (seconds)")
        pylab.xlim([0, len(methods)+1])


class Benchmark_multiconvert():
    """Convenient class to benchmark several methods for a series of converters

    ::

        c1 = Bam2Bed(infile1, outfile1)
        c2 = Bam2Bed(infile2, outfile2)
        b = Benchmark_multiconvert([c1, c2], N=5)
        b.run_methods()
        b.plot()

    """
    def __init__(self, objs, N=5, to_exclude=[], weights=None):
        """.. rubric:: constructor

        :param list objs: a list of converters
        :param int N: number of replicates
        :param list to_exclude: method to exclude from the benchmark
        :param list weights: weights to apply to timings of the converters

        *weights* should either be None or a list of numbers
        of the same length as *objs*.

        """

        self.converters = objs
        self.N = N
        self.results = None
        self.include_dummy = False
        self.to_exclude = to_exclude
        if weights is None:
            # Give equal weight to all converters
            self.weights = [1 for _ in objs]
        else:
            self.weights = weights

    def run_methods(self):
        results = {}
        # We only test the methods common to all converters
        # (The intended use is with a list of converters all
        # having the same methods, but different input files)
        methods = set(self.converters[0].available_methods[:])  # a copy !
        for converter in self.converters[1:]:
            methods &= set(converter.available_methods[:])
        methods = sorted(methods)
        if self.include_dummy:
            methods += ['dummy']

        for method in methods:
            if method in self.to_exclude:
                continue
            print("\nEvaluating method %s" % method)
            times = []
            pb = Progress(self.N)
            for i in range(self.N):
                for converter in self.converters:
                    with Timer(times):
                        converter(method=method)
                pb.animate(i+1)
            results[method] = [real_t * wt for (real_t, wt) in zip(
                times, cycle(self.weights))]
        self.results = results

    def plot(self, rerun=False):
        if self.results is None or rerun is True:
            self.run_methods()
        # an alias
        data = self.results

        methods = sorted(data, key=lambda x: pylab.mean(data[x]))
        pylab.boxplot([data[x] for x in methods])
        pylab.xticks([1+this for this in range(len(methods))], methods)
        pylab.grid(True)
        pylab.ylabel("Time (weighted seconds)")
        pylab.xlim([0, len(methods)+1])
