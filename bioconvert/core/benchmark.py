###########################################################################
# Bioconvert is a project to facilitate the interconversion               #
# of life science data from one format to another.                        #
#                                                                         #
# Copyright © 2018-2022 Institut Pasteur, Paris and CNRS.                 #
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
#                                                                         #
# Repository: https://github.com/bioconvert/bioconvert                    #
# Documentation: http://bioconvert.readthedocs.io                         #
###########################################################################
"""Tools for benchmarking"""

import statistics

import colorlog
import matplotlib.pyplot as plt
import pandas as pd
import pylab
import statsmodels.api
import statsmodels.formula.api
import statsmodels.stats.multitest
from easydev import Timer
from tqdm import tqdm

_log = colorlog.getLogger(__name__)


__all__ = ["Benchmark", "plot_multi_benchmark_max"]


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

        if self.results is None or rerun is True:
            self.run_methods()

        # an alias.
        data = self.results.copy()

        methods = sorted(data, key=lambda x: pylab.mean(data[x]))
        pylab.boxplot([data[x] for x in methods], **boxplot_args)
        # pylab.xticks([1+this for this in range(len(methods))], methods)
        if "vert" in boxplot_args and boxplot_args["vert"] is False:
            pylab.yticks(*zip(*enumerate(methods, start=1)), rotation=rot_xticks)
            pylab.xlabel(ylabel)
        else:
            pylab.xticks(*zip(*enumerate(methods, start=1)), rotation=rot_xticks)
            pylab.ylabel(ylabel)
            pylab.xlim([0, len(methods) + 1])

        pylab.grid(True)
        pylab.tight_layout()

        return data


def plot_multi_benchmark_max(path_json, output_filename="multi_benchmark.png", min_ylim=0):
    """Plotting function for the Snakefile_benchmark to be found in the doc

    The json file looks like::


        {
          "awk":{
            "0":0.777020216,
            "1":0.9638044834,
            "2":1.7623617649,
            "3":0.8348755836
          },
          "seqtk":{
            "0":1.0024843216,
            "1":0.6313509941,
            "2":1.4048073292,
            "3":1.0554351807
          },
          "Benchmark":{
            "0":1,
            "1":1,
            "2":2,
            "3":2
          }
        }

    Number of benchmark is infered from field 'Benchmark'.

    """
    # open and read JSON file
    df = pd.read_json(path_json)

    # how many runs per method ?
    N = len(df["Benchmark"]) / len(df["Benchmark"].unique())
    N = len(df["Benchmark"]) / N

    # Retrieving method names
    methods = list(df)

    # Removed the entry from the list that matches benchmark
    methods.remove("Benchmark")

    # We rotate the JSON object in relation to the benchmark number to be able to group them by methods
    df2 = df.pivot(columns="Benchmark")

    # Display of the boxplots with font at 7 and the grid is removed to make it easier to read
    df2.boxplot(fontsize=7, grid=False)

    # Initializing the x-axis title placement list
    l = []
    # Initialization of the variable which will be used to separate the different methods
    sep = 0
    # Initialization of the variable that will be used to save the lowest median
    median_best = 0
    # Initialization of the variable that will be used to save the name of the method with the lowest median
    best_method = None
    for i in methods:
        if sep != 0:
            # Creation of a separation between methods
            plt.axvline(sep + 0.5, ls="--", color="red")
        # Calculation of the median for a method
        median = statistics.median(df[i])
        if median_best == 0 or median_best >= median:
            # If the calculated median is lower than the recorded one, the old median and the old method are placed by the new one
            median_best = median
            best_method = i
        # We plot the median of each method
        plt.hlines(y=median, xmin=0.5 + sep, xmax=N + 0.5 + sep, color="orange")
        l.append(sep + (N / 2 + 0.5))
        sep += N

    # The name of each method is displayed on the x-axis
    plt.xticks(l, methods)
    # Creation of the path variable which will give the title of the output PNG image
    path = path_json.split("/")[1]
    _, max_ylim = plt.ylim()
    plt.ylim([min_ylim, max_ylim])
    # Backup of the benchmark of the different conversion in the form of a PNG image
    plt.savefig(output_filename, dpi=200)

    ############################## T-TEST ##############################
    # We recover the different times of the best method
    value_best_method = df[best_method]
    # Initialization of the dictionary which will save the results of the t-test of each method
    t_test = dict()
    for i in methods:
        if i != best_method:
            # We recover the different times of the method
            value_method = df[i]
            # Application of the t-test between the best method and all the other methods and saving these results in the dictionnary t_test
            comp = statsmodels.stats.weightstats.CompareMeans.from_data(value_best_method, value_method)
            (T_stats, P_value, degrees_f) = comp.ttest_ind()
            T_dict = {"t-stats": T_stats}
            P_dict = {"p-value": P_value}
            D_dict = {"Degree of freedom": degrees_f}
            t_test[i] = (T_dict, P_dict, D_dict)

    ############################## MULTITEST ##############################
    # Initialization of the list which will store all the p-values of the previous t-tests
    list_p_value = []
    # Initialization of the list which will store all the methods other than the best method
    list_method = []

    for i in t_test:
        # Retrieval of the p-values and the name of the different methods
        list_p_value.append(t_test[i][1]["p-value"])
        list_method.append(i)

    # Application of the multitest on the different conversion methods using the bonferroni method
    (
        areSignificant,
        correctedPvalues,
        _,
        _,
    ) = statsmodels.stats.multitest.multipletests(list_p_value, alpha=0.05, method="bonferroni")

    for i in range(len(list_method)):
        print(
            f"- By comparing the {list_method[i]} method with the best one ({best_method}), we check H0: {areSignificant[i]} with corrected P-value: {correctedPvalues[i]}"
        )
