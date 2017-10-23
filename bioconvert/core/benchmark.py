import time
from easydev import Timer, Progress
import pylab

import colorlog
_log = colorlog.getLogger(__name__)


__all__ = ["Benchmark"]


class Benchmark():
    """Convenient class to benchmark several methods for a given converter

    ::

        c = Bam2Bed(infile, outfile)
        b = Benchmark(c, N=5)
        b.run_methods()
        b.plot()

    """
    def __init__(self, obj, N=5):
        """

        :param obj: can be an instance of a converter class or a class name

        """
        if isinstance(obj, str):
            raise NotImplementedError

        self.converter = obj
        self.N = N
        self.results = None

    def run_methods(self):
        results = {}
        for method in self.converter.available_methods:
            _log.info("Evaluating method %s" % method)
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
        pylab.xticks([1+this for this in range(self.N)], methods)
        pylab.grid(True)
        pylab.ylabel("Time (seconds)")
        pylab.xlim([0, len(methods)+1])


