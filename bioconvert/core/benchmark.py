import time
from easydev import Timer
import pylab


class Benchmark():
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
            times = []
            for i in range(self.N):
                with Timer(times):
                    self.call(method=method)
            results[method] = times
        self.results = results

    def plot(self, rerun=False):
        if self.results is None or rerun is True:
            self.run_methods()
        # an alias
        data = self.results

        methods = sorted(data, key=lambda x: mean(data[x]))
        pylab.boxplot([data[x] for x in methods])
        pylab.xticks(range(self.N), methods)
