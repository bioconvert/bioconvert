import os

import pytest

from bioconvert import Benchmark
from bioconvert import bioconvert_data
from bioconvert.bam2bed import BAM2BED
from easydev import TempFile


@pytest.mark.skipif("DISPLAY" not in os.environ, reason="no DISPLAY available, will fail otherwise")
@pytest.mark.skipif("TRAVIS_OS_NAME" == "osx", reason="unknown failure on travis april 2019")
def test_benchmark():
    input_file = bioconvert_data("test_measles.sorted.bam")
    with TempFile(suffix=".bed") as fout:
        conv = BAM2BED(input_file, fout.name)
        bench = Benchmark(conv)
        bench.run_methods()
        bench.plot()

        try:
            Benchmark("BAM2BED")
            assert False
        except NotImplementedError:
            assert True
        except:
            assert False
