import os

import pytest

from bioconvert import Benchmark
from bioconvert.bam2cov import BAM2COV
from easydev import TempFile

from .. import test_dir


@pytest.mark.skipif(
    "DISPLAY" not in os.environ, reason="no DISPLAY available, will fail otherwise"
)
def test_benchmark():
    input_file = f"{test_dir}/data/bam/test_measles.sorted.bam"
    with TempFile(suffix=".cov") as fout:
        conv = BAM2COV(input_file, fout.name)
        bench = Benchmark(conv)
        bench.run_methods()
        bench.plot()

        try:
            Benchmark("BAM2COV")
            assert False
        except NotImplementedError:
            assert True
        except:
            assert False
